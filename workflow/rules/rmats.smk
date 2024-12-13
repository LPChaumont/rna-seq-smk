localrules:
    install_pairadise,


rule install_pairadise:
    output:
        outdir=directory("resources/PAIRADISE"),
        script="resources/install_r_deps.R",
    params:
        link_pairadise="https://github.com/Xinglab/PAIRADISE.git",
        link_r_deps="https://raw.githubusercontent.com/Xinglab/rmats-turbo/master/install_r_deps.R",
    log:
        "logs/rmats/install_pairadise.log",
    conda:
        "../envs/rmats.yaml"
    shell:
        """
        cd resources/
        && git clone {params.link_pairadise} 2> {log}
        && wget {params.link_r_deps} 2>> {log}
        && Rscript install_r_deps.R paired 2>> {log}
        """


rule rmats_config_prep:
    input:
        bam="results/star_align/{sample}/Aligned.sortedByCoord.out.bam",
    output:
        config="results/rmats/bam_config_prep/{sample}.txt",
    params:
        outdir=lambda w, output: os.path.dirname(output[0]),
    log:
        "logs/rmats/rmats_config_prep/{sample}.log",
    conda:
        "../envs/rmats.yaml"
    shell:
        """mkdir -p {params.outdir}
        && echo {input.bam} > {output.config} 2> {log}
        """


rule rmats_config_post:
    input:
        bam=expand(
            "results/star_align/{sample}/Aligned.sortedByCoord.out.bam",
            sample=SAMPLES,
        ),
    output:
        config="results/rmats/bam_config_post/all_samples.txt",
    params:
        outdir=lambda w, output: os.path.dirname(output[0]),
    log:
        "logs/rmats/rmats_config_post.log",
    conda:
        "../envs/rmats.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        && echo {input.bam} | tr ' ' ',' > {output.config} 2> {log}
        """


rule rmats_prep:
    input:
        gtf=get_ref_file(config["ref"]["gtf"]),
        config="results/rmats/bam_config_prep/{sample}.txt",
    output:
        prep_dir=directory("results/rmats/prep/{sample}"),
        flag=touch("results/rmats/prep/cp_with_prefix/{sample}.txt"),
    params:
        post_dir="results/rmats/post",
        read_length=config["read_length"],
        post_tmp_dir="results/rmats/post/tmp",
        prefix="{sample}_",
        lib_type=config["rmats"]["lib_type"],
    log:
        "logs/rmats/prep/{sample}.log",
    conda:
        "../envs/rmats.yaml"
    shell:
        """
        rmats.py
        --gtf {input.gtf}
        --tmp {output.prep_dir}
        --od {params.post_dir}
        --readLength {params.read_length}
        --b1 {input.config}
        -t paired
        --libType {params.lib_type}
        --nthread {resources.threads}
        --task prep
        --variable-read-length
        &> {log}
        && mkdir -p {params.post_dir} {params.post_tmp_dir}
        && python $CONDA_PREFIX/rMATS/cp_with_prefix.py
        {params.prefix}
        {params.post_tmp_dir}
        {output.prep_dir}/*.rmats
        &>> {log}
        """


rule rmats_post:
    input:
        flag=expand("results/rmats/prep/cp_with_prefix/{sample}.txt", sample=SAMPLES),
        gtf=get_ref_file(config["ref"]["gtf"]),
        config="results/rmats/bam_config_post/all_samples.txt",
    output:
        summary="results/rmats/post/summary.txt",
    log:
        "logs/rmats/rmats_post.log",
    params:
        tmp="results/rmats/post/tmp",
        outdir=lambda w, output: os.path.dirname(output[0]),
        read_length=config["read_length"],
    conda:
        "../envs/rmats.yaml"
    shell:
        """
        rmats.py"
        --gtf {input.gtf}
        --tmp {params.tmp}
        --od {params.outdir}
        --readLength {params.read_length}
        --b1 {input.config}
        -t paired
        --nthread {resources.threads}
        --task post
        --variable-read-length
        --statoff
        &> {log}
        """


rule rmats_stat:
    input:
        config_post="results/rmats/bam_config_post/all_samples.txt",
        post="results/rmats/post/summary.txt",
    output:
        summary="results/rmats/stat/{contrast}/summary.txt",
    log:
        "logs/rmats/stat/{contrast}.log",
    params:
        stat_dir=lambda w, output: os.path.dirname(output[0]),
        stat_tmp_dir="results/rmats/stat/{contrast}/tmp",
        post_dir=lambda w, input: os.path.dirname(input.post),
        group_1=lambda w: get_rmats_group_indices(w)[0],
        group_2=lambda w: get_rmats_group_indices(w)[1],
        paired_stats=use_rmats_pairadise(),
    conda:
        "../envs/rmats.yaml"
    shell:
        """
        python $CONDA_PREFIX/rMATS/rMATS_P/prepare_stat_inputs.py
        --new-output-dir {params.stat_dir}
        --old-output-dir {params.post_dir}
        --group-1-indices {params.group_1}
        --group-2-indices {params.group_2}
        &> {log}
        && rmats.py
        --nthread {resources.threads}
        --od {params.stat_dir}
        --tmp {params.stat_tmp_dir}
        --task stat
        {params.paired_stats}
        &>> {log}
        """


rule rmats_filtering:
    input:
        "results/rmats/stat/{contrast}/summary.txt",
    output:
        filtered="results/rmats/stat/{contrast}/{filter}_{event}.MATS.{junction}.txt",
    params:
        outdir=lambda w, output: os.path.dirname(output[0]),
        # Relative path based on the output directory
        script="../../../../workflow/scripts/rmats_filtering.py",
    log:
        "logs/rmats/filtering/{contrast}_{filter}_{event}_{junction}.log",
    conda:
        "../envs/rmats.yaml"
    shell:
        """
        cd {params.outdir} &&
        python {params.script} {wildcards.event}.MATS.{wildcards.junction}.txt
        """
