localrules:
    install_pairadise,


rule install_pairadise:
    output:
        outdir=directory("resources/PAIRADISE"),
        script="resources/install_r_deps.R",
    params:
        link_pairadise="https://github.com/Xinglab/PAIRADISE.git",
        link_r_deps="https://raw.githubusercontent.com/Xinglab/rmats-turbo/master/install_r_deps.R",
    conda:
        "../envs/rmats.yaml"
    shell:
        "cd resources/"
        " && git clone {params.link_pairadise}"
        " && wget {params.link_r_deps}"
        " && Rscript install_r_deps.R paired"


rule rmats_config_prep:
    input:
        bam="results/star_align/{sample}/Aligned.sortedByCoord.out.bam",
    output:
        config="results/rmats/bam_config_prep/{sample}.txt",
    params:
        outdir="results/rmats/bam_config_prep",
    shell:
        "mkdir -p {params.outdir}"
        " && echo {input.bam} > {output.config}"


rule rmats_config_post:
    input:
        bam=expand(
            "results/star_align/{sample}/Aligned.sortedByCoord.out.bam",
            sample=SAMPLES,
        ),
    output:
        config="results/rmats/bam_config_post/all_samples.txt",
    params:
        outdir="results/rmats/bam_config_post",
    shell:
        "mkdir -p {params.outdir}"
        " && echo {input.bam} | tr ' ' ',' > {output.config}"


rule rmats_prep:
    input:
        gtf=get_annotation(),
        config="results/rmats/bam_config_prep/{sample}.txt",
    output:
        prep_dir=directory("results/rmats/prep/{sample}"),
        flag=touch("results/rmats/prep/cp_with_prefix/{sample}.txt"),
    params:
        post_dir="results/rmats/post",
        read_length=config["read_length"],
        post_tmp_dir="results/rmats/post/tmp",
        prefix="{sample}_",
    log:
        "logs/rmats_prep/{sample}.log",
    conda:
        "../envs/rmats.yaml"
    shell:
        "rmats.py"
        " --gtf {input.gtf}"
        " --tmp {output.prep_dir}"
        " --od {params.post_dir}"
        " --readLength {params.read_length}"
        " --b1 {input.config}"
        " -t paired"
        " --nthread {resources.threads}"
        " --task prep"
        " --variable-read-length"
        " &> {log}"
        " && python $CONDA_PREFIX/rMATS/cp_with_prefix.py"
        " {params.prefix}"
        " {params.post_tmp_dir}"
        " {output.prep_dir}/*.rmats"
        " &>> {log}"


rule rmats_post:
    input:
        flag=expand("results/rmats/prep/cp_with_prefix/{sample}.txt", sample=SAMPLES),
        gtf=get_annotation(),
        config="results/rmats/bam_config_post/all_samples.txt",
    output:
        summary="results/rmats/post/summary.txt",
    log:
        "logs/rmats_post/rmats_post.log",
    params:
        tmp="results/rmats/prep",
        outdir="results/rmats/post",
        read_length=config["read_length"],
    conda:
        "../envs/rmats.yaml"
    shell:
        "rmats.py"
        " --gtf {input.gtf}"
        " --tmp {params.tmp}"
        " --od {params.outdir}"
        " --readLength {params.read_length}"
        " --b1 {input.config}"
        " -t paired"
        " --nthread {resources.threads}"
        " --task post"
        " --variable-read-length"
        " --statoff"
        " &> {log}"


rule rmats_stat:
    input:
        config_post="results/rmats/bam_config_post/all_samples.txt",
        post="results/rmats/post/summary.txt",
    output:
        summary="results/rmats/stat/{contrast}/summary.txt",
    log:
        "logs/rmats_stat/{contrast}.log",
    params:
        stat_dir="results/rmats/stat/{contrast}",
        stat_tmp_dir="results/rmats/stat/{contrast}/tmp",
        post_dir="results/rmats/post",
        group_1=lambda wc: get_rmats_group_indices(wc)[0],
        group_2=lambda wc: get_rmats_group_indices(wc)[1],
        paired_stats=use_rmats_pairadise(),
    conda:
        "../envs/rmats.yaml"
    shell:
        "python $CONDA_PREFIX/rMATS/rMATS_P/prepare_stat_inputs.py"
        " --new-output-dir {params.stat_dir}"
        " --old-output-dir {params.post_dir}"
        " --group-1-indices {params.group_1}"
        " --group-2-indices {params.group_2}"
        " &> {log}"
        " && rmats.py"
        " --nthread {resources.threads}"
        " --od {params.stat_dir}"
        " --tmp {params.stat_tmp_dir}"
        " --task stat"
        " {params.paired_stats}"
        " &>> {log}"


rule rmats_filtering:
    input:
        "results/rmats/stat/{contrast}/summary.txt",
    output:
        filtered="results/rmats/stat/{contrast}/{filter}_{event}.MATS.{junction}.txt",
    params:
        outdir="results/rmats/stat/{contrast}",
        # Relative path based on the output directory
        script="../../../../workflow/scripts/rmats_filtering.py",
    shell:
        "cd {params.outdir} &&"
        " python {params.script} {wildcards.event}.MATS.{wildcards.junction}.txt"
