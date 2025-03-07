localrules:
    install_pairadise,


rule install_pairadise:
    output:
        outdir=directory("PAIRADISE"),
    params:
        link="https://github.com/Xinglab/PAIRADISE.git",
        install_dir="./PAIRADISE/pairadise/src/pairadise_model/",
    log:
        "logs/rmats/install_pairadise.log",
    conda:
        "../envs/rmats.yaml"
    shell:
        """
        git clone {params.link} 2> {log}
        Rscript -e "install.packages('{params.install_dir}', repos=NULL)" 2>> {log}
        """


rule rmats_config_prep:
    input:
        bam="results/star_align/{sample}/{sample}_Aligned.sortedByCoord.out.bam",
    output:
        config="results/rmats/bam_config_prep/{sample}.txt",
    log:
        "logs/rmats/rmats_config_prep/{sample}.log",
    conda:
        "../envs/rmats.yaml"
    shell:
        "echo {input.bam} > {output.config} 2> {log}"


rule rmats_config_post:
    input:
        bam=expand(
            "results/star_align/{sample}/{sample}_Aligned.sortedByCoord.out.bam",
            sample=SAMPLES,
        ),
    output:
        config="results/rmats/bam_config_post/all_samples.txt",
    log:
        "logs/rmats/rmats_config_post.log",
    conda:
        "../envs/rmats.yaml"
    shell:
        "echo {input.bam} | tr ' ' ',' > {output.config} 2> {log}"


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
        "rmats.py"
        " --gtf {input.gtf}"
        " --tmp {output.prep_dir}"
        " --od {params.post_dir}"
        " --readLength {params.read_length}"
        " --b1 {input.config}"
        " -t paired"
        " --libType {params.lib_type}"
        " --nthread {resources.threads}"
        " --task prep"
        " --variable-read-length"
        " --allow-clipping"
        " &> {log}"
        " && mkdir -p {params.post_dir} {params.post_tmp_dir}"
        " && python $CONDA_PREFIX/rMATS/cp_with_prefix.py"
        " {params.prefix}"
        " {params.post_tmp_dir}"
        " {output.prep_dir}/*.rmats"
        " &>> {log}"


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
        " --allow-clipping"
        " --statoff"
        " &> {log}"


rule rmats_prep_stat:
    input:
        config_post="results/rmats/bam_config_post/all_samples.txt",
        post="results/rmats/post/summary.txt",
    output:
        flag=touch("results/rmats/stat/{contrast}/prep_stat.txt"),
    log:
        "logs/rmats/stat/prep_{contrast}.log",
    params:
        stat_dir=lambda w, output: os.path.dirname(output[0]),
        post_dir=lambda w, input: os.path.dirname(input[1]),
        group_1=lambda w, input: get_rmats_group_indices(input[0], w)[0],
        group_2=lambda w, input: get_rmats_group_indices(input[0], w)[1],
    conda:
        "../envs/rmats.yaml"
    shell:
        "python $CONDA_PREFIX/rMATS/rMATS_P/prepare_stat_inputs.py"
        " --new-output-dir {params.stat_dir}"
        " --old-output-dir {params.post_dir}"
        " --group-1-indices {params.group_1}"
        " --group-2-indices {params.group_2}"
        " &> {log}"


rule rmats_stat:
    input:
        flag="results/rmats/stat/{contrast}/prep_stat.txt",
    output:
        summary="results/rmats/stat/{contrast}/summary.txt",
    log:
        "logs/rmats/stat/{contrast}.log",
    params:
        stat_dir=lambda w, output: os.path.dirname(output[0]),
        paired_stats=use_rmats_pairadise(),
    conda:
        "../envs/rmats.yaml"
    shell:
        "rmats.py"
        " --nthread {resources.threads}"
        " --od {params.stat_dir}"
        " --tmp {params.stat_dir}/tmp"
        " --task stat"
        " {params.paired_stats}"
        " &>> {log}"


rule rmats_filtering:
    input:
        "results/rmats/stat/{contrast}/summary.txt",
    output:
        "results/rmats/stat/{contrast}/summary_filtered.tsv",
    params:
        rmats_dir=lambda w, input: os.path.dirname(input[0]),
    log:
        "logs/rmats/filtering/{contrast}.log",
    conda:
        "../envs/rmats.yaml"
    shell:
        "python workflow/scripts/rmats_filtering.py"
        " --rmats-dir {params.rmats_dir}"
        " --read-cov 10"
        " --sig-fdr 0.05"
        " --sig-delta-psi 0.1"
        " &> {log}"
