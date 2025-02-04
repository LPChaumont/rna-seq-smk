rule gtf2bed:
    input:
        get_ref_file(config["ref"]["gtf"]),
    output:
        "resources/gtf2bed.bed",
    log:
        "logs/gtf2bed.log",
    conda:
        "../envs/rseqc.yaml"
    shell:
        "gxf2bed --input {input} --output {output} &> {log}"


rule infer_experiment:
    input:
        bed="resources/gtf2bed.bed",
        bam="results/star_align/{sample}/{sample}_Aligned.sortedByCoord.out.bam",
    output:
        "results/rseqc/infer_experiment/{sample}.txt",
    log:
        "logs/rseqc/infer_experiment/{sample}.log",
    conda:
        "../envs/rseqc.yaml"
    shell:
        "infer_experiment.py -r {input.bed} -i {input.bam} 1> {output} 2> {log}"


rule bam_stats:
    input:
        bam="results/star_align/{sample}/{sample}_Aligned.sortedByCoord.out.bam",
    output:
        "results/rseqc/bam_stats/{sample}.txt",
    log:
        "logs/rseqc/bam_stats/{sample}.log",
    conda:
        "../envs/rseqc.yaml"
    shell:
        "bam_stat.py -i {input.bam} 1> {output} 2> {log}"
