localrules:
    get_collapse_gtf_script,


rule get_collapse_gtf_script:
    output:
        "resources/collapse_annotation.py",
    params:
        "https://raw.githubusercontent.com/broadinstitute/gtex-pipeline/master/gene_model/collapse_annotation.py",
    shell:
        "wget -P resources/ {params}"


# combining all isoforms of a gene into a single transcript
rule rnaseqc_collapse_gtf:
    input:
        script="resources/collapse_annotation.py",
        gtf=get_gtf(),
    output:
        "results/rnasesqc/collapse_rnaseqc.gtf",
    conda:
        "../envs/rnaseqc.yaml"
    shell:
        "python {input.script} {input.gtf} {output}"


rule rnaseqc:
    input:
        gtf="resources/collapse_rnaseqc.gtf",
        bam="results/star_align/{sample}/Aligned.sortedByCoord.out.bam",
    output:
        outdir=directory("results/rnasesqc/{sample}"),
    log:
        "logs/rnaseqc/{sample}.log",
    conda:
        "../envs/rnaseqc.yaml"
    shell:
        "rnaseqc"
        " {input.gtf}"
        " {input.bam}"
        " {output.outdir}"
        " &> {log}"
