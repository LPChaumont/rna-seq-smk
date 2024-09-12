rule qualimap_bamqc:
    input:
        bam="results/star_align/{sample}/Aligned.sortedByCoord.out.bam",
    output:
        outdir=directory("results/qualimap_bamqc/{sample}"),
    log:
        "logs/qualimap_bamqc/{sample}.log",
    conda:
        "../envs/qualimap.yaml"
    shell:
        "qualimap bamqc"
        " -bam {input.bam}"
        " -outdir {output.outdir}"
        " -outformat HTML"
        " -nt {resources.threads}"
        " &> {log}"


rule qualimap_rnaseq:
    input:
        bam="results/star_align/{sample}/Aligned.sortedByCoord.out.bam",
        gtf=get_annotation(),
    output:
        outdir=directory("results/qualimap_rnaseq/{sample}"),
    log:
        "logs/qualimap_rnaseq/{sample}.log",
    conda:
        "../envs/qualimap.yaml"
    shell:
        "qualimap rnaseq"
        " -bam {input.bam}"
        " -gtf {input.gtf}"
        " -outdir {output}"
        " -pe"
        " &> {log}"


rule qualimap_multiqc:
    input:
        data=expand(
            "results/qualimap_rnaseq/{sample}",
            sample=SAMPLES,
            dir=["qualimap_bamqc", "qualimap_rnaseq"],
        ),
    output:
        outdir=directory("results/multiqc/qualimap"),
        report="results/multiqc/qualimap/multiqc_qualimap_report.html",
    params:
        filename="multiqc_qualimap_report.html",
        indir=["qualimap_bamqc", "qualimap_rnaseq"],
    log:
        "logs/multiqc/qualimap.log",
    conda:
        "../envs/multiqc.yaml"
    shell:
        "multiqc"
        " {params.indir}"
        " --outdir {output.outdir}"
        " --filename {params.filename}"
        " --force"
        " &> {log}"
