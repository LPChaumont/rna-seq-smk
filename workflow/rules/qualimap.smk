rule qualimap_bamqc:
    input:
        bam="results/star_align/{sample}/Aligned.sortedByCoord.out.bam",
    output:
        outdir=directory("results/qualimap_bamqc/{sample}"),
    log:
        "logs/qualimap/bamqc/{sample}.log",
    conda:
        "../envs/qualimap.yaml"
    shell:
        "qualimap bamqc"
        " -bam {input.bam}"
        " -outdir {output.outdir}"
        " -outformat HTML"
        " -nt {resources.threads}"
        " &> {log}"


rule qualimap_multi_bamqc_input:
    output:
        "results/qualimap_bamqc/multi_bamqc_input.txt",
    log:
        "logs/qualimap/multi_bamqc/input.log",
    run:
        sys.stderr = open(log[0], "w")
        with open(output[0], "w") as f:
            for sample in SAMPLES:
                _sample = samples.loc[(sample), "sample"]
                bamqc = f"results/qualimap_bamqc/{sample}"
                condition = samples.loc[(sample), "condition"]

                line = "\t".join([_sample, bamqc, condition]) + "\n"
                f.write(line)


rule qualimap_multi_bamqc:
    input:
        bamqc=expand("results/qualimap_bamqc/{sample}", sample=SAMPLES),
        config="results/qualimap_bamqc/multi_bamqc_input.txt",
    output:
        directory("results/qualimap_multi_bamqc"),
    log:
        "logs/qualimap/multi_bamqc/multi_bamqc.log",
    conda:
        "../envs/qualimap.yaml"
    shell:
        "qualimap multi-bamqc"
        " --data {input.config}"
        " --outdir {output}"
        " &> {log}"


rule sam_sort_by_name:
    input:
        "results/star_align/{sample}/Aligned.sortedByCoord.out.bam",
    output:
        "results/star_align/{sample}/Aligned.sortedByName.out.bam",
    log:
        "logs/samtools/sort_by_name/{sample}.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools sort -n -O bam"
        " -@ {resources.threads}"
        " {input}"
        " 1> {output} 2> {log}"


# TODO fix qualimap rnaseq multiqc
rule qualimap_rnaseq:
    input:
        bam="results/star_align/{sample}/Aligned.sortedByName.out.bam",
        gtf=get_gtf(),
    output:
        outdir=directory("results/qualimap_rnaseq/{sample}"),
    log:
        "logs/qualimap/rnaseq/{sample}.log",
    conda:
        "../envs/qualimap.yaml"
    shell:
        "qualimap rnaseq"
        " -bam {input.bam}"
        " -gtf {input.gtf}"
        " -outdir {output}"
        " --paired"
        " --sorted"
        " --java-mem-size={resources.mem_gb}G"
        " &> {log}"
