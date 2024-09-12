rule sam_sort_idx:
    input:
        "results/star_align/{sample}/Aligned.sortedByCoord.out.bam",
    output:
        bam="results/sam_sort_idx/star_aln_{sample}.bam",
        bai="results/sam_sort_idx/star_aln_{sample}.bam.bai",
    log:
        "logs/sam_sort_idx/{sample}.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools sort {input} -o {output} -@ {resources.threads}"
        " && samtools index {output}"


rule sam_stats:
    input:
        "results/sam_sort_idx/star_aln_{sample}.bam",
    output:
        "results/sam_stats/{sample}.txt",
    log:
        "logs/sam_stats/{sample}.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools stats {input}"
        " -@ {resources.threads}"
        " 1> {output}"
        " 2> {log}"


rule sam_flagstat:
    input:
        "results/sam_sort_idx/star_aln_{sample}.bam",
    output:
        "results/sam_flagstat/{sample}.tsv",
    log:
        "logs/sam_flagstat/{sample}.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools flagstat {input}"
        " -@ {resources.threads}"
        " --output-fmt tsv"
        " 1> {output}"
        " 2> {log}"


rule sam_idxstats:
    input:
        "results/sam_sort_idx/star_aln_{sample}.bam.bai",
    output:
        "results/sam_idxstats/{sample}.tsv",
    log:
        "logs/sam_idxstats/{sample}.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools idxstats {input}"
        " -@ {resources.threads}"
        " 1> {output}"
        " 2> {log}"


rule multiqc_sam:
    input:
        expand("results/sam_stats/{sample}.txt", sample=SAMPLES)
        + expand("results/sam_flagstat/{sample}.tsv", sample=SAMPLES)
        + expand("results/sam_idxstats/{sample}.tsv", sample=SAMPLES),
    output:
        outdir=directory("results/multiqc/samtools"),
        report="results/multiqc/samtools/multiqc_samtools_report.html",
    params:
        indir="results/sam_stats results/sam_flagstat results/sam_idxstats",
        filename="multiqc_samtools_report.html",
    log:
        "logs/multiqc/samtools.log",
    conda:
        "../envs/multiqc.yaml"
    shell:
        "multiqc"
        " {params.indir}"
        " --outdir {output.outdir}"
        " --filename {params.filename}"
        " --force"
        " &> {log}"
