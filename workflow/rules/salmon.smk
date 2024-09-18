rule salmon_index:
    input:
        genome=get_genome(),
        transcriptome=get_full_transcriptome(),
    output:
        decoys="results/salmon_index/decoys.txt",
        gentrome="results/salmon_index/gentrome.fa",
        outdir=directory("results/salmon_index"),
    log:
        "logs/salmon/salmon_index.log",
    conda:
        "../envs/salmon.yaml"
    shell:
        "grep '^>' {input.genome} |"
        " cut -d ' ' -f 1 |"
        " sed 's/>//g' > {output.decoys}"
        " && cat {input.transcriptome} {input.genome} > {output.gentrome}"
        " && salmon index"
        " --threads {resources.threads}"
        " --decoys {output.decoys}"
        " --transcripts {output.gentrome}"
        " --index {output.outdir}"
        " &> {log}"


rule salmon_quant:
    input:
        fq1="results/fastp/{sample}_trimmed_R1.fastq.gz",
        fq2="results/fastp/{sample}_trimmed_R2.fastq.gz",
        index="results/salmon_index",
    output:
        quant="results/salmon_quant/{sample}/quant.sf",
    log:
        "logs/salmon/quant/{sample}.log",
    params:
        outdir="results/salmon_quant/{sample}",
    conda:
        "../envs/salmon.yaml"
    shell:
        "salmon quant"
        " -i {input.index}"
        " --libType A"
        " -1 {input.fq1}"
        " -2 {input.fq2}"
        " -o {params.outdir}"
        " --threads {resources.threads}"
        " &> {log}"


rule salmon_tximport:
    input:
        gtf=get_gtf(),
        quant=expand("results/salmon_quant/{sample}/quant.sf", sample=SAMPLES),
    output:
        tx_tpm="results/salmon_quant/salmon_tpm_gene.tsv",
        gene_tpm="results/salmon_quant/salmon_tpm_transcript.tsv",
        tx_count="results/salmon_quant/salmon_counts_gene.tsv",
        gene_count="results/salmon_quant/salmon_counts_transcript.tsv",
    log:
        "logs/salmon/salmon_tximport.log",
    params:
        quantdir="results/salmon_quant",
    conda:
        "../envs/deseq2.yaml"
    script:
        "../scripts/merge_salmon_quant.R"
