rule salmon_decoys:
    input:
        genome=get_ref_file(config["ref"]["dna_fasta"]),
        transcriptome=get_transcriptome(),
    output:
        decoys="results/salmon/decoys/decoys.txt",
        gentrome="results/salmon/decoys/gentrome.fa",
    log:
        "logs/salmon/decoys.log",
    conda:
        "../envs/salmon.yaml"
    shell:
        "grep '^>' {input.genome} |"
        " cut -d ' ' -f 1 |"
        " sed 's/>//g' 1> {output.decoys} 2> {log}"
        " cat {input.transcriptome} {input.genome} 1> {output.gentrome} 2>> {log}"


rule salmon_index:
    input:
        decoys="results/salmon/decoys/decoys.txt",
        gentrome="results/salmon/decoys/gentrome.fa",
    output:
        outdir=directory("results/salmon/index/"),
    log:
        "logs/salmon/index.log",
    conda:
        "../envs/salmon.yaml"
    shell:
        "salmon --no-version-check index"
        " --index {output.outdir}"
        " --threads {resources.threads}"
        " --decoys {input.decoys}"
        " --transcripts {input.gentrome}"
        " &> {log}"


rule salmon_quant:
    input:
        fq1="results/fastp/{sample}_trimmed_R1.fastq.gz",
        fq2="results/fastp/{sample}_trimmed_R2.fastq.gz",
        index="results/salmon/index",
    output:
        quant="results/salmon/quant/{sample}/quant.sf",
    log:
        "logs/salmon/quant/{sample}.log",
    params:
        outdir=lambda w, output: os.path.dirname(output[0]),
        lib_type=config["salmon"]["lib_type"],
    conda:
        "../envs/salmon.yaml"
    shell:
        "salmon --no-version-check quant"
        " -i {input.index}"
        " --libType {params.lib_type}"
        " -1 {input.fq1}"
        " -2 {input.fq2}"
        " -o {params.outdir}"
        " --threads {resources.threads}"
        " &> {log}"


rule tx2gene:
    input:
        get_ref_file(config["ref"]["gtf"]),
    output:
        "resources/tx2gene.tsv",
    log:
        "logs/salmon/tximport/tx2gene.log",
    conda:
        "../envs/tximport.yaml"
    script:
        "../scripts/tx2gene.R"


rule salmon_tximport:
    input:
        tx2gene="resources/tx2gene.tsv",
        quant=expand("results/salmon/quant/{sample}/quant.sf", sample=SAMPLES),
    output:
        tx_rds="results/salmon_tximport/salmon_tx.rds.",
        gene_rds="results/salmon_tximport/salmon_gene.rds.",
        tx_tpm="results/salmon_tximport/salmon_tx_tpm.tsv",
        gene_tpm="results/salmon_tximport/salmon_gene_tpm.tsv",
        tx_count="results/salmon_tximport/salmon_tx_counts.tsv",
        gene_count="results/salmon_tximport/salmon_gene_counts.tsv",
    log:
        "logs/salmon/salmon/tximport.log",
    params:
        quantdir=lambda w, output: os.path.dirname(output[0]),
    conda:
        "../envs/tximport.yaml"
    script:
        "../scripts/merge_salmon_quant.R"
