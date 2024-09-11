rule trimmomatic:
    input:
        unpack(get_fastqs),
    output:
        fq1="results/trimmomatic/{sample}_R1.fastq.gz",
        fq2="results/trimmomatic/{sample}_R2.fastq.gz",
        unpaired_fq1="results/trimmomatic/{sample}_R1.unpaired.fastq.gz",
        unpaired_fq2="results/trimmomatic/{sample}_R2.unpaired.fastq.gz",
    log:
        "logs/trimmomatic/{sample}.log",
    params:
        options=[
            "ILLUMINACLIP:TruSeq3-PE-2.fa:2:12:10:8:true",
            "TRAILING:30",
            "LEADING:30",
            "MINLEN:20",
        ],
    conda:
        "../envs/trimmomatic.yaml"
    shell:
        "trimmomatic PE"
        " -threads {resources.threads}"
        " -phred33 "
        " {input.fq1} {input.fq2}"
        " {output.fq1} {output.unpaired_fq1}"
        " {output.fq2} {output.unpaired_fq2}"
        " {params.options}"
        " &> {log}"
