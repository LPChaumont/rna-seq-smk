rule fastp:
    input:
        unpack(get_fastqs),
    output:
        fq1="results/fastp/{sample}_trimmed_R1.fastq.gz",
        fq2="results/fastp/{sample}_trimmed_R2.fastq.gz",
        ufq1="results/fastp/{sample}_unpaired_R1.fastq.gz",
        ufq2="results/fastp/{sample}_unpaired_R2.fastq.gz",
        json="results/fastp/{sample}_fastp.json",
        html="results/fastp/{sample}_fastp.html",
    log:
        "logs/fastp/{sample}.log",
    conda:
        "../envs/fastp.yaml"
    shell:
        "fastp"
        " -i {input.fq1} -I {input.fq2}"
        " -o {output.fq1} -O {output.fq2}"
        " --unpaired1 {output.ufq1} --unpaired2 {output.ufq2}"
        " --json {output.json}"
        " --html {output.html}"
        " --thread {resources.threads}"
        " 2> {log}"


rule fastp_multiqc:
    input:
        data=expand("results/fastp/{sample}_fastp.json", sample=SAMPLES),
    output:
        outdir=directory("results/multiqc/fastp"),
        report="results/multiqc/fastp/multiqc_fastp_report.html",
    params:
        filename="multiqc_fastp_report.html",
        indir="results/fastp",
    log:
        "logs/multiqc/fastp.log",
    conda:
        "../envs/multiqc.yaml"
    shell:
        "multiqc"
        " {params.indir}"
        " --outdir {output.outdir}"
        " --filename {params.filename}"
        " --force"
        " &> {log}"
