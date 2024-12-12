localrules:
    get_gtf,
    get_dna_fasta,
    get_cdna_fasta,
    get_ncrna_fasta,
    concat_cdna_ncrna_fasta,


rule get_gtf:
    output:
        get_ref_file(config["ref"]["gtf"]),
    params:
        lambda _: config["ref"]["gtf"],
    log:
        "logs/ref/get_gtf.log",
    conda:
        "../envs/ref.yaml"
    shell:
        "curl -L {params} | gzip -d >> {output} 2> {log}"


rule get_dna_fasta:
    output:
        get_ref_file(config["ref"]["dna_fasta"]),
    params:
        lambda _: config["ref"]["dna_fasta"],
    log:
        "logs/ref/get_dna_fasta.log",
    conda:
        "../envs/ref.yaml"
    shell:
        "curl -L {params} | gzip -d >> {output} 2> {log}"


rule get_cdna_fasta:
    output:
        get_ref_file(config["ref"]["cdna_fasta"]),
    params:
        lambda _: config["ref"]["cdna_fasta"],
    log:
        "logs/ref/get_cdna_fasta.log",
    conda:
        "../envs/ref.yaml"
    shell:
        "curl -L {params} | gzip -d >> {output} 2> {log}"


rule get_ncrna_fasta:
    output:
        get_ref_file(config["ref"]["ncrna_fasta"]),
    params:
        lambda _: config["ref"]["ncrna_fasta"],
    log:
        "logs/ref/get_ncrna_fasta.log",
    conda:
        "../envs/ref.yaml"
    shell:
        "curl -L {params} | gzip -d >> {output} 2> {log}"


rule concat_cdna_ncrna_fasta:
    input:
        cdna=get_ref_file(config["ref"]["cdna_fasta"]),
        ncrna=get_ref_file(config["ref"]["ncrna_fasta"]),
    output:
        get_transcriptome(),
    log:
        "logs/ref/concat_cdna_ncrna_fasta.log",
    conda:
        "../envs/ref.yaml"
    shell:
        "cat {input.cdna} {input.ncrna} > {output}"
