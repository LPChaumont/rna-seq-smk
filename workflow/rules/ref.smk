localrules:
    get_genome,
    get_gtf,
    get_transcriptome,
    build_full_transcriptome,


rule get_genome:
    output:
        get_genome(),
    params:
        species=SPECIES,
        release=RELEASE,
        build=BUILD,
        datatype="dna",
    log:
        "logs/ref/get_genome.log",
    wrapper:
        "v3.14.1/bio/reference/ensembl-sequence"


rule get_gtf:
    output:
        get_gtf(),
    params:
        species=SPECIES,
        release=RELEASE,
        build=BUILD,
    log:
        "logs/ref/get_gtf.log",
    wrapper:
        "v3.14.1/bio/reference/ensembl-annotation"


rule get_transcriptome:
    output:
        ".".join(["resources/ensembl", SPECIES, BUILD, RELEASE, "{datatype}", "fa"]),
    params:
        species=SPECIES,
        release=RELEASE,
        build=BUILD,
        datatype="{datatype}",
    log:
        "logs/ref/get_{datatype}_transcriptome.log",
    wrapper:
        "v1.7.1/bio/reference/ensembl-sequence"


rule build_full_transcriptome:
    input:
        expand(
            ".".join(
                [
                    "resources/ensembl",
                    SPECIES,
                    BUILD,
                    RELEASE,
                    "{datatype}",
                    "fa",
                ]
            ),
            datatype=TRANSCRIPTOME_DATATYPES,
        ),
    output:
        get_full_transcriptome(),
    log:
        "logs/ref/build_full_transcriptome.log",
    shell:
        "cat {input} > {output}"
