rule chromsize:
    input:
        get_genome(),
    output:
        "results/coco_cb/bigwig/chrom_sizes.txt",
    log:
        "logs/bedgraphtobigwig/chromsize.log",
    conda:
        "../envs/bedgraphtobigwig.yaml"
    shell:
        "chromsize --fasta {input} --output {output}"
        " --threads {resources.threads}"
        " 2> {log}"


rule bgtobw:
    input:
        chromsizes="results/coco_cb/bigwig/chrom_sizes.txt",
        bedgraph="results/coco_cb/bedgraph/{sample}.bedgraph",
    output:
        bigwig="results/coco_cb/bigwig/{sample}.bigwig",
    log:
        "logs/bedgraphtobigwig/{sample}.log",
    conda:
        "../envs/bedgraphtobigwig.yaml"
    shell:
        "bedGraphToBigWig"
        " {input.bedgraph} {input.chromsizes}"
        " 2> {log}"
