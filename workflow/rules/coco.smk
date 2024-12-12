localrules:
    install_coco,
    install_pairedBamToBed12,


rule install_coco:
    output:
        directory("resources/coco"),
        "resources/coco/bin/coco.py",
    params:
        link="https://github.com/scottgroup/coco.git",
    log:
        "logs/coco/install_coco.log",
    conda:
        "../envs/coco.yaml"
    shell:
        "cd resources/"
        " && git clone {params.link} &> {log}"


rule install_pairedBamToBed12:
    output:
        directory("resources/pairedBamToBed12/bin"),
    params:
        link="https://github.com/Population-Transcriptomics/pairedBamToBed12.git",
    log:
        "logs/coco/install_pairedBamToBed12.log",
    conda:
        "../envs/coco.yaml"
    shell:
        "cd resources/"
        " && git clone {params.link} &> {log}"
        " && cd pairedBamToBed12/"
        " && make &>> {log}"


rule coco_ca:
    input:
        coco="resources/coco/bin/coco.py",
        gtf=get_ref_file(config["ref"]["gtf"]),
    output:
        coco_gtf="results/coco_ca/correct_annotation.gtf",
    log:
        "logs/coco/coco_ca.log",
    conda:
        "../envs/coco.yaml"
    shell:
        "python {input.coco} ca {input.gtf} -o {output.coco_gtf}"
        " &> {log}"


rule coco_cc:
    input:
        coco="resources/coco/bin/coco.py",
        coco_gtf="results/coco_ca/correct_annotation.gtf",
        bam="results/star_align/{sample}/Aligned.sortedByCoord.out.bam",
    output:
        quant="results/coco_cc/{sample}.tsv",
    log:
        "logs/coco/coco_cc/{sample}.log",
    conda:
        "../envs/coco.yaml"
    shell:
        "python {input.coco} cc"
        " --countType both"
        " --thread {resources.threads}"
        " --strand 1"
        " --paired"
        " {input.coco_gtf}"
        " {input.bam}"
        " {output.quant}"
        " &> {log}"


rule merge_coco_quant:
    input:
        gtf=get_ref_file(config["ref"]["gtf"]),
        quants=expand("results/coco_cc/{sample}.tsv", sample=SAMPLES),
    output:
        counts="results/coco_cc/coco_counts.tsv",
        cpm="results/coco_cc/coco_cpm.tsv",
        tpm="results/coco_cc/coco_tpm.tsv",
    log:
        "logs/coco/merge_coco_quant.log",
    conda:
        "../envs/coco.yaml"
    script:
        "../scripts/merge_coco_quant.py"


rule coco_filter_tpm:
    input:
        tpm="results/coco_cc/coco_tpm.tsv",
        samples=config["samples"],
    output:
        filtered="results/coco_cc/filtered_coco_tpm.tsv",
    log:
        "logs/coco/filter_coco_quant.log",
    conda:
        "../envs/coco.yaml"
    shell:
        "python workflow/scripts/filter_gene_expression.py"
        " --input {input.tpm}"
        " --samples {input.samples}"
        " --min_gene_expr 1"
        " --min_condition_expr 1"
        " --save"
        " &> {log}"


rule coco_cb:
    input:
        coco="resources/coco/bin/coco.py",
        bam="results/star_align/{sample}/Aligned.sortedByCoord.out.bam",
        chrNameLength="results/star_index/chrNameLength.txt",
    output:
        unsorted_bedgraph="results/coco_cb/{sample}_unsorted.bedgraph",
    params:
        pb2b="resources/pairedBamToBed12/bin",
    log:
        "logs/coco/coco_cb/{sample}.log",
    conda:
        "../envs/coco.yaml"
    shell:
        "export PATH=$PATH:$PWD/{params.pb2b} &&"
        " python {input.coco} cb"
        " --ucsc_compatible"
        " --thread {resources.threads}"
        " {input.bam}"
        " {output.unsorted_bedgraph}"
        " {input.chrNameLength}"
        " &> {log}"


rule coco_sort_bg:
    input:
        unsorted_bedgraph="results/coco_cb/{sample}_unsorted.bedgraph",
    output:
        sorted_bedgraph="results/coco_cb/{sample}_sorted.bedgraph",
    log:
        "logs/coco/sort_bg/{sample}.log",
    conda:
        "../envs/coco.yaml"
    shell:
        "sort -k1,1 -k2,2n {input.unsorted_bedgraph}"
        " | sed 's/chrM/chrMT/g' > {output.sorted_bedgraph}"
        " && rm {input.unsorted_bedgraph}"


rule chromsize:
    input:
        fasta=get_ref_file(config["ref"]["dna_fasta"]),
    output:
        chromsize="resources/chrom_sizes.tsv",
    log:
        "logs/chromsize.log",
    conda:
        "../envs/bedgraphtobigwig.yaml"
    shell:
        """
        chromsize --fasta {input.fasta} --output {output.chromsize} &> {log}
        awk '{{print "chr"$1 "\t" $NF}}' {output.chromsize} | sed 's/\bchrM\b/chrMT/g' > {output.chromsize}.tmp
        mv {output.chromsize}.tmp {output.chromsize}
        """


rule coco_bg2bw:
    input:
        chromsize="resources/chrom_sizes.tsv",
        bg="results/coco_cb/{sample}_sorted.bedgraph",
    output:
        bw="results/coco_cb/{sample}.bigwig",
    log:
        "logs/coco/bg2bw/{sample}.log",
    conda:
        "../envs/bedgraphtobigwig.yaml"
    shell:
        "sed '$d' {input.bg} > {input.bg}.tmp && bedGraphToBigWig {input.bg}.tmp {input.chromsize} {output.bw} &> {log}"
        " && rm {input.bg}.tmp"
