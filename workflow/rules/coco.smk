localrules:
    install_coco,
    install_pairedBamToBed12,


rule install_coco:
    output:
        directory("resources/coco"),
        "resources/coco/bin/coco.py",
    params:
        link="https://github.com/scottgroup/coco.git",
    conda:
        "../envs/coco.yaml"
    shell:
        "cd resources/"
        " && git clone {params.link}"


rule install_pairedBamToBed12:
    output:
        directory("resources/pairedBamToBed12/bin"),
    params:
        link="https://github.com/Population-Transcriptomics/pairedBamToBed12.git",
    conda:
        "../envs/coco.yaml"
    shell:
        "cd resources/ "
        " && git clone {params.link}"
        " && cd pairedBamToBed12/"
        " && make "
        " && export PATH=$PWD/bin:$PATH"


rule coco_ca:
    input:
        coco="resources/coco/bin/coco.py",
        gtf=get_gtf(),
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
        gtf=get_gtf(),
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


rule coco_cb:
    input:
        coco="resources/coco/bin/coco.py",
        bam="results/star_align/{sample}/Aligned.sortedByCoord.out.bam",
        chrNameLength="results/star_index/chrNameLength.txt",
        pb2b="resources/pairedBamToBed12/bin",
    output:
        unsorted_bedgraph="results/coco_cb/bedgraph/{sample}_unsorted.bedgraph.bed12",
    log:
        "logs/coco/coco_cb/{sample}.log",
    conda:
        "../envs/coco.yaml"
    shell:
        "python {input.coco} cb"
        " --ucsc_compatible"
        " --thread {resources.threads}"
        " {input.bam}"
        " {output.unsorted_bedgraph}"
        " {input.chrNameLength}"
        " &> {log}"


rule coco_sort_bg:
    input:
        unsorted_bedgraph="results/coco_cb/bedgraph/{sample}_unsorted.bedgraph.bed12",
    output:
        sorted_bedgraph="results/coco_cb/bedgraph/{sample}.bedgraph",
    log:
        "logs/coco/sort_bg/{sample}.log",
    shell:
        "sort -k1,1 -k2,2n {input.unsorted_bedgraph}"
        " | sed 's/chrM/chrMT/g' > {output.sorted_bedgraph}"
        " && rm {input.unsorted_bedgraph}"
        " &> {log}"
