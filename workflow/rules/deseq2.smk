rule deseq2:
    input:
        raw_counts="results/coco_cc/coco_counts.tsv",
        filtered_tpm="results/coco_cc/filtered_coco_tpm.tsv",
        samples=config["samples"],
        contrasts=config["contrasts"],
    output:
        rds="results/deseq2/heatmap_samples.png",
    log:
        "logs/deseq2.log",
    params:
        outdir=lambda w, output: os.path.dirname(output[0]),
        padj=config["deseq2"]["padj"],
        lfc=config["deseq2"]["lfc"],
        min_gene_expr=config["deseq2"]["min_gene_expr"],
        min_samps_gene_expr=config["deseq2"]["min_samps_gene_expr"],
        full_model=config["deseq2"]["full_model"],
        reduced_model=config["deseq2"]["reduced_model"],
        extra=get_deseq2_extra(),
    conda:
        "../envs/deseq2.yaml"
    shell:
        "python workflow/scripts/deseq2.py"
        " --outdir {params.outdir}"
        " --counts {input.raw_counts}"
        " --samples {input.samples}"
        " --contrasts {input.contrasts}"
        " --full-model {params.full_model}"
        " --padj {params.padj}"
        " --lfc {params.lfc}"
        " --min-gene-expr {params.min_gene_expr}"
        " --min-samps-gene-expr {params.min_samps_gene_expr}"
        " {params.extra}"
