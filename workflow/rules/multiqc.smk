rule multiqc:
    input:
        get_multiqc_input,
    output:
        directory("results/multiqc"),
    log:
        "logs/multiqc/multiqc.log",
    conda:
        "../envs/multiqc.yaml"
    shell:
        """
        multiqc results
        --outdir {output}
        --force
        --verbose
        &> {log}
        """
