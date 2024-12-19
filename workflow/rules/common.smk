import pandas as pd
import os
from snakemake.utils import validate


def remove_empty_values(dictionary: dict):
    """
    Recursively removes key-value pairs from dict for which the value is None
    Works for an arbitrarily deeply nested dict.
    See: https://github.com/snakemake/snakemake/issues/1701#issuecomment-1150741131
    """
    new_dictionary = dict()
    for key, value in dictionary.items():
        if isinstance(value, dict):
            new_dictionary[key] = remove_empty_values(value)
        elif value is None:
            pass
        else:
            new_dictionary[key] = value

    return new_dictionary


config = remove_empty_values(config)
validate(config, schema="../schemas/config_schema.yaml")


samples = pd.read_table(config["samples"], header=0, dtype=str, comment="#").set_index(
    "sample", drop=False
)
validate(samples, "../schemas/samples_schema.yaml")
SAMPLES = samples["sample"]


contrasts = pd.read_table(
    config["contrasts"], header=0, dtype=str, comment="#"
).set_index("contrast", drop=False)
validate(contrasts, "../schemas/contrasts_schema.yaml")
CONTRASTS = contrasts["contrast"]


def get_ref_file(path):
    if os.path.isfile(path):
        return path

    base_name = os.path.basename(path)
    fout = os.path.join("resources", base_name)
    return os.path.splitext(fout)[0]


def get_ref_url(path):
    if os.path.isfile(path):
        return
    return path


def get_transcriptome():
    return os.path.join("resources", "cdna_ncrna.fa")


def get_fastqs(wildcards):
    fqs = samples.loc[(wildcards.sample), ["fastq_1", "fastq_2"]]

    return {"fq1": f"{fqs.fastq_1}", "fq2": f"{fqs.fastq_2}"}


def get_rmats_group_indices(wildcards):
    with open("results/rmats/bam_config_post/all_samples.txt") as f:
        bam_files = f.readline().strip().split(",")

    sample_dict = samples.groupby("condition")["sample"].apply(list).to_dict()
    contrast = contrasts.loc[
        wildcards.contrast, ["reference_level", "tested_level"]
    ].values.tolist()
    group_indices = []

    for group in contrast:
        indices = []
        for sample in sample_dict[group]:
            for i, bam in enumerate(bam_files):
                if sample in bam:
                    indices.append(i)
        group_indices.append(",".join(map(str, indices)))

    return group_indices


def use_rmats_pairadise():
    return "--paired-stats" if config["rmats"]["paired_stats"] else ""


def get_deseq2_extra():
    extra = ""

    if config["deseq2"]["tpm"]:
        extra += "--tpm results/coco_cc/filtered_coco_tpm.tsv "

    if config["deseq2"]["reduced_model"]:
        extra += "--reduced-model " + config["deseq2"]["reduced_model"]

    if config["deseq2"]["control_genes"]:
        extra += "--control-genes " + config["deseq2"]["control_genes"]

    return extra


def get_multiqc_input(wildcards):

    input = (
        expand("results/fastp/{sample}_fastp.json", sample=SAMPLES)
        + expand(
            "results/star_align/{sample}/Aligned.sortedByCoord.out.bam", sample=SAMPLES
        )
        + expand("results/sam_stats/{sample}.txt", sample=SAMPLES)
        + expand("results/sam_flagstat/{sample}.tsv", sample=SAMPLES)
        + expand("results/sam_idxstats/{sample}.tsv", sample=SAMPLES)
    )
    if config["salmon"]["activate"]:
        input += expand(
            expand("results/salmon_quant/{sample}/quant.sf", sample=SAMPLES),
            sample=SAMPLES,
        )

    return input


def all_input(wildcards):
    wanted_input = []

    # fastp
    wanted_input.extend(expand("results/fastp/{sample}_fastp.html", sample=SAMPLES))
    # STAR
    wanted_input.extend(
        expand(
            "results/star_align/{sample}/Aligned.sortedByCoord.out.bam",
            sample=SAMPLES,
        )
    )
    # rseqc
    wanted_input.extend(
        expand("results/rseqc/infer_experiment/{sample}.txt", sample=SAMPLES)
    )
    # Samtools
    wanted_input.extend(
        expand("results/sam_stats/{sample}.txt", sample=SAMPLES)
        + expand("results/sam_flagstat/{sample}.tsv", sample=SAMPLES)
        + expand("results/sam_idxstats/{sample}.tsv", sample=SAMPLES)
    )
    # salmon
    if config["salmon"]["activate"]:
        wanted_input.extend(["results/salmon_tximport/salmon_gene_tpm.tsv"])
    # CoCo
    if config["coco"]["activate"]:
        wanted_input.extend(
            expand(
                "results/coco_cc/coco_{counttype}.tsv",
                counttype=["counts", "cpm", "tpm"],
            )
            + ["results/coco_cc/filtered_coco_tpm.tsv"]
            + expand("results/coco_cb/{sample}.bigwig", sample=SAMPLES)
        )
    # rMATS
    if config["rmats"]["activate"]:
        wanted_input.extend(
            expand(
                "results/rmats/stat/{contrast}/{filter}_{event}.MATS.{junction}.txt",
                contrast=CONTRASTS,
                filter=["filtered", "up", "dn", "bg"],
                event=["SE", "MXE", "RI", "A5SS", "A3SS"],
                junction=["JC"],
            )
        )
    # DESeq2
    if config["deseq2"]["activate"]:
        wanted_input.extend(["results/deseq2/heatmap_samples.png"])
    # MultiQC
    wanted_input.extend(["results/multiqc"])

    return wanted_input


def download_input(wildcards):
    wanted_input = []

    # ref
    wanted_input.extend(
        [get_ref_file(ref) for ref in config["ref"].values()] + [get_transcriptome()]
    )
    # CoCo
    wanted_input.extend(["resources/coco", "resources/pairedBamToBed12/bin"])
    # rMATS
    wanted_input.extend(["resources/install_r_deps.R", "resources/PAIRADISE"])

    return wanted_input
