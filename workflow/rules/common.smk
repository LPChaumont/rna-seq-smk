import pandas as pd

samples = pd.read_table(config["samples"], header=0, dtype=str, comment="#").set_index(
    "sample", drop=False
)
samples.index.names = ["sample_id"]
SAMPLES = samples["sample"]


contrasts = pd.read_table(
    config["contrasts"], header=0, dtype=str, comment="#"
).set_index("contrast", drop=False)
CONTRASTS = contrasts["contrast"]

SPECIES = config["ref"]["species"]
BUILD = config["ref"]["build"]
RELEASE = config["ref"]["release"]
TRANSCRIPTOME_DATATYPES = ["cdna", "ncrna"]


def get_genome():
    if config["ref"]["custom_genome"]:
        return config["ref"]["custom_genome"]

    return f"resources/ensembl.{SPECIES}.{BUILD}.primary_assembly.fa"


def get_gtf():
    if config["ref"]["custom_gtf"]:
        return config["ref"]["custom_gtf"]

    return f"resources/ensembl.{SPECIES}.{BUILD}.{RELEASE}.gtf"


def get_transcriptome(wildcards):
    if config["ref"]["custom_transcriptome"]:
        return config["ref"]["custom_transcriptome"]

    return f"resources/ensembl.{SPECIES}.{BUILD}.{RELEASE}.{wildcards.datatype}.fa"


def get_full_transcriptome():
    if config["ref"]["custom_transcriptome"]:
        return config["ref"]["custom_transcriptome"]

    return (
        f"resources/ensembl.{SPECIES}.{BUILD}.{RELEASE}."
        + ".".join(TRANSCRIPTOME_DATATYPES)
        + ".fa"
    )


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
    # trimmomatic
    wanted_input.extend(
        expand(
            "results/trimmomatic/{sample}_{read}.fastq.gz",
            sample=SAMPLES,
            read=["R1", "R2"],
        )
    )
    # STAR
    wanted_input.extend(
        expand(
            "results/star_align/{sample}/Aligned.sortedByCoord.out.bam",
            sample=SAMPLES,
        )
    )
    # Samtools
    wanted_input.extend(
        expand("results/sam_stats/{sample}.txt", sample=SAMPLES)
        + expand("results/sam_flagstat/{sample}.tsv", sample=SAMPLES)
        + expand("results/sam_idxstats/{sample}.tsv", sample=SAMPLES)
    )
    # salmon
    if config["salmon"]["activate"]:
        wanted_input.extend(
            [
                "results/salmon_quant/salmon_tpm_gene.tsv",
            ]
        )
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
    wanted_input.extend([get_full_transcriptome(), get_genome(), get_gtf()])
    # CoCo
    wanted_input.extend(["resources/coco", "resources/pairedBamToBed12/bin"])
    # rMATS
    wanted_input.extend(["resources/install_r_deps.R", "resources/PAIRADISE"])

    return wanted_input
