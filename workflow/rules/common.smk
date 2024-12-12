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


def get_ref_file(file_path):
    def remove_gz(path):
        return path[:-3] if path.endswith(".gz") else path

    if os.path.isfile(file_path):
        return file_path

    return remove_gz(os.path.join("resources", os.path.basename(file_path)))


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
    wanted_input.extend(
        [get_ref_file(ref) for ref in config["ref"].values()] + [get_transcriptome()]
    )
    # CoCo
    wanted_input.extend(["resources/coco", "resources/pairedBamToBed12/bin"])
    # rMATS
    wanted_input.extend(["resources/install_r_deps.R", "resources/PAIRADISE"])

    return wanted_input
