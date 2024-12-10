import argparse
import os
import pandas as pd
import subprocess
import sys


def validate_file(file_path):
    if not os.path.isfile(file_path):
        raise argparse.ArgumentTypeError(f"The file does not exist: {file_path}")
    return file_path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Arguments for the analysis script.", add_help=False
    )
    required = parser.add_argument_group("required arguments")
    optional = parser.add_argument_group("optional arguments")

    optional.add_argument(
        "-h",
        "--help",
        action="help",
        default=argparse.SUPPRESS,
        help="Show this help message and exit",
    )

    required.add_argument(
        "-o",
        "--outdir",
        type=str,
        required=True,
        help="Path to the output directory where results will be saved (default: current working directory)",
    )

    required.add_argument(
        "-c",
        "--counts",
        type=validate_file,
        required=True,
        help="Path to the raw counts file containing gene count data",
    )

    required.add_argument(
        "-s",
        "--samples",
        type=validate_file,
        required=True,
        help="Path to the sample sheet which includes metadata about samples",
    )

    required.add_argument(
        "-r",
        "--contrasts",
        type=validate_file,
        required=True,
        help="Path to the contrast sheet specifying the contrasts for differential analysis",
    )

    required.add_argument(
        "-f",
        "--full-model",
        type=str,
        required=True,
        help="Full design formula with no spaces starting with '~' for the Wald test",
    )

    optional.add_argument(
        "-d",
        "--reduced-model",
        type=str,
        help="Reduced design formula with no spaces starting with '~' for the likelihood ratio test (LRT)",
    )

    optional.add_argument(
        "-t",
        "--tpm",
        type=validate_file,
        help="Path to the TPM (Transcripts Per Million) file for filtering gene expression",
    )

    optional.add_argument(
        "-g",
        "--control-genes",
        type=validate_file,
        help="Path to a file containing one gene ID per line for estimating size factors",
    )

    optional.add_argument(
        "-p",
        "--padj",
        type=str,
        default="0.05",
        help="Adjusted p-value threshold for significance (default: %(default)s)",
    )

    optional.add_argument(
        "-l",
        "--lfc",
        type=str,
        default="1",
        help="Absolute value of log2 fold change threshold for significance (default: %(default)s)",
    )

    optional.add_argument(
        "-m",
        "--min-gene-expr",
        type=str,
        default="10",
        help="Minimum gene expression level required for a gene to be considered (default: %(default)s)",
    )

    optional.add_argument(
        "-n",
        "--min-samps-gene-expr",
        type=str,
        help="Minimum number of samples in which a gene must be expressed to be considered (default: minimal number of replicates in any of the conditions)",
    )

    return parser.parse_args()


def design_starts_with_tilde(design):
    if not design.startswith("~"):
        print("Error: The formula must start with '~'", file=sys.stderr)
        sys.exit(1)


def design_in_sample_sheet(design, samples):
    with open(samples, "r") as f:
        header = f.readline().strip().split("\t")
        variables = design.split("~")[1].split("+")

    # In case reduced design is '~1'
    variables = [var for var in variables if var != "1"]
    missing_var = [var for var in variables if var not in header]

    if missing_var:
        print(
            f"Error: The following variables are missing in the sample sheet columns: {', '.join(missing_var)}",
            file=sys.stderr,
        )
        sys.exit(1)


def is_tsv_file(file):
    with open(file, "r") as f:
        first_line = f.readline().strip()

    if "\t" not in first_line:
        print(f"Error: '{file}' is not a TSV file", file=sys.stderr)
        sys.exit(1)


def check_required_columns(file, required_columns):
    with open(file, "r") as f:
        header = f.readline().strip()

    columns = header.split("\t")

    missing_columns = [col for col in required_columns if col not in columns]
    if missing_columns:
        print(
            f"Error: Missing required columns: {', '.join(missing_columns)}",
            file=sys.stderr,
        )
        sys.exit(1)


def validate_sample_ids(samples, counts):
    samples_df = pd.read_table(samples)
    counts_df = pd.read_table(counts)
    counts_df = counts_df.select_dtypes(include=["number"])
    sample_id_diff = set(samples_df["sample"]).difference(counts_df.columns)
    if len(sample_id_diff) != 0:
        print(
            "Error: Sample IDs in the sample sheet and the counts file are not the same",
            file=sys.stderr,
        )
        sys.exit(1)


def validate_contrasts(samples, contrasts):
    samples_df = pd.read_table(samples)
    countrasts_df = pd.read_table(contrasts)
    contrast_groups = pd.unique(
        countrasts_df[["reference_level", "tested_level"]].values.ravel()
    )
    samples_df_groups = pd.unique(samples_df["condition"])

    if set(contrast_groups).difference(samples_df_groups):
        print(
            "Conditions in the sample sheet and contrast sheet are not the same",
            file=sys.stderr,
        )
        sys.exit(1)


def get_min_samps_gene_expr(samples):
    conditions_df = pd.read_table(samples)["condition"]
    min_samples_per_group = str(conditions_df.value_counts().min())

    return min_samples_per_group


def get_script_dir():
    script_abs_path = os.path.abspath(sys.argv[0])
    script_dir = os.path.dirname(script_abs_path)
    return script_dir


def run_command(command):
    print(f"running: {' '.join(command)}")
    subprocess.run(command, check=True)


if __name__ == "__main__":
    args = parse_args()

    tsv_files = [args.counts, args.samples, args.contrasts, args.tpm]
    for file in tsv_files:
        if file is not None:
            is_tsv_file(file)

    required_columns = {
        args.samples: ["condition"],
        args.contrasts: ["contrast", "reference_level", "tested_level"],
    }

    for file, col in required_columns.items():
        check_required_columns(file, col)

    for design in [args.full_model, args.reduced_model]:
        if design is not None:
            design_starts_with_tilde(design)
            design_in_sample_sheet(design, args.samples)

    validate_sample_ids(args.samples, args.counts)
    validate_contrasts(args.samples, args.contrasts)

    if args.min_samps_gene_expr is None:
        args.min_samps_gene_expr = get_min_samps_gene_expr(args.samples)

    r_script_dir = get_script_dir()
    r_script_path = os.path.join(r_script_dir, "deseq2.R")

    command = [
        "Rscript",
        r_script_path,
        args.outdir,
        args.counts,
        args.samples,
        args.contrasts,
        args.full_model,
        args.reduced_model,
        args.tpm,
        args.control_genes,
        args.padj,
        args.lfc,
        args.min_gene_expr,
        args.min_samps_gene_expr,
    ]
    command = [arg if arg is not None else "" for arg in command]
    run_command(command)
