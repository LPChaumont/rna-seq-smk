import argparse
import math
import os
from rmats_class_exon import exon_SE, exon_RI, exon_AXSS, exon_MXE

# modified from https://github.com/Xinglab/rmats-turbo-tutorial/blob/main/scripts/rmats_filtering.py

EXON_CLASS_MAP = {
    "SE": exon_SE,
    "RI": exon_RI,
    "A3SS": exon_AXSS,
    "A5SS": exon_AXSS,
    "MXE": exon_MXE,
}

SUFFIX = ".MATS.JC.txt"


def parse_args():
    parser = argparse.ArgumentParser(description="Filter rMATS output")

    parser.add_argument(
        "--rmats-dir",
        type=str,
        required=True,
        help="Path to the rMATS output directory",
    )
    parser.add_argument(
        "--outdir",
        type=str,
        help="Output directory to write the filtered results (default: same as --rmats-dir)",
    )
    parser.add_argument(
        "--read-cov",
        type=int,
        default=10,
        help="Minimum read coverage (default: %(default)s)",
    )
    parser.add_argument(
        "--min-psi",
        type=float,
        default=0.05,
        help="Minimum average PSI (default: %(default)s)",
    )
    parser.add_argument(
        "--max-psi",
        type=float,
        default=0.95,
        help="Maximum average PSI (default: %(default)s)",
    )
    parser.add_argument(
        "--sig-fdr",
        type=float,
        default=0.05,
        help="Threshold for significant FDR (default: %(default)s)",
    )
    parser.add_argument(
        "--sig-delta-psi",
        type=float,
        default=0.1,
        help="Threshold for significant delta PSI (default: %(default)s)",
    )
    parser.add_argument(
        "--bg-fdr",
        type=float,
        default=0.5,
        help="Threshold for background FDR (default: %(default)s)",
    )
    parser.add_argument(
        "--bg-within-group-delta-psi",
        type=float,
        default=0.2,
        help="Threshold for background within-group delta PSI (default: %(default)s)",
    )

    return parser.parse_args()


def get_exon_class(file_name):
    event_type = os.path.basename(file_name)[: -len(SUFFIX)]

    return EXON_CLASS_MAP.get(event_type)


def filter_read_coverage_psi(event, read_cov, min_psi, max_psi):
    coverage_condition = (
        event.averageIJC_SAMPLE_1 >= read_cov or event.averageSJC_SAMPLE_1 >= read_cov
    ) and (
        event.averageIJC_SAMPLE_2 >= read_cov or event.averageSJC_SAMPLE_2 >= read_cov
    )

    psi_condition = (
        min(event.averagePsiSample1, event.averagePsiSample2) <= max_psi
        and max(event.averagePsiSample1, event.averagePsiSample2) >= min_psi
    )

    return coverage_condition and psi_condition


def filter_upregulated_sig(event, sig_fdr, sig_delta_psi):
    return event.FDR < sig_fdr and event.IncLevelDifference >= sig_delta_psi


def filter_downregulated_sig(event, sig_fdr, sig_delta_psi):
    return event.FDR < sig_fdr and event.IncLevelDifference <= -sig_delta_psi


def filter_background(event, bg_fdr, bg_within_group_delta_psi):
    return (
        event.FDR >= bg_fdr
        and max(event.IncLevel1) - min(event.IncLevel1) <= bg_within_group_delta_psi
        and max(event.IncLevel2) - min(event.IncLevel2) <= bg_within_group_delta_psi
    )


def filter_without_b2(event, read_cov, min_psi, max_psi):
    if (
        math.isnan(event.averageIJC_SAMPLE_2)
        and math.isnan(event.averageSJC_SAMPLE_2)
        and math.isnan(event.averagePsiSample2)
    ):

        return (
            event.averageIJC_SAMPLE_1 >= read_cov
            or event.averageSJC_SAMPLE_1 >= read_cov
        ) and min_psi <= event.averagePsiSample1 <= max_psi


def filter_rMATs(
    file_name,
    read_cov,
    min_psi,
    max_psi,
    sig_fdr,
    sig_delta_psi,
    bg_fdr,
    bg_within_group_delta_psi,
):
    event_dict = {"filtered": [], "up": [], "dn": [], "bg": []}

    exon = get_exon_class(file_name)

    if exon is None:
        print(f"Not filtering: {file_name}")
        return

    with open(file_name, "r") as in_handle:
        header = in_handle.readline()

        for line in in_handle:
            event = exon(line)
            
            # Filter event based on read coverage and PSI values
            if filter_read_coverage_psi(event, read_cov, min_psi, max_psi):
                event_dict["filtered"].append(event)

                if filter_downregulated_sig(event, sig_fdr, sig_delta_psi):
                    event_dict["dn"].append(event)
                elif filter_upregulated_sig(event, sig_fdr, sig_delta_psi):
                    event_dict["up"].append(event)
                elif filter_background(event, bg_fdr, bg_within_group_delta_psi):
                    event_dict["bg"].append(event)

            # Handle events where --b2 is not provided
            elif filter_without_b2(event, read_cov, min_psi, max_psi):
                event_dict["filtered"].append(event)

    return event_dict


def get_header(file_name):
    with open(file_name, "r") as in_handle:
        header = in_handle.readline()
    return header


def write_event(outdir, in_file_name, event_dict, header):
    os.makedirs(outdir, exist_ok=True)

    for category, events in event_dict.items():
        out_file_name = os.path.join(outdir, f"{category}_{in_file_name}")

        with open(out_file_name, "w") as out_handle:
            out_handle.write(header)

            for event in events:
                out_handle.write(str(event))


def write_summary(outdir, summary_dict):
    summary_file = os.path.join(outdir, "summary_filtered.tsv")
    summary_order = ["SE", "A3SS", "A5SS", "MXE", "RI"]
    summary_dict = {k: summary_dict[k] for k in summary_order if k in summary_dict}

    with open(summary_file, "w") as out_handle:
        header = "\t".join(["event_type", "filtered", "total_sig", "up", "dn"]) + "\n"
        out_handle.write(header)
        for k, v in summary_dict.items():
           out_handle.write(f"{k}\t" + "\t".join(map(str, v)) + "\n")


def main():
    args = parse_args()
    args_dict = vars(args)

    rmats_dir = args_dict.pop("rmats_dir")
    outdir = args_dict.pop("outdir") or rmats_dir
    summary_dict = {}

    file_names = sorted(os.listdir(rmats_dir))
    for name in file_names:
        if not name.endswith(SUFFIX):
            continue
        
        # filtering
        file_path = os.path.join(rmats_dir, name)
        event_dict = filter_rMATs(file_path, *args_dict.values())

        if event_dict is None:
            continue

        header = get_header(file_path)
        write_event(outdir, name, event_dict, header)

        # summary
        filtered = len(event_dict["filtered"])
        up = len(event_dict["up"])
        dn = len(event_dict["dn"])
        total_sig = up + dn

        event_type = name.split(".")[0]
        summary_dict[event_type] = [filtered, total_sig, up, dn]
        write_summary(outdir, summary_dict)


if __name__ == "__main__":
    main()
