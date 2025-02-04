#!/usr/bin/env python3

import argparse
import math
import os
import sys
from rmats_class_exon import exon_SE, exon_RI, exon_AXSS, exon_MXE

# modified from https://github.com/Xinglab/rmats-turbo-tutorial/blob/main/scripts/rmats_filtering.py

def get_exon_class(fn):
    if "SE" in fn:
        return exon_SE
    elif "RI" in fn:
        return exon_RI
    elif "A3SS" in fn:
        return exon_AXSS
    elif "A5SS" in fn:
        return exon_AXSS
    elif "MXE" in fn:
        return exon_MXE
    else:
        print("Invalid alternative event type in the input file name. Please modify the file name.")
        sys.exit()


def filter_rMATS(
    fn, read_cov, min_psi, max_psi, sig_fdr, bg_fdr, sig_delta_psi, bg_within_group_delta_psi
):

    exon = get_exon_class(fn)

    filtered_events = []
    upregulated_events = []
    downnregulated_events = []
    background_events = []

    with open(fn, "r") as f:
        header = f.readline()
        for line in f:
            event = exon(line)
            # filter events by read coverage and PSI
            if (
                (event.averageIJC_SAMPLE_1 >= read_cov or event.averageSJC_SAMPLE_1 >= read_cov)
                and (event.averageIJC_SAMPLE_2 >= read_cov or event.averageSJC_SAMPLE_2 >= read_cov)
                and min(event.averagePsiSample1, event.averagePsiSample2) <= max_psi
                and max(event.averagePsiSample1, event.averagePsiSample2) >= min_psi
            ):
                filtered_events.append(event)
                # filter events by delta PSI and FDR
                if event.FDR < sig_fdr:
                    if event.IncLevelDifference >= sig_delta_psi:
                        downnregulated_events.append(event)
                    elif event.IncLevelDifference <= -sig_delta_psi:
                        upregulated_events.append(event)
                elif (
                    event.FDR >= bg_fdr
                    and max(event.IncLevel1) - min(event.IncLevel1) <= bg_within_group_delta_psi
                    and max(event.IncLevel2) - min(event.IncLevel2) <= bg_within_group_delta_psi
                ):
                    background_events.append(event)
            # handle cases where --b2 is not provided.
            elif (
                (event.averageIJC_SAMPLE_1 >= read_cov or event.averageSJC_SAMPLE_1 >= read_cov)
                and math.isnan(event.averageIJC_SAMPLE_2)
                and min(event.averagePsiSample1, event.averagePsiSample2) <= max_psi
                and max(event.averagePsiSample1, event.averagePsiSample2) >= min_psi
            ):
                filtered_events.append(event)

    event_dict = {
        "upregulated": upregulated_events,
        "downregulated": downnregulated_events,
        "filtered": filtered_events,
    }

    return header, event_dict


def parse_args():
    parser = argparse.ArgumentParser(description="Filter rMATS output")
    
    parser.add_argument("--input", type=str, required=True, help="Path to the rMATS output file")
    parser.add_argument("--outdir", type=str, default=".", help="Output directory (default: %(default)s)")
    parser.add_argument("--read_cov", type=int, default=10, help="Minimum read coverage (default: %(default)s)")
    parser.add_argument("--min_psi", type=float, default=0.05, help="Minimum average PSI (default: %(default)s)")
    parser.add_argument("--max_psi", type=float, default=0.95, help="Maximum average PSI (default: %(default)s)")
    parser.add_argument("--sig_fdr", type=float, default=0.05, help="Threshold for significant FDR (default: %(default)s)")
    parser.add_argument("--bg_fdr", type=float, default=0.5, help="Threshold for background FDR (default: %(default)s)")
    parser.add_argument("--sig_delta_psi", type=float, default=0.1, help="Threshold for significant delta PSI (default: %(default)s)")
    parser.add_argument("--bg_within_group_delta_psi", type=float, default=0.2, help="Threshold for background within-group delta PSI (default: %(default)s)")
    
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    input_path = args.input
    input_fn = os.path.basename(input_path)

    output_files = {
        "filtered": os.path.join(args.outdir, f"filtered_{input_fn}"),
        "upregulated": os.path.join(args.outdir, f"up_{input_fn}"),
        "downregulated": os.path.join(args.outdir, f"dn_{input_fn}"),
        "background": os.path.join(args.outdir, f"bg_{input_fn}")
    }

    header, event_dict = filter_rMATS(
        args.input,
        args.read_cov,
        args.min_psi,
        args.max_psi,
        args.sig_fdr,
        args.bg_fdr,
        args.sig_delta_psi,
        args.bg_within_group_delta_psi
    )

    for category in event_dict:
        if event_dict[category]:
            with open(output_files[category], "w") as f:
                f.write(header)
                for event in event_dict[category]:
                    f.write(str(event))
