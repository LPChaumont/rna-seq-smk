import pandas as pd
import sys
import os

sys.stderr = open(snakemake.log[0], "w")

df_dict = {"counts": None, "cpm": None, "tpm": None}

for i, file in enumerate(sorted(snakemake.input.quants)):
    sample = os.path.splitext(os.path.basename(file))[0]
    df = pd.read_table(file)

    if i == 0:
        for metric in df_dict.keys():
            df_dict[metric] = df[["gene_id", "gene_name", metric]].rename(
                columns={metric: sample}
            )

    else:
        for metric in df_dict.keys():
            df_dict[metric] = pd.merge(
                df_dict[metric], df[["gene_id", metric]], on="gene_id", how="outer"
            )

for metric, df in df_dict.items():
    df.to_csv(snakemake.output[metric], sep="\t", index=False)
