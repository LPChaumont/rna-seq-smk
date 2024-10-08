import pandas as pd
import sys

sys.stderr = open(snakemake.log[0], "w")

for i, file in enumerate(sorted(snakemake.input.quants)):
    sample = file.split('/')[-1].replace('.tsv', '')
    df = pd.read_table(file)

    if i == 0:
        counts_df = df[['gene_id', 'gene_name', 'count']].rename(columns={'count': sample})
        cpm_df = df[['gene_id', 'gene_name', 'cpm']].rename(columns={'cpm': sample})
        tpm_df = df[['gene_id', 'gene_name', 'tpm']].rename(columns={'tpm': sample})

    else:
        counts_df[sample] = counts_df['gene_id'].map(dict(zip(df.gene_id, df['count'])))
        cpm_df[sample] = cpm_df['gene_id'].map(dict(zip(df.gene_id, df.cpm)))
        tpm_df[sample] = tpm_df['gene_id'].map(dict(zip(df.gene_id, df.tpm)))

counts_df.to_csv(snakemake.output.counts, sep='\t', index=False)
cpm_df.to_csv(snakemake.output.cpm, sep='\t', index=False)
tpm_df.to_csv(snakemake.output.tpm, sep='\t', index=False)
