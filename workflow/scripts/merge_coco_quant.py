import pandas as pd
import sys

sys.stderr = open(snakemake.log[0], "w")

def get_gene_biotype(gtf_path):
    """
    Parses GTF file from Ensembl and returns a dictionnary mapping
    each gene ID to its corresponding gene biotype.
    """
    biotype_dict = {}
    with open(gtf_path, 'r') as f:
        for line in f:
            if 'gene_id' in line and 'gene_biotype' in line:
                gene_id = line.split('gene_id "')[1].split('"', 1)[0]
                gene_biotype = line.split('gene_biotype "')[1].split('"', 1)[0]
                biotype_dict[gene_id] = gene_biotype

    return biotype_dict


biotype_dict = get_gene_biotype(snakemake.input.gtf)

for i, file in enumerate(sorted(snakemake.input.quants)):
    sample = file.split('/')[-1].replace('.tsv', '')
    df = pd.read_table(file)

    # Initialize counts, CPM and TPM dataframes
    if i == 0:
        counts_df = df[['gene_id', 'gene_name', 'count']].rename(columns={'count': sample})
        cpm_df = df[['gene_id', 'gene_name', 'cpm']].rename(columns={'cpm': sample})
        tpm_df = df[['gene_id', 'gene_name', 'tpm']].rename(columns={'tpm': sample})

        # Add gene_biotype column and reorder columns in place
        counts_df['gene_biotype'] = counts_df['gene_id'].map(biotype_dict)
        counts_df = counts_df[['gene_id', 'gene_name', 'gene_biotype', sample]]

        cpm_df['gene_biotype'] = cpm_df['gene_id'].map(biotype_dict)
        cpm_df = cpm_df[['gene_id', 'gene_name', 'gene_biotype', sample]]

        tpm_df['gene_biotype'] = tpm_df['gene_id'].map(biotype_dict)
        tpm_df = tpm_df[['gene_id', 'gene_name', 'gene_biotype', sample]]

    # Fill dataframes
    else:
        counts_df[sample] = counts_df['gene_id'].map(dict(zip(df.gene_id, df['count'])))
        cpm_df[sample] = cpm_df['gene_id'].map(dict(zip(df.gene_id, df.cpm)))
        tpm_df[sample] = tpm_df['gene_id'].map(dict(zip(df.gene_id, df.tpm)))

counts_df.to_csv(snakemake.output.counts, sep='\t', index=False)
cpm_df.to_csv(snakemake.output.cpm, sep='\t', index=False)
tpm_df.to_csv(snakemake.output.tpm, sep='\t', index=False)
