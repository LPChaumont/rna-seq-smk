import pandas as pd
import os
import argparse


def filter_gene_expression(file_path, condition_dict, min_gene_expr=1, min_condition_expr=1, save=False):
    df = pd.read_table(file_path)
    non_numeric_cols = df.select_dtypes(exclude='number').columns.tolist()

    melted_df = (
        df
        .melt(id_vars=non_numeric_cols, var_name='sample', value_name='abundance')
        .assign(condition=lambda x: x['sample'].map(condition_dict))
        .groupby(non_numeric_cols + ['condition'])['abundance'] 
        .mean()  
        .reset_index()  
    )

    filtered_index = (
        melted_df[melted_df['abundance'].gt(min_gene_expr)]
        .groupby(non_numeric_cols)['condition']
        .nunique()
        .loc[lambda x: x >= min_condition_expr]
        .index
    )

    filtered_df = (
        df.set_index(non_numeric_cols)
        .loc[lambda x: x.index.isin(filtered_index)]
        .reset_index()
    )

    mean_df = (
        melted_df
        .pivot(index=non_numeric_cols, columns='condition', values='abundance')
        .reset_index()
        .rename_axis(None, axis=1)
    )

    filtered_mean_df = (
        mean_df.set_index(non_numeric_cols)
        .loc[lambda x: x.index.isin(filtered_index)]
        .reset_index()
    )

    if save:
        basename = os.path.basename(file_path)
        dirname = os.path.dirname(file_path)
        filtered_fout = os.path.join(dirname, 'filtered_' + basename)
        mean_fout = os.path.join(dirname, 'mean_' + basename)
        filtered_mean_fout = os.path.join(dirname, 'filtered_mean_' + basename)

        filtered_df.to_csv(filtered_fout, index=None, sep='\t')
        mean_df.to_csv(mean_fout, index=None, sep='\t')
        filtered_mean_df.to_csv(filtered_mean_fout, index=None, sep='\t')

    return filtered_df, mean_df, filtered_mean_df


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Filter gene expression data'
    )
    parser.add_argument(
        '-i', '--input',
        type=str,
        required=True,
        help='path to the gene expression file (e.g. CoCo TPM file)'
    )
    parser.add_argument(
        '-s','--samples',
        type=str,
        required=True,
        help='path to the sample sheet which includes metadata about samples and conditions'
    )
    parser.add_argument(
        '--min_gene_expr',
        type=float,
        default=1,
        help='minimum gene expression threshold (default: %(default)s)'
    )
    parser.add_argument(
        '--min_condition_expr',
        type=int,
        default=1,
        help='minimum number of conditions a gene must be expressed in (default: %(default)s)'
    )
    parser.add_argument(
        '--save',
        action='store_true',
        help='saves the filtered dataframes at the same location as the input file'
    )

    args = parser.parse_args()

    samples_df = pd.read_table(args.samples)
    condition_dict = dict(zip(samples_df['sample'], samples_df['condition']))

    filter_gene_expression(
        file_path=args.input,
        condition_dict=condition_dict,
        min_gene_expr=args.min_gene_expr,
        min_condition_expr=args.min_condition_expr,
        save=args.save
    )
