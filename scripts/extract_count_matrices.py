import sys
from pathlib import Path

import numpy as np
import pandas as pd

import sh
from tqdm import tqdm


def main(fname_in, chrom_list, fname_matrix, fname_juicer):
    binsize = 10_000

    juicer_exec = Path(snakemake.scriptdir) / '..' / snakemake.config['tool_paths']['juicer_tools']

    for chrom in tqdm(chrom_list, desc='Chromosome'):
        # hic -> list
        print('Hi-C -> list')
        sh.java(
            '-jar', juicer_exec,
            'dump', 'observed', 'NONE',
            fname_in,
            chrom, chrom, 'BP', binsize,
            fname_juicer,
            _out=sys.stdout, _err=sys.stderr)

        # list -> matrix
        print('list -> matrix')
        df = pd.read_csv(
            fname_juicer, sep='\t',
            header=None, names=['start', 'end', 'value'])

        df_piv = pd.pivot(df, index='start', columns='end')

        print('Symmetrize matrix')
        only_in_columns = set([x[1] for x in df_piv.columns]) - set(df_piv.index.tolist())
        df_sym = df_piv.drop([('value', x) for x in only_in_columns], axis=1)

        only_in_index = set(df_piv.index.tolist()) - set([x[1] for x in df_piv.columns])
        df_sym = df_sym.drop(only_in_index, axis=0)

        mat = df_sym.values.copy()
        i_lower = np.tril_indices(mat.shape[0], -1)
        mat[i_lower] = mat.T[i_lower]

        df_final = pd.DataFrame(
            np.nan_to_num(mat), index=df_sym.index, columns=df_sym.columns)

        df_final['value'] = df_final['value'].astype(int)
        df_final = df_final['value']  # handle multilevel columns

        df_final.to_csv(fname_matrix)


if __name__ == '__main__':
    main(
        snakemake.input.fname,
        snakemake.config['chromosome_list'],
        snakemake.output.fname_matrix, snakemake.output.fname_juicer)
