import sys
from pathlib import Path

import numpy as np
import pandas as pd

import sh
from tqdm import tqdm


def make_axis_regular(df_piv, bin_size, axis):
    """Regularize axis by enforcing a distance of `bin_size` between elements."""
    assert axis in (0, 1), axis

    # get values of axis
    axis_values = (df_piv.index
                   if axis == 0
                   else df_piv.columns)

    # regularize
    idx_missing = np.where(np.diff(axis_values) != bin_size)[0]
    for idx in tqdm(idx_missing, desc=f'Adjusting axis {axis}'):
        start = axis_values[idx]
        end = axis_values[idx+1]

        new_index = np.arange(start + bin_size, end, bin_size)
        for new_idx in new_index:
            if axis == 0:
                df_piv.loc[new_idx, :] = 0
            elif axis == 1:
                df_piv.loc[:, new_idx] = 0

        print(idx, '->', new_index)

    # sort dataframe because new indices are added to the end
    return (df_piv.sort_index()
            if axis == 0
            else df_piv.T.sort_index().T)


def main(fname_in, chrom, fname_matrix, fname_juicer, zero_padding=False):
    bin_size = 10_000

    juicer_exec = Path(snakemake.scriptdir) / '..' / snakemake.config['tool_paths']['juicer_tools']

    # hic -> list
    print('Hi-C -> list')
    sh.java(
        '-jar', juicer_exec,
        'dump', 'observed', 'NONE',
        fname_in,
        chrom, chrom, 'BP', bin_size,
        fname_juicer,
        _out=sys.stdout, _err=sys.stderr)

    # list -> matrix
    print('list -> matrix')
    df = pd.read_csv(
        fname_juicer, sep='\t',
        header=None, names=['start', 'end', 'value'])

    df_piv = pd.pivot(df, index='start', columns='end')

    df_piv.columns = df_piv.columns.droplevel(0)  # drop column label
    # df_piv.columns.name = None
    # df_piv.index.name = None

    piv_shape1 = df_piv.shape

    print('Normalize matrix axes')
    df_piv = make_axis_regular(df_piv, bin_size, 0)
    df_piv = make_axis_regular(df_piv, bin_size, 1)

    piv_shape2 = df_piv.shape

    assert (np.diff(df_piv.index) == bin_size).all()
    assert (np.diff(df_piv.columns) == bin_size).all()

    # enforce same index/column number
    print('Enforce same index/column number')
    only_in_columns = set(df_piv.columns) - set(df_piv.index)
    only_in_index = set(df_piv.index) - set(df_piv.columns)
    print(only_in_columns, only_in_index)

    for idx in tqdm(only_in_columns, desc='Adding to index'):
        df_piv.loc[idx, :] = 0
    for idx in tqdm(only_in_index, desc='Adding to columns'):
        df_piv.loc[:, idx] = 0

    df_piv = df_piv.sort_index()
    df_piv = df_piv.T.sort_index().T

    print(f'Pivot shape: {piv_shape1} -> {piv_shape2} -> {df_piv.shape}')

    # copy upper to lower triangle
    print('Symmetrize matrix')
    mat = df_piv.values.copy()
    i_lower = np.tril_indices(mat.shape[0], -1)
    mat[i_lower] = mat.T[i_lower]

    df_final = pd.DataFrame(
        np.nan_to_num(mat), index=df_piv.index, columns=df_piv.columns)

    # pad with zeros
    start_coord = df_final.index[0]
    if zero_padding and start_coord > 0:
        print('Pad with zeros')
        extra_coords = np.arange(0, start_coord, bin_size)

        empty_counts_tl = np.zeros(
            shape=(len(extra_coords), len(extra_coords)))
        empty_counts_tr = np.zeros(
            shape=(len(extra_coords), df_final.shape[1]))
        empty_counts_bl = np.zeros(
            shape=(df_final.shape[0], len(extra_coords)))

        new_values = np.block([
            [empty_counts_tl, empty_counts_tr],
            [empty_counts_bl, df_final.values]
        ]).astype(int)
        new_index = np.r_[extra_coords, df_final.index]

        tmp = pd.DataFrame(new_values, index=new_index, columns=new_index)
        assert (tmp.loc[start_coord:, start_coord:].values == df_final.values).all()
        df_final = tmp

    # save result
    df_final = df_final.astype(int)
    df_final.to_csv(fname_matrix)


if __name__ == '__main__':
    main(
        snakemake.input.fname,
        snakemake.wildcards.chromosome,
        snakemake.output.fname_matrix, snakemake.output.fname_juicer)
