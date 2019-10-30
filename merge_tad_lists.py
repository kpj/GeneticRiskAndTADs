import os
from pathlib import Path

import pandas as pd
import pyranges as pr

from tqdm import tqdm
from natsort import natsorted


def read_tads(fname):
    return pr.PyRanges(pd.read_csv(
        fname, header=None, names=['Chromosome', 'Start', 'End']))


def simple_aggregation(df1, df2, df3):
    """Retain TAD if it occurs in all lists."""
    tmp = []
    for row in tqdm(df1.itertuples(), total=df1.shape[0], desc='TADs'):
        match2 = df2[
            (df2['Chromosome'] == row.Chromosome) &
            (df2['Start'] == row.Start) &
            (df2['End'] == row.End)]
        match3 = df3[
            (df3['Chromosome'] == row.Chromosome) &
            (df3['Start'] == row.Start) &
            (df3['End'] == row.End)]

        if len(match2) > 0 and len(match3) > 0:
            tmp.append(row)

    return pr.PyRanges(pd.DataFrame(tmp).drop('Index', axis=1))


def complex_aggregation(gr1, gr2, gr3):
    """Find consensus TADs.

    Algorithm:
    1) Find earliest TAD -> anchor TAD
    2) Find all overlaps with anchor TAD -> TAD cloud
    3) Generate consensus TAD: highest start position to
       lowest end position of TAD cloud
    4) Find new anchor TAD: lowest start position above highest
       end position of previous TAD cloud
    """
    # merge sources
    gr1.Origin = 'k1'
    gr2.Origin = 'k2'
    gr3.Origin = 'k3'

    gr_all = pr.concat([gr1, gr2, gr3]).sort().drop_duplicate_positions()

    # agregate TADs
    tmp = []

    for chrom, group in tqdm(
        gr_all, total=len(gr_all.chromosomes),
        desc='Chromosomes'
    ):
        gr_cur = pr.PyRanges(group)

        pbar = tqdm(total=len(gr_cur), desc='TADs')
        while not gr_cur.empty:
            # find anchor TAD
            anchor = pr.PyRanges(pd.DataFrame([
                gr_cur.df.loc[gr_cur.df['Start'].idxmin()]
            ]))
            # print(anchor)

            # find TAD cloud
            cloud = gr_cur.overlap(anchor)
            # print(cloud)

            # generate consensus
            chr_consensus = cloud.Chromosome.iloc[0]
            start_consensus = cloud.Start.max()
            end_consensus = cloud.End.min()

            assert chrom == chr_consensus
            tmp.append({
                'Chromosome': chr_consensus,
                'Start': start_consensus,
                'End': end_consensus
            })
            # print(tmp[-1])
            # print()

            # remove processed TADs
            prev_len = len(gr_cur)
            gr_cur = gr_cur[gr_cur.Start > cloud.End.max()]

            pbar.update(prev_len - len(gr_cur))
            # input()
        pbar.close()

    return pr.PyRanges(pd.DataFrame(tmp))


def main(fname1, fname2, fname3, outfile_format):
    # read input
    gr1 = read_tads(fname1)
    gr2 = read_tads(fname2)
    gr3 = read_tads(fname3)

    # generate consensus
    gr_simple = simple_aggregation(gr1.df, gr2.df, gr3.df)
    gr_complex = complex_aggregation(gr1, gr2, gr3)

    # save consensus
    gr_simple.to_csv(
        outfile_format.format(consensus_type='simple'), header=False)
    gr_complex.to_csv(
        outfile_format.format(consensus_type='complex'), header=False)


if __name__ == '__main__':
    data_source = 'newleopoldtads_Rao_IMR90_40k_hg38'

    root = Path('data/')
    outdir = root / 'generated'
    outdir.mkdir(exist_ok=True)

    file_list = natsorted(str(x) for x in root.glob(f'{data_source}*'))

    file_combinations = [(file_list[i], file_list[i+1], file_list[i+2])
                         for i in range(len(file_list)-2)]
    for fname1, fname2, fname3 in tqdm(
        file_combinations,
        desc='File combinations'
    ):
        k1 = os.path.basename(fname1).split('.')[0].split('_')[-1]
        k2 = os.path.basename(fname2).split('.')[0].split('_')[-1]
        k3 = os.path.basename(fname3).split('.')[0].split('_')[-1]

        prefix = os.path.basename(os.path.commonprefix(file_list))
        outfile_format = str(outdir / (prefix + f'{k1},{k2},{k3}_{{consensus_type}}.csv'))

        tqdm.write(outfile_format)
        main(fname1, fname2, fname3, outfile_format)
