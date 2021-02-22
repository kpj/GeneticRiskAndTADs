import os

import pandas as pd

import cooler


def main(fname_list, samplename_list, fname_out):
    tmp = []

    for samplename, fname in zip(samplename_list, fname_list):
        c = cooler.Cooler(fname)

        tmp.append(
            {
                'filename': fname,
                'source': samplename,
                'genome_assembly': c.info['genome-assembly'],
                'bin_size': c.binsize,
            }
        )

    pd.DataFrame(tmp).to_csv(fname_out, index=False)


if __name__ == '__main__':
    main(
        snakemake.input.fname_list,
        snakemake.params['samplename_list'],
        snakemake.output.fname,
    )
