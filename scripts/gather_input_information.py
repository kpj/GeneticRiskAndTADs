import os

import pandas as pd

import cooler


def main(fname_list, fname_out):
    tmp = []

    for fname in fname_list:
        c = cooler.Cooler(fname)

        tmp.append({
            'filename': fname,
            'source': os.path.splitext(os.path.basename(fname))[0],
            'genome_assembly': c.info['genome-assembly'],
            'bin_size': c.binsize
        })

    pd.DataFrame(tmp).to_csv(fname_out, index=False)


if __name__ == '__main__':
    main(snakemake.input.fname_list, snakemake.output.fname)
