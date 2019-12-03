import numpy as np
import pandas as pd

import sh


def main():
    bin_size = 10_000

    # read data
    df_count = pd.read_csv(snakemake.input.fname, index_col=0)

    # convert to TopDom compatible format
    df_count.insert(0, 'chr', f'chr{snakemake.wildcards.chromosome}')
    df_count.insert(1, 'td_start', df_count.index)
    df_count.insert(2, 'td_end', df_count.index + bin_size)

    df_count.to_csv(
        snakemake.output.topdom_input, sep='\t',
        index=False, header=False)

    # run TopDom
    cmd = """
        TopDom::TopDom('{input}', 10, outFile='{output}', debug=TRUE)
    """.format(
        input=snakemake.output.topdom_input,
        output=snakemake.params.prefix)
    sh.Rscript('-e', cmd, _fg=True)

    # extract TADs
    df_topdom = pd.read_csv(
        snakemake.output.topdom_output, sep='\t',
        header=None, names=['chrname', 'tad_start', 'tad_stop', 'type'])

    df_topdom = df_topdom[df_topdom['type'] == 'domain']

    # rescale coordinates
    df_topdom['tad_start'] /= bin_size
    df_topdom['tad_stop'] /= bin_size

    # save result
    df_topdom.drop('type', axis=1, inplace=True)

    df_topdom['tad_start'] = df_topdom['tad_start'].astype(int)
    df_topdom['tad_stop'] = df_topdom['tad_stop'].astype(int)

    df_topdom.to_csv(snakemake.output.fname, index=False)


if __name__ == '__main__':
    main()
