import pandas as pd


def main(fname_list, fname_out):
    # defined in script because 'run' directive does not work with conda envs
    pd.concat([pd.read_csv(x) for x in fname_list]).to_csv(fname_out, index=False)


if __name__ == '__main__':
    main(snakemake.input.fname_list, snakemake.output.fname)
