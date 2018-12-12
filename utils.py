import os

import yaml

import pandas as pd


def load_config(fname='config.yaml'):
    fname = os.environ.get('SNAKEMAKE__CONFIG_FILE', None) or fname

    with open(fname) as fd:
        config = yaml.load(fd)

    # add directories to config
    required_dirs = ['images', 'cache', 'results']
    for dir_ in required_dirs:
        config['output_dirs'][dir_] = os.path.join(
            config['output_dirs']['prefix'], dir_)

    # create needed directories
    for dir_ in config['output_dirs'].values():
        os.makedirs(dir_, exist_ok=True)

    return config


def split_df_row(df, column, sep=','):
    """
    split_df_row(pd.DataFrame({
        'mappi': ['a,b', 'c', 'd,e'],
        'ha': ['hui', 'bui', 'a'],
        'foo': [1,2,3]
    }), 'mappi', ',')
    """
    def splitter(row, row_accumulator, column):
        split_row = row[column].split(sep)
        for s in split_row:
            new_row = row.to_dict()
            new_row[column] = s.strip()
            row_accumulator.append(new_row)

    new_rows = []
    df.apply(splitter, axis=1, args=(new_rows, column))
    return pd.DataFrame(new_rows, columns=df.columns)
