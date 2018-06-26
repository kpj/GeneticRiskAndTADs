"""
Execute various pipeline configurations
"""

import os
import sys
import copy
import collections
from pprint import pprint

import sh
from tqdm import tqdm

from utils import load_config


def update(d, u):
    for k, v in u.items():
        if isinstance(v, collections.Mapping):
            d[k] = update(d.get(k, {}), v)
        else:
            d[k] = v
    return d


def main(execution_dir='cwd_multiconfig_execution', results_dir='toshow'):
    # load default config
    default_config = load_config()

    # setup configurations
    configurations = [
        {
            'index': 'raw'
        },
        {
            'index': 'oldgwas:e90',
            'input_files': dict(
                raw_gwascatalog='data/gwas_catalog_v1.0.1-associations_e90_r2017-11-20.tsv')
        },
        {
            'index': 'ORthres_1,3',
            'filters': dict(OR_threshold=1.3)
        },
        {
            'index': 'ORthres_1,5',
            'filters': dict(OR_threshold=1.5)
        },
        {
            'index': 'ORthres_2',
            'filters': dict(OR_threshold=2)
        },
        {
            'index': 'nonexonic',
            'filters': dict(variant_type='nonexonic')
        },
        {
            'index': 'intergenic',
            'filters': dict(variant_type='intergenic')
        },
        {
            'index': 'tads:random',
            'input_files': dict(
                tad_coordinates_hg19='data/tads_hg19_randomized.tsv')
        },
        {
            'index': 'tads:new:rao',
            'input_files': dict(
                tad_coordinates_hg19='data/TADcallsByTool_Rao.tsv')
        },
    ]

    # execute pipelines
    for conf in tqdm(configurations):
        idx = conf.pop('index')

        # merge default and run-specific configurations
        cur_conf = update(copy.deepcopy(default_config), conf)

        # setup environment
        for key, dir_value in default_config['output_dirs'].items():
            cwd = f'{execution_dir}/{idx}/{dir_value}'

            cur_conf['output_dirs'][key] = cwd
            sh.mkdir('-p', cwd)

        print(idx)
        pprint(cur_conf)
        print()

        # run
        c_list = [f'{k}={v}' for k, v in cur_conf.items()]
        sh.snakemake(
            '-pr',
            '--config', *c_list,
            _out=sys.stdout, _err=sys.stderr)

    # gather results
    sh.rm('-rf', results_dir)
    sh.mkdir('-p', results_dir)
    for conf_type in os.scandir(execution_dir):
        c_type = conf_type.name
        for img in os.scandir(f'{conf_type.path}/images'):
            img_name = img.name.replace('.pdf', f'_{c_type}.pdf')
            sh.cp(img.path, f'{results_dir}/{img_name}')

if __name__ == '__main__':
    main()
