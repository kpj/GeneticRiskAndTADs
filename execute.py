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
    configurations = [{'index': 'raw'}]

    # other GWAS versions
    configurations.append({
        'index': 'oldgwas:e90',
        'input_files': dict(
            raw_gwascatalog='data/gwas_catalog_v1.0.1-associations_e90_r2017-11-20.tsv')
    })

    # OR thresholds
    for or_thres in [1.3, 1.5, 2]:
        configurations.append({
            'index': f'ORthres_{or_thres}',
            'filters': dict(OR_threshold=or_thres)
        })

    # variant filters
    for var_filter in ['nonexonic', 'intergenic']:
        configurations.append({
            'index': var_filter,
            'filters': dict(variant_type=var_filter)
        })

    # other TAD data
    for tad_idx, tad_fname in [
        ('tads_IMR90', 'data/tads_IMR90_hg18.csv'),
        ('tads_cortex', 'data/tads_cortex_hg18.csv'),
        ('tads_mESC', 'data/tads_mESC_hg18.csv'),
        ('tads_random', 'data/tads_hg18_randomized.tsv'),
        ('DomainCaller_500M_50000', 'data/TADcallsByTool_Rao_DomainCaller_500M_50000.tsv'),
        ('HiCSeg_500M_50000', 'data/TADcallsByTool_Rao_HiCSeg_500M_50000.tsv'),
        ('TADbit_500M_50000', 'data/TADcallsByTool_Rao_TADbit_500M_50000.tsv'),
        ('TADtree_500M_50000', 'data/TADcallsByTool_Rao_TADtree_500M_50000.tsv'),
        ('TopDom_500M_50000', 'data/TADcallsByTool_Rao_TopDom_500M_50000.tsv'),
        ('armatus_500M_50000', 'data/TADcallsByTool_Rao_armatus_500M_50000.tsv'),
        ('arrowhead_500M_50000', 'data/TADcallsByTool_Rao_arrowhead_500M_50000.tsv')
    ]:
        configurations.append({
            'index': tad_idx,
            'input_files': dict(
                tad_coordinates_hg18=tad_fname)
        })

    configurations.append({
        'index': 'tads_hESChg19',
        'input_files': dict(
            tad_coordinates_hg18='data/tads_hESC_hg19.csv'),
        'parameters': dict(
            source_genomiccoordinates_version='hg19')
    })

    # execute pipelines
    for conf in tqdm(configurations):
        idx = conf.pop('index')

        print(idx, end=' ')
        if os.path.exists(f'{execution_dir}/{idx}/'):
            print('-- Skipping...')
            continue
        print()

        # merge default and run-specific configurations
        cur_conf = update(copy.deepcopy(default_config), conf)

        # setup environment
        for key, dir_value in default_config['output_dirs'].items():
            cwd = f'{execution_dir}/{idx}/{dir_value}'

            cur_conf['output_dirs'][key] = cwd
            sh.mkdir('-p', cwd)

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

        for res_type in os.scandir(conf_type.path):
            for res in os.scandir(f'{conf_type.path}/{res_type.name}'):
                suf = res.name.split('.')[-1]
                res_name = res.name[:-(len(suf)+1)] + f'_{c_type}.{suf}'

                sh.mkdir('-p', f'{results_dir}/{res_type.name}')
                sh.cp(res.path, f'{results_dir}/{res_type.name}/{res_name}')

if __name__ == '__main__':
    main()
