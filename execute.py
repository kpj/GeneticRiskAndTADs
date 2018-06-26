"""
Execute various pipeline configurations
"""

import os
import sys
import copy
import collections

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


def main():
    # load default config
    default_config = load_config()

    # setup configurations
    ow = 'exec_output'
    configurations = [
        {
            'output_dirs': dict(
                results=f'{ow}/raw/results',
                images=f'{ow}/raw/images')
        },
        {
            'output_dirs': dict(
                results=f'{ow}/oldgwas:e90/results',
                images=f'{ow}/oldgwas:e90/images'),
            'input_files': dict(
                raw_gwascatalog='data/gwas_catalog_v1.0.1-associations_e90_r2017-11-20.tsv')
        },
        {
            'output_dirs': dict(
                results=f'{ow}/tads:random/results',
                images=f'{ow}/tads:random/images'),
            'input_files': dict(
                tad_coordinates_hg19='data/tads_hg19_randomized.tsv')
        },
        {
            'output_dirs': dict(
                results=f'{ow}/tads:new:rao/results',
                images=f'{ow}/tads:new:rao/images'),
            'input_files': dict(
                tad_coordinates_hg19='data/TADcallsByTool_Rao.tsv')
        },
        {
            'output_dirs': dict(
                results=f'{ow}/ORthres_1,3/results',
                images=f'{ow}/ORthres_1,3/images'),
            'filters': dict(OR_threshold=1.3)
        },
        {
            'output_dirs': dict(
                results=f'{ow}/ORthres_1,5/results',
                images=f'{ow}/ORthres_1,5/images'),
            'filters': dict(OR_threshold=1.5)
        },
        {
            'output_dirs': dict(
                results=f'{ow}/ORthres_2/results',
                images=f'{ow}/ORthres_2/images'),
            'filters': dict(OR_threshold=2)
        },
        {
            'output_dirs': dict(
                results=f'{ow}/nonexonic/results',
                images=f'{ow}/nonexonic/images'),
            'filters': dict(variant_type='nonexonic')
        },
        {
            'output_dirs': dict(
                results=f'{ow}/intergenic/results',
                images=f'{ow}/intergenic/images'),
            'filters': dict(variant_type='intergenic')
        }
    ]

    # execute pipelines
    for conf in tqdm(configurations):
        cur_conf = update(copy.deepcopy(default_config), conf)

        for k, v in cur_conf['output_dirs'].items():
            sh.mkdir('-p', v)

        c_list = [f'{k}={v}' for k, v in cur_conf.items()]
        sh.snakemake(
            '-pr',
            '--config', *c_list,
            _out=sys.stdout, _err=sys.stderr)

    # gather results
    output_dir = 'toshow'
    sh.rm('-rf', output_dir)
    sh.mkdir('-p', output_dir)
    for conf_type in os.scandir(ow):
        c_type = conf_type.name
        for img in os.scandir(f'{conf_type.path}/images'):
            img_name = img.name.replace('.pdf', f'_{c_type}.pdf')
            sh.cp(img.path, f'{output_dir}/{img_name}')

if __name__ == '__main__':
    main()
