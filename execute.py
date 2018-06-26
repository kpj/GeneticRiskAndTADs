"""
Execute various pipeline configurations
"""

import os
import sys

import sh
from tqdm import tqdm


def main():
    ow = 'exec_output'
    configurations = [
        {
            'output_dirs': dict(
                results=f'{ow}/raw/results',
                images=f'{ow}/raw/images'),
            'filters': dict(
                OR_threshold=None, variant_type=None)
        },
        {
            'output_dirs': dict(
                results=f'{ow}/ORthres_1,3/results',
                images=f'{ow}/ORthres_1,3/images'),
            'filters': dict(
                OR_threshold=1.3, variant_type=None)
        },
        {
            'output_dirs': dict(
                results=f'{ow}/ORthres_1,5/results',
                images=f'{ow}/ORthres_1,5/images'),
            'filters': dict(
                OR_threshold=1.5, variant_type=None)
        },
        {
            'output_dirs': dict(
                results=f'{ow}/ORthres_2/results',
                images=f'{ow}/ORthres_2/images'),
            'filters': dict(
                OR_threshold=2, variant_type=None)
        },
        {
            'output_dirs': dict(
                results=f'{ow}/nonexonic/results',
                images=f'{ow}/nonexonic/images'),
            'filters': dict(
                OR_threshold=None, variant_type='nonexonic')
        },
        {
            'output_dirs': dict(
                results=f'{ow}/intergenic/results',
                images=f'{ow}/intergenic/images'),
            'filters': dict(
                OR_threshold=None, variant_type='intergenic')
        }
    ]

    # execute pipelines
    for conf in tqdm(configurations):
        c_list = [f'{k}={v}' for k, v in conf.items()]
        print(c_list)
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
