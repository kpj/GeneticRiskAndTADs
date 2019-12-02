import os
import urllib.request

import GEOparse

from tqdm import tqdm


def main(hic_source_list):
    gse_list = {x.split('_')[0] for x in hic_source_list}
    print(f'GSE selection: {gse_list}')

    for geo_idx in tqdm(gse_list):
        gse = GEOparse.get_GEO(
            geo=geo_idx,
            destdir='./hic_files/raw/', silent=True)
        for url in tqdm(gse.metadata['supplementary_file']):
            if not url.endswith('combined.hic'):
                continue

            idx = url.split('/')[-1][:-len('_combined.hic')]
            if idx not in hic_source_list:
                print(f'Skipping {idx}')
                continue

            print(f'Downloading {idx}')
            urllib.request.urlretrieve(url, f'hic_files/raw/data.{idx}.hic')


if __name__ == '__main__':
    main(snakemake.config['hic_sources'])
