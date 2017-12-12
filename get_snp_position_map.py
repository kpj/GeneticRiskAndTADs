"""
Given a list of SNP IDs retrieve their position from dbSNP
"""

import os
import urllib
import tempfile

import pandas as pd


def ensure_dir(path):
    dname = os.path.dirname(path)
    os.makedirs(dname, exist_ok=True)
    return path

def download_bed(chr_id, fname):
    _url = 'ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/' \
           f'human_9606_b149_GRCh37p13/BED/bed_chr_{chr_id}.bed.gz'
    urllib.request.urlretrieve(_url, fname)

def get_snp_list():
    df = pd.read_table('data/all_variant_disease_pmid_associations.tsv.gz')
    return df['snpId'].dropna().unique().tolist()

def main(dbSNP_dir):
    # save SNPs to file
    snp_list = get_snp_list()

    snplist_fname = tempfile.NamedTemporaryFile().name
    print('SNPs:', snplist_fname)
    with open(snplist_fname, 'w') as fd:
        fd.write('\n'.join(snp_list))

    # download chromosome BEDs if needed
    ensure_dir(dbSNP_dir)
    chroms = list(range(1, 23)) + ['X', 'Y']
    for chr_id in chroms:
        fname = os.path.join(dbSNP_dir, f'GRCh37p13_bed_chr_{chr_id}.bed.gz')
        if not os.path.exists(fname):
            print(
                f'Downloading "chr_{chr_id}" to "{fname}"...',
                end=' ', flush=True)
            try:
                download_bed(chr_id, fname)
                print('Done')
            except urllib.error.URLError:
                print('Fail')

    # subset to needed SNPs
    for entry in os.scandir(dbSNP_dir):
        if entry.name.startswith('match_'):
            continue

        pm_name = '.'.join(entry.name.split('.')[:-1])
        outfile = os.path.join(dbSNP_dir, f'match_{pm_name}')

        if not os.path.exists(f'{outfile}.gz'):
            cmd = f'zgrep -wFf {snplist_fname} {entry.path} > {outfile}'
            print(f'Running "{cmd}"...', end=' ', flush=True)
            os.system(cmd)
            os.system(f'gzip {outfile}')
            print('Done')

    # concat SNP positions
    df_map_list = []
    for entry in os.scandir(dbSNP_dir):
        if not entry.name.startswith('match_'):
            continue

        cur = pd.read_table(
            entry.path,
            header=None,
            names=('chrom', 'start', 'end', 'SNPS', 'foo', 'strand'))
        df_map_list.append(cur)
    df_map = pd.concat(df_map_list)

    # save result
    df_map.to_csv('snp_positions_hg19.csv', index=False)

if __name__ == '__main__':
    main('./cache/dbSNP/')
