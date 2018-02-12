"""
Functions for reading various data sources
"""

import os

from typing import Tuple, Dict

import pandas as pd
import networkx as nx


def load_gwas(
    fname: str='data/gwas_catalog_v1.0.1-associations_e90_r2017-11-13.tsv.gz'
) -> pd.DataFrame:
    df = pd.read_table(fname, low_memory=False)
    return df

def load_snp_data(snp_dir: str = 'data/snps') -> pd.DataFrame:
    kwargs = dict(header=None, names=[
        'EFO_term', 'disease_name', 'SNP_name', 'chromosome', 'position'])

    # load all SNP-disease associations
    df_cancer = pd.read_table(
        f'{snp_dir}/all/dSNPs_cancer2.dat', **kwargs)
    df_nocancer = pd.read_table(
        f'{snp_dir}/all/dSNPs_noncancer2.dat', **kwargs)

    df_cancer['is_cancer'] = True
    df_nocancer['is_cancer'] = False

    df_all = pd.concat([df_cancer, df_nocancer], axis=0)

    # mark SNPs in TAD-borders
    df_nocancer_tad = pd.read_table(
        f'{snp_dir}/all/dSNPs_noncancer_tad2.dat', **kwargs)
    df_cancer_tad = pd.read_table(
        f'{snp_dir}/all/dSNPs_cancer_tad2.dat', **kwargs)
    df_all_tad = pd.concat([df_nocancer_tad, df_cancer_tad], axis=0)

    snps_in_tad = set(df_all_tad.SNP_name.unique())
    df_all['is_tad'] = df_all['SNP_name'].apply(
        lambda x: x in snps_in_tad)

    # mark diseases with TAD border enrichment
    df_enr = pd.read_table(
        f'{snp_dir}/tad_enr/dSNPs_all.dat', **kwargs)

    enr_efos = set(df_enr.EFO_term.unique())
    df_all['disease_tad_enriched'] = df_all['EFO_term'].apply(
        lambda x: x in enr_efos)

    return df_all


def main():
    df = load_snp_data()
    import ipdb; ipdb.set_trace()

if __name__ == '__main__':
    main()
