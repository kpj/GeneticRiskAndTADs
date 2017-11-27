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

    enr_snps = set(df_enr.SNP_name.unique())
    df_all['disease_tad_enriched'] = df_all['SNP_name'].apply(
        lambda x: x in enr_snps)

    return df_all

def load_biogrid(
    fname: str = 'data/BIOGRID-ORGANISM-Homo_sapiens-3.4.153.tab2.txt'
) -> Tuple[nx.Graph, pd.DataFrame]:
    """ Load protein-interaction data from BioGRID
    """
    # BioGRID
    df_bgrid = pd.read_table(fname, low_memory=False)

    # subset physical interactions
    df_bgrid_phys = df_bgrid[df_bgrid['Experimental System Type'] == 'physical']

    # create network
    graph_bgrid_all = nx.convert_matrix.from_pandas_edgelist(
        df_bgrid_phys,
        source='Entrez Gene Interactor A', target='Entrez Gene Interactor B',
        edge_attr='Score')
    graph_bgrid = max(nx.connected_component_subgraphs(graph_bgrid_all), key=len)

    return graph_bgrid, df_bgrid_phys

def load_stringdb(
    fname: str = 'data/stringdb_entrez.tsv.gz',
    threshold: int = 850
) -> Tuple[nx.Graph, pd.DataFrame]:
    """ Load protein-interaction data from StringDB and convert entries to ENTREZ ids
    """
    if not os.path.exists(fname):
        print('Please execute LoadStringDB.ipynb')
        exit(-1)

    df = pd.read_table(fname)
    df_sub = df[df['combined_score']>threshold]

    graph_all = nx.convert_matrix.from_pandas_edgelist(
        df_sub,
        source='protein1', target='protein2',
        edge_attr='combined_score')
    graph = max(nx.connected_component_subgraphs(graph_all), key=len)

    return graph, df_sub


def main():
    #df = load_snp_data()
    ppi, df = load_stringdb()
    import ipdb; ipdb.set_trace()

if __name__ == '__main__':
    main()
