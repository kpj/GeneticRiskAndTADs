"""
Functions for reading various data sources
"""

import pandas as pd
import networkx as nx


def load_gwas(
    fname: str='data/gwas_catalog_v1.0.1-associations_e90_r2017-11-13.tsv'
) -> pd.DataFrame:
    df = pd.read_table(fname, low_memory=False)
    return df

def load_snp_data(snp_dir: str = 'data/snps') -> pd.DataFrame:
    kwargs = dict(header=None, names=[
        'EFO_term', 'disease_name', 'SNP_name', 'chromosome', 'position'])
    df_cancer_tad = pd.read_table(f'{snp_dir}/dSNPs_cancer_tad.dat', **kwargs)
    df_cancer_notad = pd.read_table(f'{snp_dir}/dSNPs_cancer.dat', **kwargs)
    df_nocancer_tad = pd.read_table(f'{snp_dir}/dSNPs_noncancer_tad.dat', **kwargs)
    df_nocancer_notad = pd.read_table(f'{snp_dir}/dSNPs_noncancer.dat', **kwargs)

    df_cancer_tad['is_cancer'] = True
    df_cancer_tad['is_tad'] = True
    df_cancer_notad['is_cancer'] = True
    df_cancer_notad['is_tad'] = False
    df_nocancer_tad['is_cancer'] = False
    df_nocancer_tad['is_tad'] = True
    df_nocancer_notad['is_cancer'] = False
    df_nocancer_notad['is_tad'] = False

    return pd.concat([
        df_cancer_tad, df_cancer_notad,
        df_nocancer_tad, df_nocancer_notad], axis=0)

def load_biogrid(
    fname: str = 'data/BIOGRID-ORGANISM-Homo_sapiens-3.4.153.tab2.txt'
) -> nx.Graph:
    """ Load protein-interaction data from BioGRID
    """
    # BioGRID
    df_bgrid = pd.read_csv(fname, sep='\t', low_memory=False)

    # subset physical interactions
    df_bgrid_phys = df_bgrid[df_bgrid['Experimental System Type'] == 'physical']

    # create network
    graph_bgrid_all = nx.convert_matrix.from_pandas_edgelist(
        df_bgrid_phys,
        source='Entrez Gene Interactor A', target='Entrez Gene Interactor B',
        edge_attr='Score')
    graph_bgrid = max(nx.connected_component_subgraphs(graph_bgrid_all), key=len)

    return graph_bgrid, df_bgrid_phys
