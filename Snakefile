import os
import tempfile

import yaml
import pandas as pd

from bioinf_common.tools import execute_notebook

from utils import load_config


###
# setup
#configfile: 'config.yaml'
config = load_config()

# save config for notebooks
conf_fname = tempfile.NamedTemporaryFile().name
os.environ['SNAKEMAKE__CONFIG_FILE'] = conf_fname
with open(conf_fname, 'w') as fd:
    yaml.dump(config, fd)

# create needed output directories
for wd in config['output_dirs'].values():
    os.makedirs(wd, exist_ok=True)

results = config['output_dirs']['results']
images = config['output_dirs']['images']


###
# rule definitions

rule all:
    input:
        f'{results}/TAD_enrichment.csv',
        f'{images}/tad_border_enrichment_corrected.pdf',
        lambda wc: f'{results}/disease_classification.csv' if config['do_further_investigations'] else []

rule convert_tad_coordinates:
    input:
        config['input_files']['tad_coordinates']
    output:
        f'{results}/tads_hg38.tsv'
    run:
        # if config['parameters']['source_genomiccoordinates_version'] == 'hg38':
        #     # skip conversion
        #     df = pd.read_csv(
        #         input[0],
        #         header=None, names=['chrname', 'tad_start', 'tad_stop'])
        #     df.to_csv(output[0], sep='\t', index=False)
        # else:
        #     execute_notebook('ConvertTADGenomicCoordinates.ipynb', inplace=True)

        df = pd.read_csv(
            input[0],
            header=None, names=['chrname', 'tad_start', 'tad_stop'])
        df.to_csv(output[0], sep='\t', index=False)

rule assemble_snp_database:
    input:
        config['input_files']['raw_disgenet'],
        config['input_files']['raw_gwascatalog'],
        f'{results}/tads_hg38.tsv'
    output:
        f'{results}/snpdb_enhanced.tsv',
        f'{results}/disease_cancer_classification.csv',
        f'{results}/disease_efolabels.csv'
    run:
        execute_notebook('EnhanceSNPDatabase.ipynb', inplace=True)

rule compute_enrichments:
    input:
        f'{results}/snpdb_enhanced.tsv',
        f'{results}/tads_hg38.tsv'
    output:
        f'{results}/TAD_enrichment.csv'
    run:
        execute_notebook('ComputeTADEnrichments.ipynb', inplace=True)

rule analyze_results:
    input:
        f'{results}/snpdb_enhanced.tsv',
        f'{results}/TAD_enrichment.csv',
        f'{results}/disease_cancer_classification.csv',
        f'{results}/disease_efolabels.csv',
        f'{results}/tads_hg38.tsv'
    output:
        f'{results}/final.csv',
        f'{images}/tad_border_enrichment_corrected.pdf'
    run:
        execute_notebook('PublicationReproductions.ipynb', inplace=True)

rule further_investigations:
    input:
        f'{results}/snpdb_enhanced.tsv',
        f'{results}/TAD_enrichment.csv',
    output:
        f'{images}/gwas_history.pdf',
        f'{results}/disease_classification.csv'
    run:
        execute_notebook('FurtherExperiments.ipynb', inplace=True)
        execute_notebook('DiseaseClassification.ipynb', inplace=True)
