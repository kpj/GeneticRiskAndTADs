import os
import tempfile

import yaml

from bioinf_common.tools import execute_notebook


###
# setup
configfile: 'config.yaml'

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
        f'{images}/tad_border_enrichment.pdf'

rule convert_tad_coordinates:
    input:
        config['input_files']['tad_coordinates_hg18']
    output:
        f'{results}/tads_hg38.tsv'
    run:
        execute_notebook('ConvertTADGenomicCoordinates.ipynb')

rule assemble_snp_database:
    input:
        config['input_files']['raw_disgenet'],
        config['input_files']['raw_gwascatalog'],
        f'{results}/tads_hg38.tsv'
    output:
        f'{results}/disgenet_enhanced.tsv',
        f'{results}/disease_cancer_classification.csv',
        f'{results}/disease_efolabels.csv'
    run:
        execute_notebook('EnhanceSNPDatabase.ipynb')

rule compute_enrichments:
    input:
        f'{results}/disgenet_enhanced.tsv',
        f'{results}/tads_hg38.tsv'
    output:
        f'{results}/TAD_enrichment.csv'
    run:
        execute_notebook('ComputeTADEnrichments.ipynb')

rule analyze_results:
    input:
        f'{results}/disgenet_enhanced.tsv',
        f'{results}/TAD_enrichment.csv',
        f'{results}/disease_cancer_classification.csv',
        f'{results}/disease_efolabels.csv',
        f'{results}/tads_hg38.tsv'
    output:
        f'{images}/tad_border_enrichment.pdf'
    run:
        execute_notebook('PublicationReproductions.ipynb')
        execute_notebook('FurtherExperiments.ipynb')
