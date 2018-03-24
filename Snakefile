import nbformat
from nbconvert.preprocessors import ExecutePreprocessor


###
# setup

def execute_notebook(nb_path):
    with open(nb_path) as fd:
        nb = nbformat.read(fd, as_version=4)

    ep = ExecutePreprocessor(timeout=600)
    ep.preprocess(nb, {})


###
# rule definitions

rule all:
    input:
        'results/TAD_enrichment.csv',
        'images/tad_border_enrichment.pdf'

rule convert_tad_coordinates:
    input:
        'data/tads_hESC_hg19_with_ids.txt'
    output:
        'results/tads_hESC_hg38.tsv'
    run:
        execute_notebook('ConvertTADGenomicCoordinates.ipynb')

rule assemble_snp_database:
    input:
        'data/curated_variant_disease_associations.tsv.gz',
        'results/tads_hESC_hg38.tsv'
    output:
        'results/disgenet_enhanced_hg38.tsv',
        'results/disease_terms.csv'
    run:
        execute_notebook('LoadDisGeNET.ipynb')

rule compute_enrichments:
    input:
        'results/disgenet_enhanced_hg38.tsv',
        'results/tads_hESC_hg38.tsv'
    output:
        'results/TAD_enrichment.csv'
    run:
        execute_notebook('ComputeTADEnrichments.ipynb')

rule produce_figures:
    input:
        'results/disgenet_enhanced_hg38.tsv',
        'results/TAD_enrichment.csv'
    output:
        'images/tad_border_enrichment.pdf'
    run:
        execute_notebook('PublicationReproductions.ipynb')
