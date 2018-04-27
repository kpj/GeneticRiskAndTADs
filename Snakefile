from bioinf_common.tools import execute_notebook


###
# setup
configfile: 'config.yaml'


###
# rule definitions

rule all:
    input:
        'results/TAD_enrichment.csv',
        'images/tad_border_enrichment.pdf'

rule convert_tad_coordinates:
    input:
        config['input_files']['tad_coordinates_hg19']
    output:
        'results/tads_hESC_hg38.tsv'
    run:
        execute_notebook('ConvertTADGenomicCoordinates.ipynb')

rule assemble_snp_database:
    input:
        config['input_files']['raw_disgenet'],
        config['input_files']['raw_gwascatalog'],
        'results/tads_hESC_hg38.tsv'
    output:
        'results/disgenet_enhanced_hg38.tsv'
    run:
        execute_notebook('LoadDisGeNET.ipynb')

rule compute_enrichments:
    input:
        'results/disgenet_enhanced_hg38.tsv',
        'results/tads_hESC_hg38.tsv'
    output:
        'results/TAD_enrichment.csv',
        'results/disease_cancer_classification.csv'
    run:
        execute_notebook('ComputeTADEnrichments.ipynb')

rule produce_figures:
    input:
        'results/disgenet_enhanced_hg38.tsv',
        'results/TAD_enrichment.csv',
        'results/disease_cancer_classification.csv'
    output:
        'images/tad_border_enrichment.pdf'
    run:
        execute_notebook('PublicationReproductions.ipynb')
