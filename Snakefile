from bioinf_common.tools import execute_notebook


###
# setup
configfile: 'config.yaml'

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
        config['input_files']['tad_coordinates_hg19']
    output:
        f'{results}/tads_hESC_hg38.tsv'
    run:
        execute_notebook('ConvertTADGenomicCoordinates.ipynb')

rule assemble_snp_database:
    input:
        config['input_files']['raw_disgenet'],
        config['input_files']['raw_gwascatalog'],
        f'{results}/tads_hESC_hg38.tsv'
    output:
        f'{results}/disgenet_enhanced.tsv',
        f'{results}/disease_cancer_classification.csv'
    run:
        execute_notebook('LoadDisGeNET.ipynb')

rule compute_enrichments:
    input:
        f'{results}/disgenet_enhanced.tsv',
        f'{results}/tads_hESC_hg38.tsv'
    output:
        f'{results}/TAD_enrichment.csv'
    run:
        execute_notebook('ComputeTADEnrichments.ipynb')

rule analyze_results:
    input:
        f'{results}/disgenet_enhanced.tsv',
        f'{results}/TAD_enrichment.csv',
        f'{results}/disease_cancer_classification.csv',
        f'{results}/tads_hESC_hg38.tsv'
    output:
        f'{images}/tad_border_enrichment.pdf'
    run:
        execute_notebook('PublicationReproductions.ipynb')
        execute_notebook('FurtherExperiments.ipynb')
