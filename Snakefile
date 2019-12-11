from pathlib import Path


# setup
configfile: 'config.yaml'
workdir: config['workdir']

localrules: all, download_hic_files, aggregate_tads, filter_database, compute_database_statistics, include_tad_relations, compute_enrichments, create_figures


# rule definitions
rule all:
    input:
        'results/final.csv.gz',
        expand(
            'databases/statistics/{filter}/',
            filter=config['snp_filters'].keys()),
        'tads/plots/',
        'post_analysis/',
        expand(
            'plots/{source}/{tad_parameter}/{filter}/',
            source=config['hic_sources'],
            tad_parameter=config['window_size_list'],
            filter=config['snp_filters'].keys())


rule download_hic_files:
    output:
        fname = 'hic_files/raw/data.{source}.hic'
    script:
        'scripts/download_hic_files.py'


rule extract_count_matrices:
    input:
        fname = 'hic_files/raw/data.{source}.hic'
    output:
        fname_matrix = 'hic_files/counts/{source}/{chromosome}/matrix.csv',
        fname_juicer = 'hic_files/counts/{source}/{chromosome}/juicer.tsv'
    script:
        'scripts/extract_count_matrices.py'


rule compute_tads:
    input:
        fname = 'hic_files/counts/{source}/{chromosome}/matrix.csv'
    output:
        fname = 'tads/{source}/{tad_parameter}/tads.chr{chromosome}.csv',
        topdom_input = 'tads/{source}/{tad_parameter}/topdom/topdom_input.chr{chromosome}.tsv',
        topdom_output = 'tads/{source}/{tad_parameter}/topdom/topdom.chr{chromosome}.bed'
    params:
        prefix = 'tads/{source}/{tad_parameter}/topdom/topdom.chr{chromosome}'
    script:
        'scripts/compute_tads.py'


rule aggregate_tads:
    input:
        fname_list = expand(
            'tads/{{source}}/{{tad_parameter}}/tads.chr{chromosome}.csv',
            chromosome=config['chromosome_list'])
    output:
        fname = 'tads/tads.{source}.{tad_parameter}.csv'
    run:
        import pandas as pd

        pd.concat([
            pd.read_csv(x) for x in input.fname_list
        ]).to_csv(output.fname, index=False)


rule compare_tad_lists:
    input:
        tad_fname_list = expand(
            'tads/tads.{source}.{tad_parameter}.csv',
            source=config['hic_sources'],
            tad_parameter=config['window_size_list'])
    output:
        tad_similarity_cache = 'tads/statistics/tad_similarities.csv',
        outdir = directory('tads/plots/'),
        notebook_output = 'notebooks/TADListComparison.ipynb'
    notebook:
        'notebooks/TADListComparison.ipynb'


rule assemble_input_databases:
    input:
        disgenet_fname = srcdir(config['input_files']['raw_disgenet']),
        gwascatalog_fname = srcdir(config['input_files']['raw_gwascatalog']),
        efo_fname = srcdir(config['input_files']['exp_factor_ontology']),
        so_fname = srcdir(config['input_files']['sequence_ontology'])
    output:
        db_fname = 'databases/initial.csv',
        raw_veps = 'databases/vep.csv',
        notebook_output = 'notebooks/AssembleInputDatabases.ipynb'
    notebook:
        'notebooks/AssembleInputDatabases.ipynb'


rule filter_database:
    input:
        db_fname = 'databases/initial.csv'
    output:
        db_fname = 'databases/initial_filtered.{filter}.csv',
        notebook_output = 'notebooks/FilterDatabase.{filter}.ipynb'
    notebook:
        'notebooks/FilterDatabase.ipynb'


rule compute_database_statistics:
    input:
        fname = 'databases/initial_filtered.{filter}.csv'
    output:
        outdir = directory('databases/statistics/{filter}/'),
        notebook_output = 'notebooks/DatabaseStatistics.{filter}.ipynb'
    notebook:
        'notebooks/DatabaseStatistics.ipynb'


rule include_tad_relations:
    input:
        tads_fname = 'tads/tads.{source}.{tad_parameter}.csv',
        db_fname = 'databases/initial_filtered.{filter}.csv'
    output:
        db_fname = 'databases/per_source/snpdb.{source}.{tad_parameter}.{filter}.csv',
        tad_length_plot = 'tads/length_plots/tad_length_histogram.{source}.{tad_parameter}.{filter}.pdf',
        notebook_output = 'notebooks/IncludeTADRelations.{source}.{tad_parameter}.{filter}.ipynb'
    notebook:
        'notebooks/IncludeTADRelations.ipynb'


rule compute_enrichments:
    input:
        db_fname = 'databases/per_source/snpdb.{source}.{tad_parameter}.{filter}.csv',
        tads_fname = 'tads/tads.{source}.{tad_parameter}.csv'
    output:
        fname = 'enrichments/results.{source}.{tad_parameter}.{filter}.csv',
        notebook_output = 'notebooks/ComputeTADEnrichments.{source}.{tad_parameter}.{filter}.ipynb'
    notebook:
        'notebooks/ComputeTADEnrichments.ipynb'


rule create_figures:
    input:
        db_fname = 'databases/per_source/snpdb.{source}.{tad_parameter}.{filter}.csv',
        enr_fname = 'enrichments/results.{source}.{tad_parameter}.{filter}.csv',
    output:
        outdir = directory('plots/{source}/{tad_parameter}/{filter}/'),
        notebook_output = 'notebooks/CreateFigures.{source}.{tad_parameter}.{filter}.ipynb'
    notebook:
        'notebooks/CreateFigures.ipynb'


rule aggregate_results:
    input:
        database_files = expand(
            'databases/per_source/snpdb.{source}.{tad_parameter}.{filter}.csv',
            source=config['hic_sources'],
            tad_parameter=config['window_size_list'],
            filter=config['snp_filters'].keys()),
        enrichment_files = expand(
            'enrichments/results.{source}.{tad_parameter}.{filter}.csv',
            source=config['hic_sources'],
            tad_parameter=config['window_size_list'],
            filter=config['snp_filters'].keys())
    output:
        fname = 'results/final.csv.gz',
        notebook_output = 'notebooks/AggregateResults.ipynb'
    notebook:
        'notebooks/AggregateResults.ipynb'


rule multi_run_post_analysis:
    input:
        db_fname = 'results/final.csv.gz'
    output:
        outdir = directory('post_analysis/'),
        notebook_output = 'notebooks/MultiRunPostAnalysis.ipynb'
    notebook:
        'notebooks/MultiRunPostAnalysis.ipynb'
