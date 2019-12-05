from pathlib import Path


# setup
configfile: 'config.yaml'
workdir: config['workdir']

localrules: all, download_hic_files


# rule definitions
rule all:
    input:
        'results/final.csv',
        'databases/statistics/',
        'tads/plots/',
        'post_analysis/',
        expand(
            'plots/subset/{source}/{tad_parameter}/',
            source=config['hic_sources'],
            tad_parameter=config['window_size_list'])


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
    script:
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
    script:
        'notebooks/AssembleInputDatabases.ipynb'


rule compute_database_statistics:
    input:
        fname = 'databases/initial.csv'
    output:
        outdir = directory('databases/statistics/'),
        notebook_output = 'notebooks/DatabaseStatistics.ipynb'
    script:
        'notebooks/DatabaseStatistics.ipynb'


rule include_tad_relations:
    input:
        tads_fname = 'tads/tads.{source}.{tad_parameter}.csv',
        db_fname = 'databases/initial.csv'
    output:
        db_fname = 'databases/per_source/snpdb.{source}.{tad_parameter}.csv',
        tad_length_plot = 'tads/tad_length_histogram.{source}.{tad_parameter}.pdf',
        notebook_output = 'notebooks/IncludeTADRelations.{source}.{tad_parameter}.ipynb'
    script:
        'notebooks/IncludeTADRelations.ipynb'


rule compute_enrichments:
    input:
        db_fname = 'databases/per_source/snpdb.{source}.{tad_parameter}.csv',
        tads_fname = 'tads/tads.{source}.{tad_parameter}.csv'
    output:
        fname = 'enrichments/results.{source}.{tad_parameter}.csv',
        notebook_output = 'notebooks/ComputeTADEnrichments.{source}.{tad_parameter}.ipynb'
    script:
        'notebooks/ComputeTADEnrichments.ipynb'


rule create_figures:
    input:
        db_fname = 'databases/per_source/snpdb.{source}.{tad_parameter}.csv',
        enr_fname = 'enrichments/results.{source}.{tad_parameter}.csv',
    output:
        outdir = directory('plots/subset/{source}/{tad_parameter}/'),
        notebook_output = 'notebooks/CreateFigures.{source}.{tad_parameter}.ipynb'
    script:
        'notebooks/CreateFigures.ipynb'


rule aggregate_results:
    input:
        database_files = expand(
            'databases/per_source/snpdb.{source}.{tad_parameter}.csv',
            source=config['hic_sources'],
            tad_parameter=config['window_size_list']),
        enrichment_files = expand(
            'enrichments/results.{source}.{tad_parameter}.csv',
            source=config['hic_sources'],
            tad_parameter=config['window_size_list'])
    output:
        fname = 'results/final.csv',
        notebook_output = 'notebooks/AggregateResults.ipynb'
    script:
        'notebooks/AggregateResults.ipynb'


rule multi_run_post_analysis:
    input:
        db_fname = 'results/final.csv'
    output:
        outdir = directory('post_analysis/'),
        notebook_output = 'notebooks/AggregateResults.ipynb'
    script:
        'notebooks/MultiRunPostAnalysis.ipynb'
