from pathlib import Path


# setup
configfile: 'config.yaml'
workdir: config['workdir']


# helper functions
def aggregate_tad_files(wildcards):
    checkpoint_output = checkpoints.extract_count_matrices.get(**wildcards).output.out_dir

    chroms = glob_wildcards(Path(checkpoint_output) / '{chromosome}' / 'matrix.csv').chromosome

    return expand(
        'tads/{source}/tads.chr{chromosome}.csv',
        source=wildcards.source,
        chromosome=chroms
    )


# rule definitions
rule all:
    input:
        'results/final.csv',
        'databases/statistics/',
        'tads/plots/',
        'post_analysis/',
        expand('plots/subset/{source}/', source=config['hic_sources'])


checkpoint download_hic_files:
    output:
        fname = 'hic_files/raw/data.{source}.hic'
    script:
        'scripts/download_hic_files.py'


checkpoint extract_count_matrices:
    input:
        fname = 'hic_files/raw/data.{source}.hic'
    output:
        out_dir = directory('hic_files/counts/{source}/')
    script:
        'scripts/extract_count_matrices.py'


rule compute_tads:
    input:
        fname = 'hic_files/counts/{source}/{chromosome}/matrix.csv'
    output:
        fname = 'tads/{source}/tads.chr{chromosome}.csv',
        topdom_input = 'tads/{source}/topdom_input.chr{chromosome}.tsv',
        topdom_output = 'tads/{source}/topdom.chr{chromosome}.bed'
    params:
        prefix = 'tads/{source}/topdom.chr{chromosome}'
    script:
        'scripts/compute_tads.py'


rule aggregate_tads:
    input:
        fname_list = aggregate_tad_files
    output:
        fname = 'tads/tads.{source}.csv'
    run:
        import pandas as pd

        pd.concat([
            pd.read_csv(x) for x in input.fname_list
        ]).to_csv(output.fname, index=False)


rule compare_tad_lists:
    input:
        tad_fname_list = expand(
            'tads/tads.{source}.csv', source=config['hic_sources'])
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
        tads_fname = 'tads/tads.{source}.csv',
        db_fname = 'databases/initial.csv'
    output:
        db_fname = 'databases/per_source/snpdb.{source}.csv',
        tad_length_plot = 'tads/tad_length_histogram.{source}.pdf',
        notebook_output = 'notebooks/IncludeTADRelations.{source}.ipynb'
    script:
        'notebooks/IncludeTADRelations.ipynb'


rule compute_enrichments:
    input:
        db_fname = 'databases/per_source/snpdb.{source}.csv',
        tads_fname = 'tads/tads.{source}.csv'
    output:
        fname = 'enrichments/results.{source}.csv',
        notebook_output = 'notebooks/ComputeTADEnrichments.{source}.ipynb'
    script:
        'notebooks/ComputeTADEnrichments.ipynb'


rule create_figures:
    input:
        db_fname = 'databases/per_source/snpdb.{source}.csv',
        enr_fname = 'enrichments/results.{source}.csv',
    output:
        outdir = directory('plots/subset/{source}/'),
        notebook_output = 'notebooks/CreateFigures.{source}.ipynb'
    script:
        'notebooks/CreateFigures.ipynb'


rule aggregate_results:
    input:
        database_files = expand(
            'databases/per_source/snpdb.{source}.csv',
            source=config['hic_sources']),
        enrichment_files = expand(
            'enrichments/results.{source}.csv',
            source=config['hic_sources'])
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
