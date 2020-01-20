from urllib.parse import urlparse

from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider

HTTP = HTTPRemoteProvider()
FTP = FTPRemoteProvider()


# setup
configfile: 'config.yaml'
workdir: config['workdir']

localrules: all, download_hic_files, aggregate_tads, provide_input_files, filter_database, compute_database_statistics, include_tad_relations, compute_enrichments, create_figures, create_report


def url_wrapper(url):
    if os.path.isfile(srcdir(url)):
        # is local
        return srcdir(url)
    else:
        # is remote
        o = urlparse(url)
        if o.scheme in ('http', 'https'):
            return HTTP.remote(url)
        elif o.scheme == 'ftp':
            return FTP.remote(url)
        else:
            raise RuntimeError(f'Invalid url: "{url}"')


# rule definitions
rule all:
    input:
        'results/final_enr.csv.gz',
        expand(
            'hic_files/plots/{source}/heatmap.{chromosome}.pdf',
            source=config['hic_sources'], chromosome=config['chromosome_list']),
        expand(
            'databases/statistics/{filter}/',
            filter=config['snp_filters'].keys()),
        'tads/plots/',
        'post_analysis/',
        expand(
            'plots/{source}/{tad_parameter}/{filter}/',
            source=config['hic_sources'],
            tad_parameter=config['window_size_list'],
            filter=config['snp_filters'].keys()),
        expand(
            'reports/report.{source}.{tad_parameter}.{filter}.pdf',
            source=config['hic_sources'],
            tad_parameter=config['window_size_list'],
            filter=config['snp_filters'].keys())


rule download_hic_files:
    output:
        fname = 'hic_files/raw/data.{source}.hic'
    conda:
        'envs/python_stack.yaml'
    script:
        'scripts/download_hic_files.py'


rule extract_count_matrices:
    input:
        fname = 'hic_files/raw/data.{source}.hic'
    output:
        fname_matrix = 'hic_files/counts/{source}/{chromosome}/matrix.csv',
        fname_juicer = 'hic_files/counts/{source}/{chromosome}/juicer.tsv'
    conda:
        'envs/python_stack.yaml'
    script:
        'scripts/extract_count_matrices.py'


rule visualize_count_matrix:
    input:
        fname_matrix = 'hic_files/counts/{source}/{chromosome}/matrix.csv'
    output:
        fname_heatmap = 'hic_files/plots/{source}/heatmap.{chromosome}.pdf'
    log:
        notebook = 'notebooks/VisualizeContactMatrix.{source}.{chromosome}.ipynb'
    conda:
        'envs/python_stack.yaml'
    notebook:
        'nnotebooks/VisualizeContactMatrix.{source}.{chromosome}.ipynb'


rule compute_tads:
    input:
        fname = 'hic_files/counts/{source}/{chromosome}/matrix.csv'
    output:
        fname = 'tads/{source}/{tad_parameter}/tads.chr{chromosome}.csv',
        topdom_input = 'tads/{source}/{tad_parameter}/topdom/topdom_input.chr{chromosome}.tsv',
        topdom_output = 'tads/{source}/{tad_parameter}/topdom/topdom.chr{chromosome}.bed'
    params:
        prefix = 'tads/{source}/{tad_parameter}/topdom/topdom.chr{chromosome}'
    conda:
        'envs/python_stack.yaml'
    script:
        'scripts/compute_tads.py'


rule aggregate_tads:
    input:
        fname_list = expand(
            'tads/{{source}}/{{tad_parameter}}/tads.chr{chromosome}.csv',
            chromosome=config['chromosome_list'])
    output:
        fname = 'tads/tads.{source}.{tad_parameter}.csv'
    conda:
        'envs/python_stack.yaml'
    script:
        'scripts/aggregate_tads.py'


rule compare_tad_lists:
    input:
        tad_fname_list = expand(
            'tads/tads.{source}.{tad_parameter}.csv',
            source=config['hic_sources'],
            tad_parameter=config['window_size_list'])
    output:
        tad_similarity_cache = 'tads/statistics/tad_similarities.csv',
        outdir = directory('tads/plots/')
    log:
        notebook = 'notebooks/TADListComparison.ipynb'
    conda:
        'envs/python_stack.yaml'
    notebook:
        'notebooks/TADListComparison.ipynb'


rule provide_input_files:
    input:
        disgenet_fname = url_wrapper(config['input_files']['disgenet']),
        gwascatalog_fname = url_wrapper(config['input_files']['gwascatalog']),
        efo_fname = url_wrapper(config['input_files']['exp_factor_ontology']),
        so_fname = url_wrapper(config['input_files']['sequence_ontology'])
    output:
        disgenet_fname = 'input/disgenet.tsv.gz',
        gwascatalog_fname = 'input/gwas_catalog.tsv.gz',
        efo_fname = 'input/efo.owl',
        so_fname = 'input/so.owl'
    conda:
        'envs/python_stack.yaml'
    script:
        'scripts/provide_files.py'


rule assemble_input_databases:
    input:
        disgenet_fname = 'input/disgenet.tsv.gz',
        gwascatalog_fname = 'input/gwas_catalog.tsv.gz',
        efo_fname = 'input/efo.owl',
        so_fname = 'input/so.owl'
    output:
        db_fname = 'databases/initial.csv',
        raw_veps = 'databases/vep.csv'
    log:
        notebook = 'notebooks/AssembleInputDatabases.ipynb'
    conda:
        'envs/python_stack.yaml'
    notebook:
        'notebooks/AssembleInputDatabases.ipynb'


rule filter_database:
    input:
        db_fname = 'databases/initial.csv'
    output:
        db_fname = 'databases/initial_filtered.{filter}.csv'
    log:
        notebook = 'notebooks/FilterDatabase.{filter}.ipynb'
    conda:
        'envs/python_stack.yaml'
    notebook:
        'notebooks/FilterDatabase.ipynb'


rule compute_database_statistics:
    input:
        fname = 'databases/initial_filtered.{filter}.csv'
    output:
        outdir = directory('databases/statistics/{filter}/')
    log:
        notebook = 'notebooks/DatabaseStatistics.{filter}.ipynb'
    conda:
        'envs/python_stack.yaml'
    notebook:
        'notebooks/DatabaseStatistics.ipynb'


rule include_tad_relations:
    input:
        tads_fname = 'tads/tads.{source}.{tad_parameter}.csv',
        db_fname = 'databases/initial_filtered.{filter}.csv'
    output:
        db_fname = 'databases/per_source/snpdb.{source}.{tad_parameter}.{filter}.csv',
        tad_length_plot = 'tads/length_plots/tad_length_histogram.{source}.{tad_parameter}.{filter}.pdf'
    log:
        notebook = 'notebooks/IncludeTADRelations.{source}.{tad_parameter}.{filter}.ipynb'
    conda:
        'envs/python_stack.yaml'
    notebook:
        'notebooks/IncludeTADRelations.ipynb'


rule compute_enrichments:
    input:
        db_fname = 'databases/per_source/snpdb.{source}.{tad_parameter}.{filter}.csv',
        tads_fname = 'tads/tads.{source}.{tad_parameter}.csv'
    output:
        fname = 'enrichments/results.{source}.{tad_parameter}.{filter}.csv'
    log:
        notebook = 'notebooks/ComputeTADEnrichments.{source}.{tad_parameter}.{filter}.ipynb'
    conda:
        'envs/python_stack.yaml'
    notebook:
        'notebooks/ComputeTADEnrichments.ipynb'


rule create_figures:
    input:
        db_fname = 'databases/per_source/snpdb.{source}.{tad_parameter}.{filter}.csv',
        enr_fname = 'enrichments/results.{source}.{tad_parameter}.{filter}.csv',
    output:
        outdir = directory('plots/{source}/{tad_parameter}/{filter}/')
    log:
        notebook = 'notebooks/CreateFigures.{source}.{tad_parameter}.{filter}.ipynb'
    conda:
        'envs/python_stack.yaml'
    notebook:
        'notebooks/CreateFigures.ipynb'


rule create_report:
    input:
        fname_enr = 'enrichments/results.{source}.{tad_parameter}.{filter}.csv'
    output:
        'reports/report.{source}.{tad_parameter}.{filter}.pdf'
    conda:
        'envs/r_stack.yaml'
    script:
        'report/report.Rmd'


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
        fname_data = 'results/final_data.csv.gz',
        fname_enr = 'results/final_enr.csv.gz'
    log:
        notebook = 'notebooks/AggregateResults.ipynb'
    conda:
        'envs/python_stack.yaml'
    notebook:
        'notebooks/AggregateResults.ipynb'


rule multi_run_post_analysis:
    input:
        fname_data = 'results/final_data.csv.gz',
        fname_enr = 'results/final_enr.csv.gz',
    output:
        outdir = directory('post_analysis/')
    log:
        notebook = 'notebooks/MultiRunPostAnalysis.ipynb'
    conda:
        'envs/python_stack.yaml'
    notebook:
        'notebooks/MultiRunPostAnalysis.ipynb'
