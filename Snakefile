from snakemake.utils import validate


# setup
configfile: 'config.yaml'
validate(config, 'config.schema.yaml')

workdir: config['workdir']
include: 'rules/common.smk'

hic_sources = config['samples'].keys()

wildcard_constraints:
    tad_parameter = r'\d+'  # contains window size which must be an integer

localrules: all, gather_input_information, aggregate_tads, provide_input_files, compute_database_statistics, include_tad_relations, compute_enrichments, create_figures, create_report


# rule definitions
rule all:
    input:
        'results/final_enr.csv.gz',
        expand(
            'hic_files/plots/{source}/heatmap.chr{chromosome}.pdf',
            source=hic_sources, chromosome=config['chromosome_list']),
        expand('databases/statistics/'),
        'tads/plots/',
        'post_analysis/',
        'publication_figures/main/',
        'publication_figures/tad_plots/',
        expand(
            'plots/{source}/{tad_parameter}/{filter}/',
            source=hic_sources,
            tad_parameter=config['window_size_list'],
            filter=config['snp_filters'].keys()),
        # expand(
        #     'reports/report.{source}.{tad_parameter}.{filter}.pdf',
        #     source=hic_sources,
        #     tad_parameter=config['window_size_list'],
        #     filter=config['snp_filters'].keys())


rule extract_count_matrices:
    input:
        unpack(
            lambda wildcards:
                {'fname': url_wrapper(config['samples'][wildcards.source])}
        )
    output:
        fname_matrix = 'hic_files/counts/{source}/{chromosome}/matrix.csv'
    conda:
        'envs/python_stack.yaml'
    script:
        'scripts/extract_count_matrices.py'


rule gather_input_information:
    input:
        fname_list = [url_wrapper(config['samples'][hic])
                      for hic in hic_sources]
    output:
        fname = 'hic_files/info.csv'
    params:
        samplename_list = list(hic_sources)
    conda:
        'envs/python_stack.yaml'
    script:
        'scripts/gather_input_information.py'


rule visualize_count_matrix:
    input:
        fname_matrix = 'hic_files/counts/{source}/{chromosome}/matrix.csv'
    output:
        fname_heatmap = 'hic_files/plots/{source}/heatmap.chr{chromosome}.pdf'
    log:
        notebook = 'notebooks/VisualizeContactMatrix.{source}.{chromosome}.ipynb'
    conda:
        'envs/python_stack.yaml'
    notebook:
        'notebooks/VisualizeContactMatrix.ipynb'


rule compute_tads:
    input:
        fname = 'hic_files/counts/{source}/{chromosome}/matrix.csv',
        fname_info = 'hic_files/info.csv'
    output:
        fname = 'tads/data/{source}/{tad_parameter}/tads.chr{chromosome}.csv',
        topdom_input = temp('tads/data/{source}/{tad_parameter}/topdom/topdom_input.chr{chromosome}.tsv'),
        topdom_output = 'tads/data/{source}/{tad_parameter}/topdom/topdom.chr{chromosome}.bed'
    params:
        prefix = 'tads/data/{source}/{tad_parameter}/topdom/topdom.chr{chromosome}'
    conda:
        'envs/python_stack.yaml'
    script:
        'scripts/compute_tads.py'


rule aggregate_tads:
    input:
        fname_list = expand(
            'tads/data/{{source}}/{{tad_parameter}}/tads.chr{chromosome}.csv',
            chromosome=config['chromosome_list'])
    output:
        fname = 'tads/data/tads.{source}.{tad_parameter}.csv'
    conda:
        'envs/python_stack.yaml'
    script:
        'scripts/aggregate_tads.py'


rule compare_tad_lists:
    input:
        tad_fname_list = expand(
            'tads/data/tads.{source}.{tad_parameter}.csv',
            source=hic_sources,
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


rule include_tad_relations:
    input:
        tads_fname = 'tads/data/tads.{source}.{tad_parameter}.csv',
        db_fname = 'databases/initial.csv',
        info_fname = 'hic_files/info.csv'
    output:
        db_fname = 'databases/per_source/snpdb.{source}.{tad_parameter}.csv',
        tad_length_plot = 'tads/length_plots/tad_length_histogram.{source}.{tad_parameter}.pdf'
    log:
        notebook = 'notebooks/IncludeTADRelations.{source}.{tad_parameter}.ipynb'
    conda:
        'envs/python_stack.yaml'
    notebook:
        'notebooks/IncludeTADRelations.ipynb'


rule compute_enrichments:
    input:
        db_fname = 'databases/per_source/snpdb.{source}.{tad_parameter}.csv',
        tads_fname = 'tads/data/tads.{source}.{tad_parameter}.csv',
        info_fname = 'hic_files/info.csv'
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
        db_fname = 'databases/per_source/snpdb.{source}.{tad_parameter}.csv',
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
            'databases/per_source/snpdb.{source}.{tad_parameter}.csv',
            source=hic_sources,
            tad_parameter=config['window_size_list']),
        enrichment_files = expand(
            'enrichments/results.{source}.{tad_parameter}.{filter}.csv',
            source=hic_sources,
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


rule compute_database_statistics:
    input:
        fname = 'results/final_data.csv.gz'
    output:
        outdir = directory('results/database_statistics/')
    log:
        notebook = 'notebooks/DatabaseStatistics.ipynb'
    conda:
        'envs/python_stack.yaml'
    notebook:
        'notebooks/DatabaseStatistics.ipynb'


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


rule publication_figures:
    input:
        fname_data = 'results/final_data.csv.gz',
        fname_enr = 'results/final_enr.csv.gz',

        sketch_hicfile = url_wrapper(
            config['samples'][config['sketch']['data_source']]),
        sketch_tadfile = 'tads/data/tads.{source}.{tad_parameter}.csv'.format(
            source=config['sketch']['data_source'],
            tad_parameter=config['sketch']['window_size'])
    output:
        outdir = directory('publication_figures/main/')
    log:
        notebook = 'notebooks/PublicationFigures.ipynb'
    conda:
        'envs/python_stack.yaml'
    notebook:
        'notebooks/PublicationFigures.ipynb'


rule supplementary_tadplots:
    input:
        fname_data = 'results/final_data.csv.gz',

        sketch_hicfile = url_wrapper(
            config['samples'][config['sketch']['data_source']]),
        sketch_tadfile = 'tads/data/tads.{source}.{tad_parameter}.csv'.format(
            source=config['sketch']['data_source'],
            tad_parameter=config['sketch']['window_size']),

        tad_fname_list = expand(
            'tads/data/tads.{source}.{tad_parameter}.csv',
            source=hic_sources,
            tad_parameter=config['window_size_list'])
    output:
        outdir = directory('publication_figures/tad_plots')
    log:
        notebook = 'notebooks/Supplementaries_TADPlots.ipynb'
    conda:
        'envs/python_stack.yaml'
    notebook:
        'notebooks/Supplementaries_TADPlots.ipynb'
