rule include_tad_relations:
    """Map SNPs to TAD interiors and boundary borders and plot length distributions."""
    input:
        tads_fname=RESULTS_DIR + "/tads/data/tads.{source}.{tad_parameter}.csv",
        db_fname=RESULTS_DIR + "/databases/initial.csv",
        info_fname=RESULTS_DIR + "/hic_files/info.csv",
    output:
        db_fname=(
            RESULTS_DIR + "/databases/per_source/snpdb.{source}.{tad_parameter}.csv"
        ),
        tad_length_plot=(
            RESULTS_DIR
            + "/tads/length_plots/tad_length_histogram.{source}.{tad_parameter}.pdf"
        ),
    log:
        notebook=(
            RESULTS_DIR
            + "/notebooks/IncludeTADRelations.{source}.{tad_parameter}.ipynb"
        ),
    conda:
        "../envs/python_stack.yaml"
    resources:
        notebook_slots=1,
    notebook:
        "../notebooks/IncludeTADRelations.ipynb"


rule snp_majority_vote:
    """Derive broad consensus TAD assignments across all evaluated window sizes."""
    input:
        fname_list=expand(
            RESULTS_DIR + "/databases/per_source/snpdb.{source}.{tad_parameter}.csv",
            tad_parameter=actual_window_size_list,
            allow_missing=True,
        ),
    output:
        fname=RESULTS_DIR + "/databases/per_source/snpdb.{source}.0.csv",
        fname_tads=RESULTS_DIR + "/tads/data/tads.{source}.0.csv",  # empty dummy
    log:
        notebook=RESULTS_DIR + "/notebooks/SNPMajorityVote.{source}.0.ipynb",
    conda:
        "../envs/python_stack.yaml"
    resources:
        mem_mb=3_000,
    params:
        tad_parameter_range=None,
        border_fraction_threshold=config["parameters"][
            "snp_majority_vote_border_fraction_threshold"
        ],
    notebook:
        "../notebooks/SNPMajorityVote.ipynb"


rule snp_majority_vote_narrow:
    """Derive narrow consensus TAD assignments across focused window size range."""
    input:
        fname_list=expand(
            RESULTS_DIR + "/databases/per_source/snpdb.{source}.{tad_parameter}.csv",
            tad_parameter=actual_window_size_list,
            allow_missing=True,
        ),
    output:
        fname=RESULTS_DIR + "/databases/per_source/snpdb.{source}.1.csv",
        fname_tads=RESULTS_DIR + "/tads/data/tads.{source}.1.csv",  # empty dummy
    log:
        notebook=RESULTS_DIR + "/notebooks/SNPMajorityVote.{source}.1.ipynb",
    conda:
        "../envs/python_stack.yaml"
    resources:
        mem_mb=3_000,
    params:
        tad_parameter_range=config["parameters"]["snp_majority_vote_narrow"][
            "tad_parameter_range"
        ],
        border_fraction_threshold=config["parameters"]["snp_majority_vote_narrow"][
            "border_fraction_threshold"
        ],
    notebook:
        "../notebooks/SNPMajorityVote.ipynb"


rule compute_enrichments:
    """Calculate statistical enrichment of disease SNPs inside TAD borders."""
    input:
        db_fname=(
            RESULTS_DIR + "/databases/per_source/snpdb.{source}.{tad_parameter}.csv"
        ),
        tads_fname=RESULTS_DIR + "/tads/data/tads.{source}.{tad_parameter}.csv",
        info_fname=RESULTS_DIR + "/hic_files/info.csv",
    output:
        fname=(
            RESULTS_DIR + "/enrichments/results.{source}.{tad_parameter}.{filter}.csv"
        ),
    log:
        notebook=(
            RESULTS_DIR
            + "/notebooks/ComputeTADEnrichments.{source}.{tad_parameter}.{filter}.ipynb"
        ),
    conda:
        "../envs/python_stack.yaml"
    resources:
        notebook_slots=1,
    notebook:
        "../notebooks/ComputeTADEnrichments.ipynb"


rule aggregate_results:
    """Aggregate all database annotations and enrichment statistics into compressed tables."""
    input:
        database_files=expand(
            RESULTS_DIR + "/databases/per_source/snpdb.{source}.{tad_parameter}.csv",
            source=hic_sources,
            tad_parameter=config["window_size_list"],
        ),
        enrichment_files=expand(
            RESULTS_DIR + "/enrichments/results.{source}.{tad_parameter}.{filter}.csv",
            source=hic_sources,
            tad_parameter=config["window_size_list"],
            filter=config["snp_filters"].keys(),
        ),
    output:
        fname_data=RESULTS_DIR + "/results/final_data.csv.gz",
        fname_enr=report(
            RESULTS_DIR + "/results/final_enr.csv.gz",
            caption="../report/final_enr.rst",
            category="Enrichment Results",
        ),
    log:
        notebook=RESULTS_DIR + "/notebooks/AggregateResults.ipynb",
    conda:
        "../envs/python_stack.yaml"
    resources:
        mem_mb=lambda wildcards, attempt: (attempt + 2) * 10_000,
    notebook:
        "../notebooks/AggregateResults.ipynb"
