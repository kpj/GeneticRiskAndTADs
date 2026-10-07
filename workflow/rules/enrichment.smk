rule include_tad_relations:
    """Map SNPs to TAD interiors and boundary borders and plot length distributions."""
    input:
        tads_fname=RESULTS_DIR + "/tads/data/tads.{source}.{caller_config}.csv",
        db_fname=RESULTS_DIR + "/databases/initial.csv",
        info_fname=RESULTS_DIR + "/hic_files/info.csv",
    output:
        db_fname=(
            RESULTS_DIR + "/databases/per_source/snpdb.{source}.{caller_config}.csv"
        ),
        tad_length_plot=(
            RESULTS_DIR
            + "/tads/length_plots/tad_length_histogram.{source}.{caller_config}.pdf"
        ),
    log:
        notebook=(
            RESULTS_DIR
            + "/notebooks/IncludeTADRelations.{source}.{caller_config}.ipynb"
        ),
    conda:
        "../envs/python_stack.yaml"
    resources:
        notebook_slots=1,
    notebook:
        "../notebooks/IncludeTADRelations.ipynb"


def get_majority_vote_inputs(source: str, consensus_id: str) -> list[str]:
    cfg = resolved_majority_votes[consensus_id]
    participating = cfg["resolved_caller_configs"]
    return [
        RESULTS_DIR + f"/databases/per_source/snpdb.{source}.{c}.csv"
        for c in participating
    ]


def get_majority_vote_threshold(consensus_id: str) -> float:
    return float(
        resolved_majority_votes[consensus_id].get("border_fraction_threshold", 0.3)
    )


rule snp_majority_vote:
    """Derive consensus TAD assignments across specified caller configs."""
    input:
        fname_list=lambda wildcards: get_majority_vote_inputs(
            wildcards.source, wildcards.consensus_id
        ),
    output:
        fname=RESULTS_DIR + "/databases/per_source/snpdb.{source}.{consensus_id}.csv",
        fname_tads=RESULTS_DIR + "/tads/data/tads.{source}.{consensus_id}.csv",  # empty dummy
    log:
        notebook=RESULTS_DIR + "/notebooks/SNPMajorityVote.{source}.{consensus_id}.ipynb",
    conda:
        "../envs/python_stack.yaml"
    resources:
        mem_mb=3_000,
    params:
        border_fraction_threshold=lambda wildcards: get_majority_vote_threshold(
            wildcards.consensus_id
        ),
    notebook:
        "../notebooks/SNPMajorityVote.ipynb"


rule compute_enrichments:
    """Calculate statistical enrichment of disease SNPs inside TAD borders."""
    input:
        db_fname=(
            RESULTS_DIR + "/databases/per_source/snpdb.{source}.{caller_config}.csv"
        ),
        tads_fname=RESULTS_DIR + "/tads/data/tads.{source}.{caller_config}.csv",
        info_fname=RESULTS_DIR + "/hic_files/info.csv",
    output:
        fname=(
            RESULTS_DIR + "/enrichments/results.{source}.{caller_config}.{filter}.csv"
        ),
    log:
        notebook=(
            RESULTS_DIR
            + "/notebooks/ComputeTADEnrichments.{source}.{caller_config}.{filter}.ipynb"
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
            RESULTS_DIR + "/databases/per_source/snpdb.{source}.{caller_config}.csv",
            source=hic_sources,
            caller_config=all_caller_configs,
        ),
        enrichment_files=expand(
            RESULTS_DIR + "/enrichments/results.{source}.{caller_config}.{filter}.csv",
            source=hic_sources,
            caller_config=all_caller_configs,
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
