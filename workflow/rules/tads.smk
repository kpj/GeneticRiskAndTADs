rule install_extra_dependencies:
    """Ensure non-conda dependencies (e.g. TopDom R package) are available."""
    output:
        touch(RESULTS_DIR + "/dependencies/extra_dependencies_installed.marker"),
    conda:
        "../envs/python_stack.yaml"
    script:
        "../scripts/install_extra_dependencies.py"


rule compute_tads:
    """Identify TAD boundaries using TopDom for each chromosome contact matrix."""
    input:
        fname=RESULTS_DIR + "/hic_files/counts/{source}/{chromosome}/matrix.csv",
        fname_info=RESULTS_DIR + "/hic_files/info.csv",
        fname_dep_marker=(
            RESULTS_DIR + "/dependencies/extra_dependencies_installed.marker"
        ),
    output:
        fname=(
            RESULTS_DIR
            + "/tads/data/{source}/{tad_parameter}/tads.chr{chromosome}.csv"
        ),
        topdom_input=temp(
            RESULTS_DIR
            + "/tads/data/{source}/{tad_parameter}/topdom/topdom_input.chr{chromosome}.tsv"
        ),
        topdom_output=(
            RESULTS_DIR
            + "/tads/data/{source}/{tad_parameter}/topdom/topdom.chr{chromosome}.bed"
        ),
    conda:
        "../envs/python_stack.yaml"
    resources:
        mem_mb=30_000,
    params:
        prefix=(
            RESULTS_DIR
            + "/tads/data/{source}/{tad_parameter}/topdom/topdom.chr{chromosome}"
        ),
    script:
        "../scripts/compute_tads.py"


rule aggregate_tads:
    """Combine per-chromosome TAD predictions into genome-wide domain files."""
    input:
        fname_list=expand(
            RESULTS_DIR
            + "/tads/data/{source}/{tad_parameter}/tads.chr{chromosome}.csv",
            chromosome=config["chromosome_list"],
            allow_missing=True,
        ),
    output:
        fname=RESULTS_DIR + "/tads/data/tads.{source}.{tad_parameter}.csv",
    conda:
        "../envs/python_stack.yaml"
    script:
        "../scripts/aggregate_tads.py"


rule compare_tad_lists:
    """Evaluate concordance and overlaps between TAD calls across window sizes."""
    input:
        tad_fname_list=expand(
            RESULTS_DIR + "/tads/data/tads.{source}.{tad_parameter}.csv",
            source=hic_sources,
            tad_parameter=actual_window_size_list,
        ),
    output:
        outdir=directory(RESULTS_DIR + "/tads/plots/"),
    log:
        notebook=RESULTS_DIR + "/notebooks/TADListComparison.ipynb",
    conda:
        "../envs/python_stack.yaml"
    resources:
        mem_mb=4_000,
    notebook:
        "../notebooks/TADListComparison.ipynb"
