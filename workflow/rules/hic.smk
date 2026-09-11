rule extract_cool_from_mcool:
    """Extract a single resolution .cool matrix from a multi-resolution .mcool file."""
    input:
        mcool=lambda wc: sample_url_from_source_wildcard(wc.source),
    output:
        cool=RESULTS_DIR + "/hic_files/cool/{source}.cool",
    conda:
        "../envs/python_stack.yaml"
    resources:
        mem_mb=10_000,
        hdf5_lock=1,
    script:
        "../scripts/extract_cool_from_mcool.py"


rule extract_count_matrices:
    """Extract chromosome contact count matrix from a .cool file."""
    input:
        fname=lambda wc: cool_input_for_source_wildcard(wc.source),
    output:
        fname_matrix=RESULTS_DIR + "/hic_files/counts/{source}/{chromosome}/matrix.csv",
    conda:
        "../envs/python_stack.yaml"
    resources:
        mem_mb=10_000,
        hdf5_lock=1,
    script:
        "../scripts/extract_count_matrices.py"


rule gather_input_information:
    """Extract metadata (bin size, chromosome lengths, assembly) from Hi-C datasets."""
    input:
        fname_list=[url_wrapper(config["samples"][hic]) for hic in hic_sources_raw],
    output:
        fname=RESULTS_DIR + "/hic_files/info.csv",
    conda:
        "../envs/python_stack.yaml"
    resources:
        hdf5_lock=1,
    params:
        samplename_list=list(hic_sources_raw),
    script:
        "../scripts/gather_input_information.py"


rule visualize_count_matrix:
    """Generate interactive contact matrix heatmap per chromosome using papermill notebook."""
    input:
        fname_matrix=RESULTS_DIR + "/hic_files/counts/{source}/{chromosome}/matrix.csv",
    output:
        fname_heatmap=RESULTS_DIR
        + "/hic_files/plots/{source}/heatmap.chr{chromosome}.pdf",
    log:
        notebook=(
            RESULTS_DIR
            + "/notebooks/VisualizeContactMatrix.{source}.{chromosome}.ipynb"
        ),
    conda:
        "../envs/python_stack.yaml"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 20_000,
    notebook:
        "../notebooks/VisualizeContactMatrix.ipynb"
