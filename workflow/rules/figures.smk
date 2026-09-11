rule create_figures:
    """Generate per-dataset and per-filter diagnostic enrichment plots."""
    input:
        db_fname=(
            RESULTS_DIR + "/databases/per_source/snpdb.{source}.{tad_parameter}.csv"
        ),
        enr_fname=(
            RESULTS_DIR + "/enrichments/results.{source}.{tad_parameter}.{filter}.csv"
        ),
    output:
        outdir=directory(RESULTS_DIR + "/plots/{source}/{tad_parameter}/{filter}/"),
    log:
        notebook=(
            RESULTS_DIR
            + "/notebooks/CreateFigures.{source}.{tad_parameter}.{filter}.ipynb"
        ),
    conda:
        "../envs/python_stack.yaml"
    resources:
        notebook_slots=1,
    notebook:
        "../notebooks/CreateFigures.ipynb"


rule create_report:
    """Compile PDF summary report for specific source and TAD parameters."""
    input:
        fname_enr=(
            RESULTS_DIR + "/enrichments/results.{source}.{tad_parameter}.{filter}.csv"
        ),
    output:
        RESULTS_DIR + "/reports/report.{source}.{tad_parameter}.{filter}.pdf",
    conda:
        "../envs/r_stack.yaml"
    script:
        "../report/report.Rmd"


rule compute_database_statistics:
    """Generate summary distributions and disease coverage plots for assembled databases."""
    input:
        fname=RESULTS_DIR + "/results/final_data.csv.gz",
        tad_fname_list=expand(
            RESULTS_DIR + "/tads/data/tads.{source}.{tad_parameter}.csv",
            source=hic_sources,
            tad_parameter=actual_window_size_list,
        ),
    output:
        outdir=report(
            directory(RESULTS_DIR + "/results/database_statistics/"),
            patterns=["{name}.pdf"],
            caption="../report/database_statistics.rst",
            category="Database Statistics",
        ),
    log:
        notebook=RESULTS_DIR + "/notebooks/DatabaseStatistics.ipynb",
    conda:
        "../envs/python_stack.yaml"
    resources:
        mem_mb=lambda wildcards, attempt: (attempt + 2) * 10_000,
    notebook:
        "../notebooks/DatabaseStatistics.ipynb"


rule multi_run_post_analysis:
    """Perform meta-analysis and comparative evaluations across multiple pipeline runs."""
    input:
        fname_data=RESULTS_DIR + "/results/final_data.csv.gz",
        fname_enr=RESULTS_DIR + "/results/final_enr.csv.gz",
    output:
        outdir=directory(RESULTS_DIR + "/post_analysis/"),
    log:
        notebook=RESULTS_DIR + "/notebooks/MultiRunPostAnalysis.ipynb",
    conda:
        "../envs/python_stack.yaml"
    resources:
        mem_mb=lambda wildcards, attempt: (attempt + 2) * 10_000,
    notebook:
        "../notebooks/MultiRunPostAnalysis.ipynb"


rule publication_figures:
    """Render publication-grade manuscript figures combining Hi-C contact maps with TAD calls and enrichments."""
    input:
        fname_data=RESULTS_DIR + "/results/final_data.csv.gz",
        fname_enr=RESULTS_DIR + "/results/final_enr.csv.gz",
        sketch_hicfile=cool_input_for_source_wildcard(config["sketch"]["data_source"]),
        sketch_tadfile=(
            RESULTS_DIR + "/tads/data/tads.{source}.{tad_parameter}.csv"
        ).format(
            source=config["sketch"]["data_source"],
            tad_parameter=config["sketch"]["window_size"],
        ),
    output:
        outdir=report(
            directory(RESULTS_DIR + "/publication_figures/main/"),
            patterns=["{name}.pdf"],
            caption="../report/publication_figures.rst",
            category="Publication Figures",
        ),
    log:
        notebook=RESULTS_DIR + "/notebooks/PublicationFigures.ipynb",
    conda:
        "../envs/python_stack.yaml"
    resources:
        mem_mb=lambda wildcards, attempt: (attempt + 5) * 10_000,
        hdf5_lock=1,
    notebook:
        "../notebooks/PublicationFigures.ipynb"


rule supplementary_tadplots_multidataset:
    """Generate multi-dataset locus comparison figures across diverse Hi-C sources."""
    input:
        fname_data=RESULTS_DIR + "/results/final_data.csv.gz",
        sketch_hicfile=cool_input_for_source_wildcard(config["sketch"]["data_source"]),
        sketch_tadfile=(
            RESULTS_DIR + "/tads/data/tads.{source}.{tad_parameter}.csv"
        ).format(
            source=config["sketch"]["data_source"],
            tad_parameter=config["sketch"]["window_size"],
        ),
        tad_fname_list=expand(
            RESULTS_DIR + "/tads/data/tads.{source}.{tad_parameter}.csv",
            source=hic_sources,
            tad_parameter=config["parameters"].get("supplementary_tad_parameter", 10),
        ),
    output:
        outdir=report(
            directory(RESULTS_DIR + "/publication_figures/tad_plots_multidataset/"),
            patterns=["{name}.pdf"],
            caption="../report/supplementary_plots.rst",
            category="Supplementary TAD Plots",
            subcategory="Multi-Dataset",
        ),
    log:
        notebook=RESULTS_DIR + "/notebooks/Supplementaries_TADPlots.ipynb",
    conda:
        "../envs/python_stack.yaml"
    resources:
        mem_mb=lambda wildcards, attempt: (attempt + 2) * 10_000,
    params:
        multitad_plot_type="multidataset",
    notebook:
        "../notebooks/Supplementaries_TADPlots.ipynb"


rule supplementary_tadplots_multiwindowsize:
    """Generate parameter sweep TAD figures across varied TopDom window sizes."""
    input:
        fname_data=RESULTS_DIR + "/results/final_data.csv.gz",
        sketch_hicfile=cool_input_for_source_wildcard(config["sketch"]["data_source"]),
        sketch_tadfile=(
            RESULTS_DIR + "/tads/data/tads.{source}.{tad_parameter}.csv"
        ).format(
            source=config["sketch"]["data_source"],
            tad_parameter=config["sketch"]["window_size"],
        ),
        tad_fname_list=expand(
            RESULTS_DIR + "/tads/data/tads.{source}.{tad_parameter}.csv",
            source=config["parameters"]["main_dataset"],
            tad_parameter=actual_window_size_list,
        ),
    output:
        outdir=report(
            directory(RESULTS_DIR + "/publication_figures/tad_plots_multiwindowsize/"),
            patterns=["{name}.pdf"],
            caption="../report/supplementary_plots.rst",
            category="Supplementary TAD Plots",
            subcategory="Multi-Window-Size",
        ),
    log:
        notebook=RESULTS_DIR + "/notebooks/Supplementaries_TADPlots.ipynb",
    conda:
        "../envs/python_stack.yaml"
    resources:
        mem_mb=lambda wildcards, attempt: (attempt + 2) * 10_000,
    params:
        multitad_plot_type="multiwindowsize",
    notebook:
        "../notebooks/Supplementaries_TADPlots.ipynb"
