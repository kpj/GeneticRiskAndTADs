rule provide_input_files:
    """Stage and decompress input association databases and ontologies into experiment workspace."""
    input:
        disgenet_fname=url_wrapper(config["input_files"]["disgenet"]),
        gwascatalog_fname=url_wrapper(config["input_files"]["gwascatalog"]),
        efo_fname=url_wrapper(config["input_files"]["exp_factor_ontology"]),
        so_fname=url_wrapper(config["input_files"]["sequence_ontology"]),
        dbsnp_hg19_fname=url_wrapper(config["input_files"]["dbsnp_hg19"]),
    output:
        disgenet_fname=RESULTS_DIR + "/input/disgenet.tsv",
        gwascatalog_fname=RESULTS_DIR + "/input/gwas_catalog.tsv",
        efo_fname=RESULTS_DIR + "/input/efo.owl",
        so_fname=RESULTS_DIR + "/input/so.owl",
        dbsnp_hg19_fname=RESULTS_DIR + "/input/dbSnp155Common_hg19.bb",
    conda:
        "../envs/python_stack.yaml"
    script:
        "../scripts/provide_files.py"


rule assemble_input_databases:
    """Harmonize variant-trait associations across GWAS Catalog, DisGeNET, and ontologies."""
    input:
        disgenet_fname=RESULTS_DIR + "/input/disgenet.tsv",
        gwascatalog_fname=RESULTS_DIR + "/input/gwas_catalog.tsv",
        efo_fname=RESULTS_DIR + "/input/efo.owl",
        so_fname=RESULTS_DIR + "/input/so.owl",
        dbsnp_hg19_fname=RESULTS_DIR + "/input/dbSnp155Common_hg19.bb",
    output:
        db_fname=RESULTS_DIR + "/databases/initial.csv",
        raw_veps=RESULTS_DIR + "/databases/vep.csv",
    log:
        notebook=RESULTS_DIR + "/notebooks/AssembleInputDatabases.ipynb",
    conda:
        "../envs/python_stack.yaml"
    resources:
        mem_mb=3_000,
    notebook:
        "../notebooks/AssembleInputDatabases.ipynb"
