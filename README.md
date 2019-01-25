# TADs

## Workflow

Force execution as follows:
```bash
$ snakemake -pr -R convert_tad_coordinates
```

Post-processing of multiple pipeline runs (generated using `project_manager`) may be done using `MultiRunPostAnalysis.ipynb`.


## Important notes

Install required dependencies with `pip install -r requirements.txt`.
Additional R dependencies:
    * `biomaRt`

Files in `./cache` are not automatically recreated.
