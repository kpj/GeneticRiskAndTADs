# TADs

## Workflow

Force execution as follows:
```bash
$ snakemake -pr -R convert_tad_coordinates
```

or run `python execute.py` to auto-run different configurations.


## Important notes

Install required dependencies with `pip install -r requirements.txt`.
Additional R dependencies:
    * `biomaRt`

Files in `./cache` are not automatically recreated.
