# TADs

## Workflow

Force execution as follows:
```bash
$ snakemake -pr -R convert_tad_coordinates
```

or run `python execute.py` to auto-run different configurations.


## Important notes

Files that are not automatically recreated:
* cache/doid_graph.edgelist.gz
* cache/efo_graph.edgelist.gz
* cache/snp_annotations.json.gz
