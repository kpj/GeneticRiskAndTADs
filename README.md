# TADs

## Workflow

* Setup
    * `ConvertTADGenomicCoordinates`: convert TAD-coordinates from hg18 to hg38
    * `LoadDisGeNET`: assemble custom DisGeNET (TAD associations, etc)

* Analyses
    * `TAD_Borders`: Compute network-coherences
    * `ComputeTADEnrichments`: compute TAD (-border) enrichments

Force execution as follows:
```bash
$ snakemake -pr -R convert_tad_coordinates
```
