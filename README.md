# TADs

## Workflow

* Setup
    * `ConvertTADGenomicCoordinates`: convert TAD-coordinates from hg19 to hg38
    * `LoadDisGeNET`: assemble custom DisGeNET (TAD associations, etc)

* Analyses
    * `TAD_Borders`: Compute network-coherences
    * `ComputeTADEnrichments`: compute TAD (-border) enrichments

Execute as follows:
```bash
$ snakemake -pr
```
