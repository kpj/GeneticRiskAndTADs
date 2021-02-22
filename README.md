# Genetic Risks and Topologically-Associating Domains

A Snakemake workflow for the investigation of disease-associated SNP enrichments in topologically-associating domains.


## Usage

Run pipeline as follows:
```bash
$ snakemake -pr -j 1 --use-conda
```


## Tests

Run tests as follows:
```bash
$ ./tests/run.sh
```


# Data sources

* Hi-C data: ftp://cooler.csail.mit.edu/coolers/hg19/
* More data sources can be found in `config/config.yaml`
