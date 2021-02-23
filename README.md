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
* [GWAS Catalog](https://www.ebi.ac.uk/gwas/home): `v1.0.2-associations_e100_r2021-02-10`
* [DisGeNET](https://www.disgenet.org/): `DisGeNET 2020 - v7.0`
* [HumanDiseaseOntology](https://disease-ontology.org/): `r2021-01-28`
* [Experimental Factor Ontology](https://www.ebi.ac.uk/efo/): `v3.27.0`
* [SequenceOntology](http://www.sequenceontology.org/): `v3.1`
