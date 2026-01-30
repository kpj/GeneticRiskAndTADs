# Genetic Risks and Topologically-Associating Domains

A Snakemake workflow for the investigation of disease-associated SNP enrichments in topologically-associating domains.

## Usage

Run pipeline as follows:

```bash
snakemake --jobs 1 --software-deployment-method conda --resources hdf5_lock=1
```

## Tests

Run tests as follows:

```bash
./tests/run.sh
```

## Data sources

* Hi-C data: <https://data.4dnucleome.org/hic-data-overview>, <https://usgs2.osn.mghpcc.org/cooler01/index.html>
* [GWAS Catalog](https://www.ebi.ac.uk/gwas/home)
* [DisGeNET](https://www.disgenet.org/)
* [HumanDiseaseOntology](https://disease-ontology.org/)
* [Experimental Factor Ontology](https://www.ebi.ac.uk/efo/)
* [SequenceOntology](http://www.sequenceontology.org/)
