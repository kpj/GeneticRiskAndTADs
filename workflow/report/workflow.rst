Genetic Risk and Topologically-Associating Domains
==================================================

This workflow investigates the enrichment of disease-associated single nucleotide
polymorphisms (SNPs) from GWAS Catalog and DisGeNET within topologically-associating
domains (TADs) across human chromatin conformation capture (Hi-C / Micro-C) datasets.

Pipeline Structure
------------------
* **Hi-C Data Ingestion**: Extraction of contact count matrices from .cool/.mcool datasets.
* **TAD Calling**: Domain identification using TopDom across multiple window size parameters.
* **Database Aggregation**: Integration of SNP-disease annotations with Sequence Ontology and Disease Ontology mappings.
* **Enrichment Analysis**: Hypergeometric / binomial enrichment testing within TAD boundaries and borders.
* **Majority Voting**: Consensus aggregation across multiple TAD resolution parameters.
* **Visualization**: Generation of contact matrix heatmaps, enrichment distributions, and publication figures.
