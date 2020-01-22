#!/usr/bin/env bash

cd "$( dirname "${BASH_SOURCE[0]}" )"

# setup environment
rm -rf pipeline_test

mkdir -p pipeline_test/tads/
cp -v data/dummy_tads.csv pipeline_test/tads/tads.dummy.42.csv

mkdir -p pipeline_test/hic_files/plots/dummy/
cp -v data/empty_hic_heatmap.pdf pipeline_test/hic_files/plots/dummy/heatmap.chr1.pdf

# execute pipeline
cd ..
snakemake -pr --configfile "tests/config_dummy.yaml" --use-conda
