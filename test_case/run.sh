#!/usr/bin/env bash

cd "$( dirname "${BASH_SOURCE[0]}" )"

# setup environment
rm -rf pipeline_test

mkdir -p pipeline_test/tads/
cp -v data/dummy_tads.csv pipeline_test/tads/tads.dummy.42.csv

# execute pipeline
cd ..
snakemake -pr --configfile "test_case/config_dummy.yaml"
