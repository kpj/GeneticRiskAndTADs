#!/usr/bin/env bash

cd "$( dirname "${BASH_SOURCE[0]}" )"

# setup environment
rm -rf pipeline_test

# execute pipeline
cd ..
snakemake -pr --configfile "tests/config_dummy.yaml" --use-conda
