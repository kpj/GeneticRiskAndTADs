#!/usr/bin/env bash

cd "$( dirname "${BASH_SOURCE[0]}" )"

# setup environment
rm -rf pipeline_test
cd ..

# plot pipeline overview
snakemake --config configfile="tests/config_dummy.yaml" --forceall --dag | dot -Tpdf > tests/dag_test.pdf

# execute pipeline
snakemake -pr --config configfile="tests/config_dummy.yaml" -j 1 --use-conda "$@"
