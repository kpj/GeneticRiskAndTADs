#!/usr/bin/env bash

cd "$( dirname "${BASH_SOURCE[0]}" )"

# setup environment
cd ..

# plot pipeline overview(s)
for graph_type in dag rulegraph filegraph; do
    snakemake --config configfile="tests/config_dummy.yaml" --forceall --$graph_type | dot -Tpdf > "tests/test_$graph_type.pdf"
done

# execute pipeline
snakemake -p --config configfile="tests/config_dummy.yaml" -j 1 --software-deployment-method conda "$@"
