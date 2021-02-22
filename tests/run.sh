#!/usr/bin/env bash

cd "$( dirname "${BASH_SOURCE[0]}" )"

# setup environment
rm -rf pipeline_test
cd ..

# plot pipeline overview(s)
for graph_type in dag rulegraph filegraph; do
    snakemake --config configfile="tests/config_dummy.yaml" --forceall --$graph_type | dot -Tpdf > "tests/test_$graph_type.pdf"
done

# execute pipeline
snakemake -pr --config configfile="tests/config_dummy.yaml" -j 1 --use-conda "$@"
