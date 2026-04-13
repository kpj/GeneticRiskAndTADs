#!/usr/bin/env bash

cd "$( dirname "${BASH_SOURCE[0]}" )"

# parse arguments
PLOT=false
for arg in "$@"; do
    if [ "$arg" == "--plot" ]; then
        PLOT=true
        shift
    fi
done

# setup environment
cd ..

# plot pipeline overview(s)
if [ "$PLOT" = true ]; then
    for graph_type in dag rulegraph filegraph; do
        uv run snakemake --config configfile="tests/config_dummy.yaml" --forceall --$graph_type | dot -Tpdf > "tests/test_$graph_type.pdf"
    done
fi

# execute pipeline
uv run snakemake --config configfile="tests/config_dummy.yaml" --jobs 1 --software-deployment-method conda --resources hdf5_lock=1 "$@"
