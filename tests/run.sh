#!/usr/bin/env bash
set -e

# ensure clean conda environment state for nested snakemake activations
unset CONDA_SHLVL CONDA_PREFIX CONDA_DEFAULT_ENV

cd "$( dirname "${BASH_SOURCE[0]}" )"

# parse arguments
PLOT=false
SNAKEMAKE_ARGS=()
for arg in "$@"; do
    if [ "$arg" == "--plot" ]; then
        PLOT=true
    else
        SNAKEMAKE_ARGS+=("$arg")
    fi
done

# setup environment
cd ..

# plot pipeline overview(s)
if [ "$PLOT" = true ]; then
    for graph_type in dag rulegraph filegraph; do
        uv run snakemake --config experiment=test_dummy --forceall --$graph_type | dot -Tpdf > "tests/test_$graph_type.pdf"
    done
fi

# execute pipeline
uv run snakemake --config experiment=test_dummy --jobs 1 --software-deployment-method conda --resources hdf5_lock=1 "${SNAKEMAKE_ARGS[@]}"

# run output assertions and schema validations
uv run pytest tests/

