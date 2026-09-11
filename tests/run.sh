#!/usr/bin/env bash

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

# verify that post-cutoff SNP (rs999) was filtered out
if [ -f "results/test_dummy/databases/initial.csv" ]; then
    uv run python -c "
import csv
with open('results/test_dummy/databases/initial.csv') as f:
    reader = csv.DictReader(f)
    snp_ids = {row['snpId'] for row in reader}
assert 'rs999' not in snp_ids, 'rs999 should have been filtered out by gwas_date_range!'
print('Test assertion passed: rs999 was correctly filtered out by gwas_date_range.')
"
fi
