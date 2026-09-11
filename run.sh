#!/usr/bin/env bash
set -euo pipefail

if [[ $# -eq 0 || "$1" =~ ^- ]]; then
    echo "Error: An experiment name must be provided as the first argument." >&2
    echo "Available experiments in config/experiments/:" >&2
    for f in config/experiments/*.yaml; do
        echo "  - $(basename "$f" .yaml)" >&2
    done
    echo "" >&2
    echo "Usage: ./run.sh <experiment> [snakemake options...]" >&2
    exit 1
fi

EXPERIMENT="$1"
shift

if [[ ! -f "config/experiments/${EXPERIMENT}.yaml" ]]; then
    echo "Error: Experiment config 'config/experiments/${EXPERIMENT}.yaml' not found." >&2
    echo "Available experiments in config/experiments/:" >&2
    for f in config/experiments/*.yaml; do
        echo "  - $(basename "$f" .yaml)" >&2
    done
    exit 1
fi

N_JOBS=$(uv run python -c "import os; print(max(1, (os.cpu_count() or 1) - 2))")
uv run snakemake --jobs "$N_JOBS" --software-deployment-method conda --resources hdf5_lock=1 notebook_slots=1 --config "experiment=$EXPERIMENT" "$@"
