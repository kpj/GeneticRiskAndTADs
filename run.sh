#!/usr/bin/env bash

N_JOBS=$(uv run python -c "import os; print(max(1, (os.cpu_count() or 1) - 2))")
uv run snakemake --jobs "$N_JOBS" --software-deployment-method conda --resources hdf5_lock=1 notebook_slots=1 "$@"
