"""Dispatcher script for chromosome-level TAD calling across modular callers."""

import sys
from pathlib import Path

# Add callers directory to sys.path
script_dir = Path(__file__).resolve().parent
if str(script_dir) not in sys.path:
    sys.path.insert(0, str(script_dir))

import pandas as pd
from callers import parse_caller_config
from callers.topdom import run_topdom
from callers.spectraltad import run_spectraltad
from callers.cooltools import run_cooltools


def main():
    caller_config = getattr(
        snakemake.wildcards,
        "caller_config",
        getattr(snakemake.wildcards, "tad_parameter", None),
    )
    if not caller_config:
        raise ValueError("Missing caller configuration wildcard (caller_config or tad_parameter)")

    parsed = parse_caller_config(caller_config)
    caller = parsed["caller"]
    source = snakemake.wildcards.source
    chromosome = snakemake.wildcards.chromosome
    fname_matrix = snakemake.input.fname
    fname_info = snakemake.input.fname_info
    fname_out = snakemake.output.fname

    df_info = pd.read_csv(fname_info, index_col=1)
    bin_size = int(df_info.loc[source, "bin_size"])

    print(
        f"[compute_tads] Running caller '{caller}' with config '{caller_config}' "
        f"on source '{source}', chr{chromosome}"
    )

    if caller == "topdom":
        run_topdom(
            fname_matrix=fname_matrix,
            bin_size=bin_size,
            chromosome=chromosome,
            window_size=parsed["window_size"],
            fname_out=fname_out,
        )
    elif caller == "spectraltad":
        run_spectraltad(
            fname_matrix=fname_matrix,
            bin_size=bin_size,
            chromosome=chromosome,
            levels=parsed["levels"],
            fname_out=fname_out,
        )
    elif caller == "cooltools":
        fname_cool = getattr(snakemake.input, "cool", None)
        if not fname_cool:
            raise ValueError(
                "Cooler input required for cooltools caller but not provided in snakemake.input.cool"
            )
        run_cooltools(
            fname_cool=fname_cool,
            chromosome=chromosome,
            window_bp=parsed["window_bp"],
            fname_out=fname_out,
        )
    else:
        raise ValueError(
            f"Unsupported TAD caller '{caller}' (parsed from caller_config='{caller_config}')"
        )


if __name__ == "__main__":
    main()
