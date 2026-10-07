"""cooltools insulation caller adapter."""

from pathlib import Path

import cooler
import cooltools
import pandas as pd


def run_cooltools(
    fname_cool: str,
    chromosome: str,
    window_bp: int,
    fname_out: str,
) -> None:
    """Execute cooltools diamond insulation and extract domain intervals."""
    clr = cooler.Cooler(fname_cool)

    chrom_str = str(chromosome)
    if chrom_str not in clr.chromnames:
        if f"chr{chromosome}" in clr.chromnames:
            chrom_str = f"chr{chromosome}"
        elif chromosome.startswith("chr") and chromosome[3:] in clr.chromnames:
            chrom_str = chromosome[3:]
        else:
            raise ValueError(
                f"Chromosome '{chromosome}' not found in cooler: {clr.chromnames}"
            )

    window_int = int(window_bp)
    print(
        f"[cooltools] Computing insulation score for {chrom_str} (window_bp={window_int})"
    )

    try:
        ins_df = cooltools.insulation(clr, [window_int])
    except ValueError:
        ins_df = cooltools.insulation(clr, [window_int], ignore_diags=2)
    ins_chrom = ins_df[ins_df["chrom"] == chrom_str].copy()

    boundary_col = f"is_boundary_{window_int}"
    if boundary_col not in ins_chrom.columns:
        matching_cols = [c for c in ins_chrom.columns if c.startswith("is_boundary_")]
        if matching_cols:
            boundary_col = matching_cols[0]
        else:
            raise KeyError(
                f"Expected boundary column '{boundary_col}' in cooltools output: {ins_chrom.columns}"
            )

    boundaries = (
        ins_chrom[ins_chrom[boundary_col]].sort_values("start").reset_index(drop=True)
    )

    domains = []
    out_chrom = (
        f"chr{chromosome}" if not str(chromosome).startswith("chr") else str(chromosome)
    )
    for i in range(len(boundaries) - 1):
        d_start = int(boundaries.loc[i, "end"])
        d_stop = int(boundaries.loc[i + 1, "start"])
        if d_stop > d_start:
            domains.append(
                {
                    "chrname": out_chrom,
                    "tad_start": d_start,
                    "tad_stop": d_stop,
                }
            )

    df_domains = pd.DataFrame(domains)
    if df_domains.empty:
        df_domains = pd.DataFrame(columns=["chrname", "tad_start", "tad_stop"])

    Path(fname_out).parent.mkdir(parents=True, exist_ok=True)
    df_domains.to_csv(fname_out, index=False)
    print(f"[cooltools] Saved {len(df_domains)} domains to {fname_out}")
