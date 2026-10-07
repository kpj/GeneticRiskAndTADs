"""ConsensusTAD aggregator using Weighted Interval Scheduling (WIS)."""

from bisect import bisect_right
from pathlib import Path
from typing import List


def weighted_interval_scheduling(intervals: List[dict]) -> List[dict]:
    """Find the optimal non-overlapping subset of intervals maximizing total weight."""
    if not intervals:
        return []

    # Sort intervals by end position
    sorted_intervals = sorted(intervals, key=lambda x: (x["tad_stop"], x["tad_start"]))
    n = len(sorted_intervals)

    # Precompute p(j): the rightmost interval i < j that does not overlap with j
    end_times = [iv["tad_stop"] for iv in sorted_intervals]
    p = []
    for j in range(n):
        start_j = sorted_intervals[j]["tad_start"]
        idx = bisect_right(end_times, start_j) - 1
        p.append(idx)

    # Dynamic programming table
    opt = [0.0] * (n + 1)
    for j in range(1, n + 1):
        weight_j = sorted_intervals[j - 1]["weight"]
        p_j = p[j - 1]
        weight_with_j = weight_j + (opt[p_j + 1] if p_j >= 0 else 0.0)
        weight_without_j = opt[j - 1]
        opt[j] = max(weight_without_j, weight_with_j)

    # Backtrack to reconstruct the optimal solution
    selected = []
    curr = n
    while curr > 0:
        weight_curr = sorted_intervals[curr - 1]["weight"]
        p_curr = p[curr - 1]
        weight_with = weight_curr + (opt[p_curr + 1] if p_curr >= 0 else 0.0)
        if weight_with >= opt[curr - 1]:
            selected.append(sorted_intervals[curr - 1])
            curr = p_curr + 1
        else:
            curr = curr - 1

    selected.reverse()
    return selected


def run_consensustad(
    tad_files: List[str],
    fname_out: str,
) -> None:
    """Aggregate predictions across multiple caller outputs using WIS."""
    import pandas as pd

    print(f"[ConsensusTAD] Aggregating {len(tad_files)} TAD caller files")
    all_intervals = []
    for f in tad_files:
        if not Path(f).exists():
            continue
        df = pd.read_csv(f)
        for _, row in df.iterrows():
            all_intervals.append(
                {
                    "chrname": str(row["chrname"]),
                    "tad_start": int(row["tad_start"]),
                    "tad_stop": int(row["tad_stop"]),
                    "weight": 1.0,
                }
            )

    df_all = pd.DataFrame(all_intervals)
    final_domains = []

    if not df_all.empty:
        for chrom, group in df_all.groupby("chrname"):
            grouped_intervals = (
                group.groupby(["tad_start", "tad_stop"])
                .size()
                .reset_index(name="recurrence")
            )
            chrom_intervals = []
            for _, r in grouped_intervals.iterrows():
                chrom_intervals.append(
                    {
                        "chrname": chrom,
                        "tad_start": int(r["tad_start"]),
                        "tad_stop": int(r["tad_stop"]),
                        "weight": float(r["recurrence"]),
                    }
                )

            chosen = weighted_interval_scheduling(chrom_intervals)
            final_domains.extend(chosen)

    df_out = pd.DataFrame(final_domains)
    if not df_out.empty:
        df_out = (
            df_out[["chrname", "tad_start", "tad_stop"]]
            .sort_values(["chrname", "tad_start"])
            .reset_index(drop=True)
        )
    else:
        df_out = pd.DataFrame(columns=["chrname", "tad_start", "tad_stop"])

    Path(fname_out).parent.mkdir(parents=True, exist_ok=True)
    df_out.to_csv(fname_out, index=False)
    print(f"[ConsensusTAD] Saved {len(df_out)} consensus domains to {fname_out}")
