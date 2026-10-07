"""TopDom TAD caller adapter."""

import os
import tempfile
from pathlib import Path
import pandas as pd
import sh


def run_topdom(
    fname_matrix: str,
    bin_size: int,
    chromosome: str,
    window_size: int,
    fname_out: str,
) -> None:
    """Execute TopDom boundary caller and format domain predictions."""
    chrom_str = f"chr{chromosome}" if not str(chromosome).startswith("chr") else str(chromosome)

    print(f"[TopDom] Preparing input matrix for {chrom_str} (bin_size={bin_size}, w={window_size})")
    df_count = pd.read_csv(fname_matrix, index_col=0)

    # Convert to TopDom compatible format (chr, start, end, ... contacts ...)
    df_count.insert(0, "chr", chrom_str)
    df_count.insert(1, "td_start", df_count.index)
    df_count.insert(2, "td_end", df_count.index + bin_size)

    with tempfile.TemporaryDirectory() as tmpdir:
        input_tsv = os.path.join(tmpdir, "topdom_input.tsv")
        out_prefix = os.path.join(tmpdir, "topdom_out")
        out_bed = f"{out_prefix}.bed"

        df_count.to_csv(input_tsv, sep="\t", index=False, header=False)

        print(f"[TopDom] Running TopDom R script with window_size={window_size}")
        cmd = f"""
            TopDom::TopDom('{input_tsv}', {window_size}, outFile='{out_prefix}', debug=FALSE)
        """
        sh.Rscript("--vanilla", "-e", cmd, _fg=True)

        if not os.path.exists(out_bed):
            raise RuntimeError(f"[TopDom] Expected output file not found: {out_bed}")

        df_topdom = pd.read_csv(
            out_bed,
            sep="\t",
            header=None,
            names=["chrname", "tad_start", "tad_stop", "type"],
        )

        df_topdom = df_topdom[df_topdom["type"] == "domain"].copy()
        df_topdom.drop("type", axis=1, inplace=True)

        df_topdom["tad_start"] = df_topdom["tad_start"].astype(int)
        df_topdom["tad_stop"] = df_topdom["tad_stop"].astype(int)

        Path(fname_out).parent.mkdir(parents=True, exist_ok=True)
        df_topdom.to_csv(fname_out, index=False)
        print(f"[TopDom] Saved {len(df_topdom)} domains to {fname_out}")
