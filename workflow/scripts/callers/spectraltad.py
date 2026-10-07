"""SpectralTAD caller adapter."""

import os
import tempfile
from pathlib import Path

import pandas as pd
import sh


def run_spectraltad(
    fname_matrix: str,
    bin_size: int,
    chromosome: str,
    levels: int,
    fname_out: str,
) -> None:
    """Execute SpectralTAD clustering caller and extract disjoint Level 1 domains."""
    chrom_str = (
        f"chr{chromosome}" if not str(chromosome).startswith("chr") else str(chromosome)
    )

    print(
        f"[SpectralTAD] Preparing input matrix for {chrom_str} (bin_size={bin_size}, levels={levels})"
    )
    df_count = pd.read_csv(fname_matrix, index_col=0)

    # Convert to n x (n+3) format: chr, start, end, ... contacts ...
    df_count.insert(0, "chr", chrom_str)
    df_count.insert(1, "td_start", df_count.index)
    df_count.insert(2, "td_end", df_count.index + bin_size)

    with tempfile.TemporaryDirectory() as tmpdir:
        input_tsv = os.path.join(tmpdir, "spectraltad_input.tsv")
        out_tsv = os.path.join(tmpdir, "spectraltad_out.tsv")

        df_count.to_csv(input_tsv, sep="\t", index=False, header=False)

        print(f"[SpectralTAD] Running SpectralTAD R script with levels={levels}")
        cmd = f"""
            suppressPackageStartupMessages(library(SpectralTAD))
            df_mat <- read.table('{input_tsv}', sep='\\t', header=FALSE)
            res <- tryCatch(
                SpectralTAD(df_mat, chr = '{chrom_str}', levels = {levels}, qual_filter = FALSE),
                error = function(e) {{
                    warning(paste("SpectralTAD error:", e$message))
                    data.frame(chr = character(0), start = integer(0), end = integer(0))
                }}
            )
            if (is.data.frame(res)) {{
                df_out <- res
            }} else if (is.list(res) && length(res) > 0) {{
                df_out <- res[[1]]
            }} else {{
                df_out <- data.frame(chr=character(0), start=integer(0), end=integer(0))
            }}
            write.table(df_out, file='{out_tsv}', sep='\\t', row.names=FALSE, col.names=TRUE, quote=FALSE)
        """
        sh.Rscript("--vanilla", "-e", cmd, _fg=True)

        if not os.path.exists(out_tsv):
            raise RuntimeError(
                f"[SpectralTAD] Expected output file not found: {out_tsv}"
            )

        df_spectral = pd.read_csv(out_tsv, sep="\t")

        # Standardize columns: find chr, start, end
        col_map = {}
        for col in df_spectral.columns:
            cl = str(col).lower()
            if cl in ("chr", "chrname", "chrom", "chromosome"):
                col_map[col] = "chrname"
            elif cl in ("start", "tad_start", "td_start"):
                col_map[col] = "tad_start"
            elif cl in ("end", "stop", "tad_stop", "td_end"):
                col_map[col] = "tad_stop"

        if len(col_map) == 3:
            df_spectral = df_spectral[list(col_map.keys())].rename(columns=col_map)
        else:
            if df_spectral.shape[1] >= 3:
                df_spectral = df_spectral.iloc[:, :3]
                df_spectral.columns = ["chrname", "tad_start", "tad_stop"]
            else:
                df_spectral = pd.DataFrame(columns=["chrname", "tad_start", "tad_stop"])

        if not df_spectral.empty:
            df_spectral["chrname"] = chrom_str
            df_spectral["tad_start"] = df_spectral["tad_start"].astype(int)
            df_spectral["tad_stop"] = df_spectral["tad_stop"].astype(int)
            df_spectral = df_spectral[
                df_spectral["tad_stop"] > df_spectral["tad_start"]
            ].copy()
            df_spectral = df_spectral.sort_values("tad_start").reset_index(drop=True)

        Path(fname_out).parent.mkdir(parents=True, exist_ok=True)
        df_spectral.to_csv(fname_out, index=False)
        print(f"[SpectralTAD] Saved {len(df_spectral)} domains to {fname_out}")
