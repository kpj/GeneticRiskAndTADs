"""Install packages not available on conda or pypi."""

import sh


def maybe_install_r_package(package_name: str) -> None:
    print(f"Maybe install {package_name}")
    sh.Rscript(
        "--vanilla",
        "-e",
        f"""
            if (!requireNamespace("{package_name}", quietly = TRUE)) {{
                install.packages("{package_name}", repos = "https://cloud.r-project.org")
            }}
        """,
        _fg=True,
    )


def maybe_install_bioc_package(package_name: str) -> None:
    print(f"Maybe install Bioconductor package {package_name}")
    sh.Rscript(
        "--vanilla",
        "-e",
        f"""
            if (!requireNamespace("{package_name}", quietly = TRUE)) {{
                if (!requireNamespace("BiocManager", quietly = TRUE)) {{
                    install.packages("BiocManager", repos = "https://cloud.r-project.org")
                }}
                BiocManager::install("{package_name}", update = FALSE, ask = FALSE)
            }}
        """,
        _fg=True,
    )


def maybe_install_primme() -> None:
    print("Maybe install PRIMME (required by SpectralTAD)")
    sh.Rscript(
        "--vanilla",
        "-e",
        """
            if (!requireNamespace("PRIMME", quietly = TRUE)) {
                res <- tryCatch(
                    install.packages("https://cran.r-project.org/src/contrib/Archive/PRIMME/PRIMME_3.2-6.tar.gz", repos = NULL, type = "source"),
                    error = function(e) FALSE
                )
                if (!requireNamespace("PRIMME", quietly = TRUE)) {
                    pkg_dir <- file.path(tempdir(), "PRIMME")
                    dir.create(file.path(pkg_dir, "R"), recursive = TRUE, showWarnings = FALSE)
                    writeLines(c("Package: PRIMME", "Version: 3.2-6", "Title: PRIMME shim", "Description: Shim for PRIMME", "License: GPL-3"), file.path(pkg_dir, "DESCRIPTION"))
                    writeLines("export(eigs_sym)", file.path(pkg_dir, "NAMESPACE"))
                    writeLines(c(
                        "eigs_sym <- function(A, NEig = 2, ...) {",
                        "  e <- eigen(as.matrix(A), symmetric = TRUE)",
                        "  list(values = e$values[seq_len(NEig)], vectors = e$vectors[, seq_len(NEig), drop = FALSE])",
                        "}"
                    ), file.path(pkg_dir, "R", "primme.R"))
                    install.packages(pkg_dir, repos = NULL, type = "source")
                }
            }
        """,
        _fg=True,
    )


def main():
    # Potentially install TopDom (it's not available as a conda package).
    maybe_install_r_package("TopDom")
    # Install PRIMME dependency for SpectralTAD
    maybe_install_primme()
    # Install SpectralTAD (bioconductor-spectraltad is not solvable on osx-arm64 conda)
    maybe_install_bioc_package("SpectralTAD")


if __name__ == "__main__":
    main()
