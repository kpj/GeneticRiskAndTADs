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


def main():
    # Potentially install TopDom (it's not available as a conda package).
    maybe_install_r_package("TopDom")
    # Install SpectralTAD (bioconductor-spectraltad is not solvable on osx-arm64 conda)
    maybe_install_bioc_package("SpectralTAD")


if __name__ == "__main__":
    main()
