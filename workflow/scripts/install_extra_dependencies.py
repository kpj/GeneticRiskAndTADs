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
    )


def main():
    # Potentially install TopDom (it's not available as a conda package).
    maybe_install_r_package("TopDom")


if __name__ == "__main__":
    main()
