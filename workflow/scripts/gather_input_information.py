import cooler
import pandas as pd


def get_genome_version(samplename, genome_version_by_sample_name, cooler_info):
    implicit_genome_version = cooler_info["genome-assembly"]

    has_explicit_genome_version = samplename in genome_version_by_sample_name
    explicit_genome_version = genome_version_by_sample_name.get(samplename, "undefined")

    if implicit_genome_version == "unknown" and has_explicit_genome_version:
        return explicit_genome_version

    if implicit_genome_version != "unknown" and has_explicit_genome_version:
        print(
            f"Warning: Both implicit ({implicit_genome_version}) and explicit ({explicit_genome_version}) genome versions are available for sample {samplename}. Using explicit version."
        )
        return explicit_genome_version

    if implicit_genome_version == "unknown" and not has_explicit_genome_version:
        raise ValueError(
            f"Error: No genome version information available for sample {samplename}."
        )

    if implicit_genome_version != "unknown" and not has_explicit_genome_version:
        return implicit_genome_version

    raise RuntimeError("Unreachable code reached in get_genome_version.")


def main(fname_list, samplename_list, genome_version_by_sample_name, fname_out):
    tmp = []

    for samplename, fname in zip(samplename_list, fname_list):
        if fname.endswith(".mcool"):
            available_resolutions = cooler.fileops.list_coolers(fname)

            for resolution_id in available_resolutions:
                res = resolution_id.split("/")[-1]
                c = cooler.Cooler(f"{fname}::{resolution_id}")

                tmp.append(
                    {
                        "filename": fname,
                        "source": f"{samplename}__{res}",
                        "genome_assembly": get_genome_version(
                            samplename, genome_version_by_sample_name, c.info
                        ),
                        "bin_size": c.binsize,
                    }
                )
        else:
            c = cooler.Cooler(fname)

            tmp.append(
                {
                    "filename": fname,
                    "source": samplename,
                    "genome_assembly": get_genome_version(
                        samplename, genome_version_by_sample_name, c.info
                    ),
                    "bin_size": c.binsize,
                }
            )

    pd.DataFrame(tmp).to_csv(fname_out, index=False)


if __name__ == "__main__":
    main(
        snakemake.input.fname_list,
        snakemake.params["samplename_list"],
        snakemake.config.get("genome_versions", {}),
        snakemake.output.fname,
    )
