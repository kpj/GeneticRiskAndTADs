import cooler
import pandas as pd


def main(fname_list, samplename_list, fname_out):
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
                        "genome_assembly": c.info["genome-assembly"],
                        "bin_size": c.binsize,
                    }
                )
        else:
            c = cooler.Cooler(fname)

            tmp.append(
                {
                    "filename": fname,
                    "source": samplename,
                    "genome_assembly": c.info["genome-assembly"],
                    "bin_size": c.binsize,
                }
            )

    pd.DataFrame(tmp).to_csv(fname_out, index=False)


if __name__ == "__main__":
    main(
        snakemake.input.fname_list,
        snakemake.params["samplename_list"],
        snakemake.output.fname,
    )
