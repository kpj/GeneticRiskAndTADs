import re
import cooler


# Copied over from "rules/common.smk".
def parse_source_wildcard(source: str):
    """Convert source wildcard value into (base, resolution) tuple."""
    m = re.match(r"^(?P<base>.+)__((?P<res>\d+))$", source)
    if m:
        return m.group("base"), int(m.group("res"))
    return source, None


def main(input_mcool, source_wildcard, output_cool):
    _, res = parse_source_wildcard(source_wildcard)
    resolution_id = f"/resolutions/{res}"

    available_resolutions = cooler.fileops.list_coolers(input_mcool)
    assert resolution_id in available_resolutions, (
        f"{resolution_id} not in {available_resolutions}"
    )

    in_path = f"{input_mcool}::{resolution_id}"
    cooler.fileops.cp(in_path, output_cool)


if __name__ == "__main__":
    main(
        snakemake.input.mcool,
        snakemake.wildcards.source,
        snakemake.output.cool,
    )
