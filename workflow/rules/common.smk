import os
import sys
from urllib.parse import urlparse


def url_wrapper(url, remote_kwargs={"keep_local": True}, use_basedir=False):
    """Wrap URL to local or remote storage as appropriate."""
    if use_basedir:
        src_url = os.path.join(workflow.basedir, url)
    else:
        src_url = url

    if src_url is not None and os.path.isfile(src_url):
        # is local
        return src_url
    else:
        # is remote
        o = urlparse(url)
        if o.scheme in ("http", "https"):
            return storage.http(url)
        elif o.scheme == "ftp":
            return storage.ftp(url)
        else:
            print(
                f'Invalid url "{url}", returning input without transformation',
                file=sys.stderr,
            )
            return url


def is_mcool(path: str) -> bool:
    """Check if the given path corresponds to a .mcool file."""
    return str(path).endswith(".mcool")


def generate_data_source_list(input_samples, mcool_resolutions_by_input):
    """Generate list of data source wildcards based on input samples and resolutions."""
    data_source_list = []
    for base in input_samples.keys():
        src = input_samples[base]
        if is_mcool(src):
            for res in mcool_resolutions_by_input[base]:
                data_source_list.append(f"{base}__{res}")
        else:
            data_source_list.append(base)
    return data_source_list


def sample_url_from_source_wildcard(source: str):
    """Get path to retrieve original sample from for given source wildcard."""
    base, _ = parse_source_wildcard(source)
    return url_wrapper(config["samples"][base])


def parse_source_wildcard(source: str):
    """Convert source wildcard value into (base, resolution) tuple."""
    m = re.match(r"^(?P<base>.+)__((?P<res>\d+))$", source)
    if m:
        return m.group("base"), int(m.group("res"))
    return source, None


def cool_input_for_source_wildcard(source: str):
    """Get path to (potentially processed) cool input file for given source wildcard."""
    base, res = parse_source_wildcard(source)
    src = config["samples"][base]

    if is_mcool(src):
        if res is None:
            raise ValueError(
                f"Sample '{base}' is .mcool but '{source}' has no '__<resolution>' suffix."
            )
        return f"results/hic_files/cool/{source}.cool"
    else:
        return url_wrapper(src)
