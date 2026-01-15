import os
import sys
from urllib.parse import urlparse


def url_wrapper(url, remote_kwargs={"keep_local": True}, use_basedir=False):
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
