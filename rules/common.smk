import os
from urllib.parse import urlparse

from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider

HTTP = HTTPRemoteProvider()
FTP = FTPRemoteProvider()


def url_wrapper(url):
    if os.path.isfile(srcdir(url)):
        # is local
        return srcdir(url)
    else:
        # is remote
        o = urlparse(url)
        if o.scheme in ('http', 'https'):
            return HTTP.remote(url)
        elif o.scheme == 'ftp':
            return FTP.remote(url)
        else:
            raise RuntimeError(f'Invalid url: "{url}"')
