import os
from urllib.parse import urlparse

from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider

HTTP = HTTPRemoteProvider()
FTP = FTPRemoteProvider()


def url_wrapper(url, remote_kwargs={'keep_local': True}):
    src_url = os.path.join(workflow.basedir, url)
    if src_url is not None and os.path.isfile(src_url):
        # is local
        return src_url
    else:
        # is remote
        o = urlparse(url)
        if o.scheme in ('http', 'https'):
            return HTTP.remote(url, **remote_kwargs)
        elif o.scheme == 'ftp':
            return FTP.remote(url, **remote_kwargs)
        else:
            raise RuntimeError(f'Invalid url: "{url}"')
