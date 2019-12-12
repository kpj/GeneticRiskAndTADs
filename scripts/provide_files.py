import shutil
from pathlib import Path

import sh


def main(fname_pair_list):
    for source_url, target_url in fname_pair_list:
        source_url = Path(source_url)
        target_url = Path(target_url)

        # copy file
        shutil.copy(source_url, target_url)

        # handle format conversions
        source_suffix = source_url.suffix
        target_suffix = target_url.suffix
        if source_suffix != target_suffix:
            if target_suffix == '.gz':
                # remove suffix
                new_fname = target_url.parent / target_url.name[:-len('.gz')]
                target_url.rename(new_fname)

                # gzip
                sh.gzip(new_fname)
            elif source_suffix == '.gz':
                # add suffix
                new_fname = target_url.parent / f'{target_url.name}.gz'
                target_url.rename(new_fname)

                # gunzip
                sh.gunzip(new_fname)
            else:
                raise RuntimeError(
                    'Unsupported suffixes: '
                    f'"{source_suffix}", "{target_suffix}"')


if __name__ == '__main__':
    assert len(snakemake.input) == len(snakemake.output)
    main(zip(snakemake.input, snakemake.output))
