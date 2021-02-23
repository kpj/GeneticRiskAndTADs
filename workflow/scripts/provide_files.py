import sh


def main(fname_pair_list):
    for source_url, target_url in fname_pair_list:
        sh.zstd('-d', source_url, '-o', target_url)


if __name__ == '__main__':
    assert len(snakemake.input) == len(snakemake.output)
    main(zip(snakemake.input, snakemake.output))
