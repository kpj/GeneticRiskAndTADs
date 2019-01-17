cd "$( dirname "${BASH_SOURCE[0]}" )"

rm -rf pipeline_run
mkdir -p pipeline_run/cache
cp -v data/snp_annotations.csv pipeline_run/cache

cd ..
SNAKEMAKE__CONFIG_FILE="test_case/config_dummy.yaml" snakemake -pr
