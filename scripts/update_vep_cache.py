import io
import os
import time

import pandas as pd
import pybiomart
import yaml
import zstandard as zstd
from tqdm import tqdm

config_path = "config/config.yaml"
if not os.path.exists(config_path):
    raise FileNotFoundError(
        f"Could not find config at {config_path}. Please run this script from the project root."
    )

with open(config_path) as f:
    config = yaml.safe_load(f)

MAX_RETRIES = 5
BACKOFF_FACTOR = 2


def format_error(e, max_len=150):
    err_str = str(e)
    return err_str[:max_len] + "..." if len(err_str) > max_len else err_str


def load_snps(path, col, prefix=""):
    print(f"Loading SNPs from {path}...")
    with open(path, "rb") as f:
        data = zstd.ZstdDecompressor().decompress(f.read())
    df = pd.read_table(io.BytesIO(data), low_memory=False)
    return {f"{prefix}{s}" for s in df[col].dropna().unique()}


print("Identifying unique SNPs across databases...")
snps = load_snps(config["input_files"]["gwascatalog"], "SNP_ID_CURRENT", prefix="rs")
snps |= load_snps(config["input_files"]["disgenet"], "snpId")
print(f"Total unique SNPs identified: {len(snps):,}")


def update_cache(build, host, output_path, snp_list):
    print(f"\n--- Updating {build} cache ---")
    print(f"Target file: {output_path}")
    print(f"Ensembl Host: {host}")

    target_snps = list(snp_list)
    print(f"Fetching {len(target_snps):,} SNPs from Ensembl.")

    dataset = None
    for attempt in range(MAX_RETRIES):
        try:
            server = pybiomart.Server(host=host, use_cache=False)
            dataset = server.marts["ENSEMBL_MART_SNP"].datasets["hsapiens_snp"]
            break
        except Exception as e:
            if attempt < MAX_RETRIES - 1:
                sleep_time = BACKOFF_FACTOR**attempt
                print(
                    f"\nTransient error initializing BioMart server/dataset "
                    f"(attempt {attempt + 1}/{MAX_RETRIES}): {format_error(e)}. Retrying in {sleep_time}s..."
                )
                time.sleep(sleep_time)
            else:
                print(
                    f"\nPersistent error initializing BioMart server/dataset "
                    f"after {MAX_RETRIES} attempts: {format_error(e)}"
                )
                raise e

    new_results = []
    chunk_size = 500
    for i in tqdm(range(0, len(target_snps), chunk_size), desc="Fetching from BioMart"):
        chunk = target_snps[i : i + chunk_size]

        res = None
        for attempt in range(MAX_RETRIES):
            try:
                res = dataset.query(
                    attributes=[
                        "refsnp_id",
                        "chr_name",
                        "chrom_start",
                        "consequence_type_tv",
                        "ensembl_transcript_stable_id",
                    ],
                    filters={"snp_filter": chunk},
                    use_attr_names=True,
                )
                break
            except Exception as e:
                if attempt < MAX_RETRIES - 1:
                    sleep_time = BACKOFF_FACTOR**attempt
                    print(
                        f"\nTransient error querying chunk starting with {chunk[0]} "
                        f"(attempt {attempt + 1}/{MAX_RETRIES}): {format_error(e)}. Retrying in {sleep_time}s..."
                    )
                    time.sleep(sleep_time)
                else:
                    print(
                        f"\nPersistent error querying chunk starting with {chunk[0]} "
                        f"after {MAX_RETRIES} attempts: {format_error(e)}"
                    )
                    raise e

        if res is not None:
            new_results.append(res)

    if new_results:
        updated_df = pd.concat(new_results)
        updated_df.drop_duplicates(
            subset=["refsnp_id", "chr_name", "chrom_start", "consequence_type_tv"],
            inplace=True,
        )

        # Sanity check to ensure target SNPs have been annotated.
        annotated_snps = set(updated_df["refsnp_id"].dropna().unique())
        missing_snps = set(target_snps) - annotated_snps
        if missing_snps:
            print(
                f"\n[WARNING] {len(missing_snps):,} out of {len(target_snps):,} target SNPs "
                f"({len(missing_snps) / len(target_snps):.2%}) were not found or annotated in Ensembl {build}."
            )
            sample_missing = sorted(list(missing_snps))[:10]
            print(f"Sample of missing SNPs: {sample_missing}")
        else:
            print(
                f"\n[INFO] Sanity check passed: All {len(target_snps):,} SNPs were successfully annotated in Ensembl {build}!"
            )

        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        with open(output_path, "wb") as f:
            f.write(zstd.ZstdCompressor().compress(updated_df.to_csv().encode()))

        print(
            f"Successfully saved {len(updated_df)} total annotations to {output_path}"
        )
    else:
        print("No data retrieved.")


sources = [
    ("hg19", "http://grch37.ensembl.org"),
    ("hg38", "http://www.ensembl.org"),
]

for build, host in sources:
    update_cache(build, host, config["annotation_sources"][build], snps)
