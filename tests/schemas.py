import pandera.pandas as pa
from pandera.pandas import Check, Column

InitialDatabaseSchema = pa.DataFrameSchema(
    columns={
        "diseaseId": Column(pa.String, Check.str_matches(r"^(EFO|DOID)_\d+$")),
        "snpId": Column(pa.String, Check.str_startswith("rs")),
        "snp_source": Column(pa.String),
        "diseaseIdType": Column(pa.String),
        "odds_ratio": Column(pa.Float64, Check.ge(0.0), nullable=True),
        "diseaseLabel": Column(pa.String),
        "is_cancer": Column(pa.Bool),
        "chromosome_hg19": Column(pa.Int64),
        "chromosome_hg38": Column(pa.Int64),
        "position_hg19": Column(pa.Int64, Check.gt(0)),
        "position_hg38": Column(pa.Int64, Check.gt(0)),
        "variant_group_hg19": Column(pa.String),
        "variant_group_hg38": Column(pa.String),
        "variant_type_hg19": Column(pa.String),
        "variant_type_hg38": Column(pa.String),
    },
    checks=[
        Check(
            lambda df: len(df) > 0,
            name="initial_db_not_empty",
            error="Initial database should not be empty",
        ),
        Check(
            lambda df: "rs999" not in df["snpId"].values,
            name="rs999_filtered_by_gwas_date_range",
            error="rs999 should have been filtered out by gwas_date_range!",
        ),
    ],
    strict=False,  # Allow additional filter_* columns
)

FinalEnrichmentSchema = pa.DataFrameSchema(
    columns={
        "diseaseId": Column(pa.String),
        "#snp": Column(pa.Int64, Check.ge(0)),
        "#border_snp": Column(pa.Int64, Check.ge(0)),
        "pval_tad": Column(pa.Float64, Check.in_range(0.0, 1.0), nullable=True),
        "pval_border": Column(pa.Float64, Check.in_range(0.0, 1.0), nullable=True),
        "pval_outside": Column(pa.Float64, Check.in_range(0.0, 1.0), nullable=True),
        "pval_tad__notcorrected": Column(
            pa.Float64, Check.in_range(0.0, 1.0), nullable=True
        ),
        "pval_border__notcorrected": Column(
            pa.Float64, Check.in_range(0.0, 1.0), nullable=True
        ),
        "pval_outside__notcorrected": Column(
            pa.Float64, Check.in_range(0.0, 1.0), nullable=True
        ),
        "TAD_type": Column(pa.String),
        "is_cancer": Column(pa.Bool),
        "tad_source": Column(pa.String),
        "window_size": Column(pa.Int64, Check.ge(0)),
        "filter": Column(pa.String),
    },
    checks=[
        Check(
            lambda df: len(df) > 0,
            name="final_enr_not_empty",
            error="Final enrichment dataframe should not be empty",
        ),
        Check(
            lambda df: (df["#snp"] >= df["#border_snp"]).all(),
            name="border_snp_lte_total_snp",
            error="Number of border SNPs (#border_snp) cannot exceed total SNPs (#snp)",
        ),
        Check(
            lambda df: (df["tad_source"] == "dummy_name").all(),
            name="tad_source_matches_dummy",
            error="TAD source must match configured sample 'dummy_name'",
        ),
    ],
    strict=False,
)
