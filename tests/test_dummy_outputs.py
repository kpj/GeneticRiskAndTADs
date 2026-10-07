from pathlib import Path

import pandas as pd
import pytest

from tests.schemas import FinalEnrichmentSchema, InitialDatabaseSchema

RESULTS_DIR = Path("results/test_dummy")


@pytest.fixture(scope="session")
def check_results_exist():
    """Ensure pipeline output directory exists before running assertions."""
    if not RESULTS_DIR.exists():
        pytest.skip(
            f"Pipeline output directory '{RESULTS_DIR}' does not exist. Run tests/run.sh first."
        )


def test_required_deliverables_exist(check_results_exist):
    """Ensure all core target files and figure directories are generated and non-empty."""
    expected_files = [
        RESULTS_DIR / "databases" / "initial.csv",
        RESULTS_DIR / "results" / "final_enr.csv.gz",
        RESULTS_DIR / "results" / "database_statistics" / "gene_counts.pdf",
        RESULTS_DIR
        / "results"
        / "database_statistics"
        / "disease_count_distribution.pdf",
        RESULTS_DIR / "publication_figures" / "main" / "figure1.pdf",
        RESULTS_DIR / "publication_figures" / "main" / "figure2.pdf",
        RESULTS_DIR / "publication_figures" / "main" / "figure3.pdf",
        RESULTS_DIR / "publication_figures" / "main" / "figure4.pdf",
    ]
    for path in expected_files:
        assert path.exists(), f"Expected deliverable not found: {path}"
        assert path.stat().st_size > 0, f"Expected deliverable is empty: {path}"


def test_initial_database_schema_and_invariants(check_results_exist):
    """Validate initial database table schema, types, and domain invariants."""
    db_path = RESULTS_DIR / "databases" / "initial.csv"
    assert db_path.exists(), f"Initial database missing: {db_path}"

    df = pd.read_csv(db_path)
    # Pandera lazy validation collects and reports all schema/check violations
    InitialDatabaseSchema.validate(df, lazy=True)

    # Explicit check for the GWAS date cutoff filter
    assert "rs999" not in df["snpId"].values, (
        "rs999 should have been filtered out by gwas_date_range cutoff!"
    )


def test_final_enrichments_schema_and_invariants(check_results_exist):
    """Validate final enrichment table schema, probability ranges, and relational constraints."""
    enr_path = RESULTS_DIR / "results" / "final_enr.csv.gz"
    assert enr_path.exists(), f"Final enrichments file missing: {enr_path}"

    df = pd.read_csv(enr_path)
    FinalEnrichmentSchema.validate(df, lazy=True)


def test_all_configured_callers_and_majority_votes_present(check_results_exist):
    """Verify presence of multiple callers (TopDom + cooltools) and semantic majority votes."""
    enr_path = RESULTS_DIR / "results" / "final_enr.csv.gz"
    df = pd.read_csv(enr_path)

    configs_present = set(df["caller_config"].unique())
    expected_configs = {
        "majority_vote",
        "majority_vote_narrow",
        "topdom_majority_vote",
        "cooltools_majority_vote",
        "topdom_w9",
        "topdom_w10",
        "topdom_w11",
        "cooltools_w3mb",
    }
    missing = expected_configs - configs_present
    assert not missing, f"Expected caller configs missing from final enrichments: {missing}"

    # Verify per-caller TAD files exist and have non-empty domain predictions
    for cfg in ["topdom_w9", "topdom_w10", "topdom_w11", "cooltools_w3mb"]:
        tad_file = RESULTS_DIR / "tads" / "data" / f"tads.dummy_name.{cfg}.csv"
        assert tad_file.exists(), f"TAD file missing for caller {cfg}: {tad_file}"
        df_tad = pd.read_csv(tad_file)
        assert len(df_tad) > 0, f"No domains called by {cfg} in {tad_file}"
        assert {"chrname", "tad_start", "tad_stop"}.issubset(df_tad.columns)


def test_symmetric_boundaries_computed(check_results_exist):
    """Verify that both asymmetric (*in) and symmetric (*sym) boundary geometries are computed."""
    enr_path = RESULTS_DIR / "results" / "final_enr.csv.gz"
    df = pd.read_csv(enr_path)

    tad_types = set(df["TAD_type"].unique())
    assert {"5in", "20in"}.issubset(tad_types), f"Missing asymmetric boundaries in {tad_types}"
    assert {"10sym", "20sym"}.issubset(tad_types), f"Missing symmetric boundaries in {tad_types}"


def test_no_legacy_artifacts_exist(check_results_exist):
    """Verify that legacy numeric tokens (.0. / .1.) do not exist in generated outputs."""
    legacy_files = list(RESULTS_DIR.glob("**/*.*.0.*")) + list(
        RESULTS_DIR.glob("**/*.*.1.*")
    )
    assert not legacy_files, f"Legacy numeric token files found: {legacy_files}"


def test_publication_figures_are_valid_pdfs(check_results_exist):
    """Verify that all generated figures exist and have valid PDF headers."""
    pdf_files = list((RESULTS_DIR / "publication_figures").glob("**/*.pdf")) + list(
        (RESULTS_DIR / "results" / "database_statistics").glob("*.pdf")
    )

    assert len(pdf_files) > 0, "No PDF figures found to validate"

    for pdf_path in pdf_files:
        assert pdf_path.stat().st_size > 500, (
            f"PDF file suspiciously small (<500B): {pdf_path}"
        )
        with open(pdf_path, "rb") as f:
            header = f.read(4)
            assert header == b"%PDF", (
                f"File does not have valid PDF magic bytes: {pdf_path}"
            )
