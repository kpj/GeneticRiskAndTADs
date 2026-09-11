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
