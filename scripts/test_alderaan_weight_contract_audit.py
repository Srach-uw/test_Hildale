from pathlib import Path

from alderaan_weight_contract_audit import audit_file


REPO_ROOT = Path(__file__).resolve().parents[1]


def test_existing_alderaan_fits_expose_nonuniform_dynesty_weights():
    fits_files = sorted(
        (REPO_ROOT / "data" / "alderaan_factorial_validation_20260715" / "results").rglob(
            "*-results.fits"
        )
    )
    assert fits_files
    row = audit_file(fits_files[0])
    assert row["status"] == "ok"
    assert row["weight_ess"] < row["rows"]
