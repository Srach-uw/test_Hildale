from __future__ import annotations

from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


def test_combined_confirmation_runner_freezes_the_intended_interaction() -> None:
    text = (REPO_ROOT / "cloud" / "ld_validation" / "run_paper_priors_reference_lcsc.sh").read_text(
        encoding="utf-8"
    )
    assert 'TARGET_CSV="targets_short_cadence_validation.csv"' in text
    assert 'CATALOG_SOURCE="sagear_ld_reference_catalog.csv"' in text
    assert 'RUN_ID="sagear_validation_paper_priors_reference_lcsc"' in text
    assert 'ALDERAAN_REPO="${ALDERAAN_PAPER_PRIOR_REPO:-$HOME/alderaan_sagear_paper_priors}"' in text
    assert 'CADENCE_MODE="both"' in text


def test_combined_confirmation_is_not_added_to_the_completed_matrix() -> None:
    matrix = (REPO_ROOT / "cloud" / "ld_validation" / "run_validation_matrix.sh").read_text(
        encoding="utf-8"
    )
    assert "paper_priors_reference_lcsc" not in matrix
