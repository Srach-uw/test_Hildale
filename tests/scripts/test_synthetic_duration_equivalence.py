from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

import synthetic_duration_equivalence as equivalence  # noqa: E402


@pytest.mark.parametrize(
    ("truth_e", "truth_omega", "truth_b"),
    [
        (0.0, 0.0, 0.40),
        (0.25, 1.1, 0.65),
    ],
)
def test_duration_reweighting_matches_direct_evaluation_on_synthetic_cases(
    truth_e: float, truth_omega: float, truth_b: float
):
    experiment = equivalence.DurationExperiment()
    truth = equivalence.synthetic_observation(
        experiment,
        eccentricity=truth_e,
        omega_radians=truth_omega,
        impact_parameter=truth_b,
        seed=71,
    )
    result = equivalence.run_comparison(experiment, truth["observed_duration_days"], sample_count=80_000, seed=101)

    assert result["support"]["duration_grid_normalized_mass"] == pytest.approx(1.0, abs=1e-10)
    assert result["support"]["prior_draw_fraction_within_duration_support"] > 0.99
    assert result["direct_evaluation"]["effective_sample_size"] > 1_000
    assert result["duration_posterior_importance_reweighting"]["effective_sample_size"] > 1_000
    assert result["equivalence"]["total_variation_distance"] < 2e-5
    assert "does not reproduce ALDERAAN" in result["model"]["limitation"]


def test_omitting_duration_prior_factor_fails_the_negative_control():
    experiment = equivalence.DurationExperiment(sigma_days=0.012)
    truth = equivalence.synthetic_observation(
        experiment,
        eccentricity=0.22,
        omega_radians=1.3,
        impact_parameter=0.55,
        seed=13,
    )
    result = equivalence.run_comparison(experiment, truth["observed_duration_days"], sample_count=100_000, seed=103)

    assert result["equivalence"]["total_variation_distance"] < 2e-5
    assert result["negative_control_omitted_factor"]["total_variation_distance"] > 0.01
    direct_median = result["direct_evaluation"]["eccentricity_p16_p50_p84"][1]
    broken_median = result["intentionally_omitted_duration_prior_factor"]["eccentricity_p16_p50_p84"][1]
    assert abs(direct_median - broken_median) > 0.003


def test_results_are_deterministic_and_include_support_diagnostics():
    experiment = equivalence.DurationExperiment()
    first = equivalence.run_comparison(experiment, 0.10, sample_count=20_000, seed=42)
    second = equivalence.run_comparison(experiment, 0.10, sample_count=20_000, seed=42)
    assert first == second
    assert set(first["support"]) == {
        "prior_draw_fraction_within_duration_support",
        "duration_grid_normalized_mass",
        "duration_posterior_edge_mass",
    }
    assert np.isfinite(first["support"]["duration_posterior_edge_mass"])


def test_cli_writes_self_describing_json(tmp_path: Path):
    output = tmp_path / "comparison.json"
    command = [
        sys.executable,
        str(ROOT / "scripts" / "synthetic_duration_equivalence.py"),
        "--output",
        str(output),
        "--truth-e",
        "0.2",
        "--truth-omega",
        "1.0",
        "--truth-b",
        "0.6",
        "--sample-count",
        "12000",
        "--seed",
        "4",
    ]
    completed = subprocess.run(command, check=True, capture_output=True, text=True)
    result = json.loads(output.read_text(encoding="utf-8"))
    assert "wrote" in completed.stdout
    assert result["synthetic_truth"]["truth_eccentricity"] == pytest.approx(0.2)
    assert result["equivalence"]["total_variation_distance"] < 2e-5


def test_invalid_support_is_rejected():
    with pytest.raises(equivalence.SyntheticDurationError, match="duration_max_days"):
        equivalence.DurationExperiment(duration_min_days=0.2, duration_max_days=0.1).validate()
