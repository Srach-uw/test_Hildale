"""Compare two equivalent toy photoeccentricity inference routes.

This module is a bounded synthetic check, not an ALDERAAN replacement.  It
uses a Gaussian likelihood for an exposure-integrated box-transit duration and
an intentionally small orbital model.  It does not implement ALDERAAN's light
curve, timing, limb-darkening, Gaussian-process, or nuisance-parameter
likelihood.  A successful comparison therefore checks only the stated change
of measure between a duration posterior and direct evaluation of this toy
model.

Both routes use the same prior: uniform eccentricity, argument of periastron,
and impact parameter, with a log-uniform interim prior on the apparent
duration.  The direct route evaluates the Gaussian duration likelihood at
prior draws.  The duration route evaluates the normalized interim duration
posterior and divides by its log-uniform prior.  Those weights are proportional
to the direct likelihood when the support and likelihood are the same.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
from typing import Any

import numpy as np


LIMITATION = (
    "This is a toy duration-likelihood equivalence check. It does not reproduce "
    "ALDERAAN's full transit, timing, limb-darkening, Gaussian-process, or "
    "nuisance-parameter likelihood."
)


class SyntheticDurationError(ValueError):
    """Raised when a synthetic duration experiment is not well specified."""


@dataclass(frozen=True)
class DurationExperiment:
    """Fully stated ingredients of the toy duration experiment."""

    base_duration_days: float = 0.12
    exposure_days: float = 0.020416666666666666
    sigma_days: float = 0.003
    duration_min_days: float = 0.015
    duration_max_days: float = 0.25
    eccentricity_max: float = 0.60
    impact_max: float = 0.90
    duration_grid_size: int = 32769

    def validate(self) -> None:
        positive = {
            "base_duration_days": self.base_duration_days,
            "exposure_days": self.exposure_days,
            "sigma_days": self.sigma_days,
            "duration_min_days": self.duration_min_days,
        }
        for name, value in positive.items():
            if not np.isfinite(value) or value <= 0:
                raise SyntheticDurationError(f"{name} must be finite and positive")
        if not np.isfinite(self.duration_max_days) or self.duration_max_days <= self.duration_min_days:
            raise SyntheticDurationError("duration_max_days must exceed duration_min_days")
        if not 0 < self.eccentricity_max < 1:
            raise SyntheticDurationError("eccentricity_max must lie in (0, 1)")
        if not 0 < self.impact_max < 1:
            raise SyntheticDurationError("impact_max must lie in (0, 1)")
        if self.duration_grid_size < 1025:
            raise SyntheticDurationError("duration_grid_size must be at least 1025")


def _normal_pdf(value: np.ndarray, mean: float, sigma: float) -> np.ndarray:
    standardized = (value - mean) / sigma
    return np.exp(-0.5 * standardized**2) / (sigma * np.sqrt(2.0 * np.pi))


def apparent_duration_days(
    eccentricity: np.ndarray | float,
    omega_radians: np.ndarray | float,
    impact_parameter: np.ndarray | float,
    experiment: DurationExperiment,
) -> np.ndarray:
    """Return the box-model duration after one exposure-time broadening.

    The intrinsic chord duration scales as ``sqrt(1-b^2) / g``, where
    ``g = (1 + e sin(omega)) / sqrt(1-e^2)``.  Adding one cadence is exact for
    the apparent width of a box transit convolved with a box exposure.  This
    is a declared synthetic observable, not a physical replacement for the
    integrated transit profile used by ALDERAAN.
    """
    experiment.validate()
    e = np.asarray(eccentricity, dtype=np.float64)
    omega = np.asarray(omega_radians, dtype=np.float64)
    impact = np.asarray(impact_parameter, dtype=np.float64)
    if np.any(~np.isfinite(e)) or np.any(~np.isfinite(omega)) or np.any(~np.isfinite(impact)):
        raise SyntheticDurationError("orbital inputs must be finite")
    if np.any(e < 0) or np.any(e > experiment.eccentricity_max):
        raise SyntheticDurationError("eccentricity falls outside the declared prior support")
    if np.any(impact < 0) or np.any(impact > experiment.impact_max):
        raise SyntheticDurationError("impact parameter falls outside the declared prior support")
    g = (1.0 + e * np.sin(omega)) / np.sqrt(1.0 - e**2)
    intrinsic = experiment.base_duration_days * np.sqrt(1.0 - impact**2) / g
    return intrinsic + experiment.exposure_days


def duration_prior_density(duration_days: np.ndarray | float, experiment: DurationExperiment) -> np.ndarray:
    """Evaluate the declared log-uniform interim duration density."""
    duration = np.asarray(duration_days, dtype=np.float64)
    result = np.zeros_like(duration, dtype=np.float64)
    support = (duration >= experiment.duration_min_days) & (duration <= experiment.duration_max_days)
    normalizer = np.log(experiment.duration_max_days / experiment.duration_min_days)
    result[support] = 1.0 / (duration[support] * normalizer)
    return result


def normalized_duration_posterior(
    observed_duration_days: float, experiment: DurationExperiment
) -> tuple[np.ndarray, np.ndarray, float]:
    """Construct the normalized interim duration posterior on a fixed grid."""
    if not np.isfinite(observed_duration_days):
        raise SyntheticDurationError("observed_duration_days must be finite")
    grid = np.linspace(experiment.duration_min_days, experiment.duration_max_days, experiment.duration_grid_size)
    unnormalized = _normal_pdf(grid, observed_duration_days, experiment.sigma_days) * duration_prior_density(grid, experiment)
    normalizer = float(np.trapezoid(unnormalized, grid))
    if not np.isfinite(normalizer) or normalizer <= 0:
        raise SyntheticDurationError("duration posterior has zero or nonfinite normalization")
    return grid, unnormalized / normalizer, normalizer


def _prior_draws(experiment: DurationExperiment, sample_count: int, seed: int) -> dict[str, np.ndarray]:
    if sample_count < 1000:
        raise SyntheticDurationError("sample_count must be at least 1000")
    rng = np.random.default_rng(seed)
    return {
        "eccentricity": rng.uniform(0.0, experiment.eccentricity_max, sample_count),
        "omega_radians": rng.uniform(0.0, 2.0 * np.pi, sample_count),
        "impact_parameter": rng.uniform(0.0, experiment.impact_max, sample_count),
    }


def _normalize_weights(raw_weights: np.ndarray) -> np.ndarray:
    weights = np.asarray(raw_weights, dtype=np.float64)
    if np.any(~np.isfinite(weights)) or np.any(weights < 0):
        raise SyntheticDurationError("weights must be finite and nonnegative")
    total = float(weights.sum())
    if not np.isfinite(total) or total <= 0:
        raise SyntheticDurationError("weights have zero total mass")
    return weights / total


def _weighted_quantile(values: np.ndarray, weights: np.ndarray, quantiles: tuple[float, ...]) -> list[float]:
    order = np.argsort(values)
    ordered_values = np.asarray(values)[order]
    ordered_weights = np.asarray(weights)[order]
    cdf = np.cumsum(ordered_weights) - 0.5 * ordered_weights
    cdf /= cdf[-1]
    return [float(np.interp(q, cdf, ordered_values)) for q in quantiles]


def _posterior_summary(draws: dict[str, np.ndarray], weights: np.ndarray) -> dict[str, Any]:
    e = draws["eccentricity"]
    impact = draws["impact_parameter"]
    duration = draws["apparent_duration_days"]
    return {
        "eccentricity_p16_p50_p84": _weighted_quantile(e, weights, (0.16, 0.50, 0.84)),
        "impact_parameter_p16_p50_p84": _weighted_quantile(impact, weights, (0.16, 0.50, 0.84)),
        "apparent_duration_days_p16_p50_p84": _weighted_quantile(duration, weights, (0.16, 0.50, 0.84)),
        "effective_sample_size": float(1.0 / np.sum(weights**2)),
    }


def _comparison_metrics(reference: np.ndarray, candidate: np.ndarray) -> dict[str, float]:
    return {
        "total_variation_distance": float(0.5 * np.sum(np.abs(reference - candidate))),
        "maximum_absolute_weight_difference": float(np.max(np.abs(reference - candidate))),
    }


def run_comparison(
    experiment: DurationExperiment,
    observed_duration_days: float,
    *,
    sample_count: int = 100_000,
    seed: int = 20260925,
) -> dict[str, Any]:
    """Compare direct and duration-posterior weights on identical prior draws.

    Reusing the same draws is intentional: it isolates the change-of-measure
    calculation from Monte Carlo variation.  Callers can run this function at
    different deterministic seeds to estimate sampling variation separately.
    """
    experiment.validate()
    draws = _prior_draws(experiment, sample_count, seed)
    duration = apparent_duration_days(
        draws["eccentricity"], draws["omega_radians"], draws["impact_parameter"], experiment
    )
    draws["apparent_duration_days"] = duration
    support = (duration >= experiment.duration_min_days) & (duration <= experiment.duration_max_days)
    likelihood = _normal_pdf(duration, observed_duration_days, experiment.sigma_days)
    direct_weights = _normalize_weights(likelihood * support)

    grid, duration_posterior, duration_normalizer = normalized_duration_posterior(observed_duration_days, experiment)
    posterior_at_duration = np.interp(duration, grid, duration_posterior, left=0.0, right=0.0)
    interim_prior_at_duration = duration_prior_density(duration, experiment)
    corrected_raw = np.divide(
        posterior_at_duration,
        interim_prior_at_duration,
        out=np.zeros_like(posterior_at_duration),
        where=interim_prior_at_duration > 0,
    )
    corrected_weights = _normalize_weights(corrected_raw)
    omitted_factor_weights = _normalize_weights(posterior_at_duration)

    edge_cells = max(1, experiment.duration_grid_size // 1000)
    grid_mass = float(np.trapezoid(duration_posterior, grid))
    return {
        "schema_version": 1,
        "model": {
            "prior": {
                "eccentricity": f"Uniform(0, {experiment.eccentricity_max})",
                "omega_radians": "Uniform(0, 2*pi)",
                "impact_parameter": f"Uniform(0, {experiment.impact_max})",
                "interim_apparent_duration_days": (
                    f"LogUniform({experiment.duration_min_days}, {experiment.duration_max_days})"
                ),
            },
            "likelihood": "Gaussian(observed apparent duration | box-exposure duration, sigma_days)",
            "exposure_model": "apparent duration = intrinsic box duration + one exposure",
            "limitation": LIMITATION,
        },
        "experiment": asdict(experiment),
        "observed_duration_days": float(observed_duration_days),
        "seed": int(seed),
        "sample_count": int(sample_count),
        "duration_posterior_normalizer": duration_normalizer,
        "support": {
            "prior_draw_fraction_within_duration_support": float(np.mean(support)),
            "duration_grid_normalized_mass": grid_mass,
            "duration_posterior_edge_mass": float(
                np.trapezoid(duration_posterior[:edge_cells], grid[:edge_cells])
                + np.trapezoid(duration_posterior[-edge_cells:], grid[-edge_cells:])
            ),
        },
        "direct_evaluation": _posterior_summary(draws, direct_weights),
        "duration_posterior_importance_reweighting": _posterior_summary(draws, corrected_weights),
        "intentionally_omitted_duration_prior_factor": _posterior_summary(draws, omitted_factor_weights),
        "equivalence": _comparison_metrics(direct_weights, corrected_weights),
        "negative_control_omitted_factor": _comparison_metrics(direct_weights, omitted_factor_weights),
    }


def synthetic_observation(
    experiment: DurationExperiment,
    *,
    eccentricity: float,
    omega_radians: float,
    impact_parameter: float,
    seed: int,
) -> dict[str, float]:
    """Generate one deterministic noisy duration observation from declared truth."""
    truth_duration = float(apparent_duration_days(eccentricity, omega_radians, impact_parameter, experiment))
    observed = float(truth_duration + np.random.default_rng(seed).normal(0.0, experiment.sigma_days))
    return {
        "truth_eccentricity": float(eccentricity),
        "truth_omega_radians": float(omega_radians),
        "truth_impact_parameter": float(impact_parameter),
        "truth_apparent_duration_days": truth_duration,
        "observed_duration_days": observed,
        "observation_seed": int(seed),
    }


def _json_default(value: Any) -> Any:
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(f"cannot serialize {type(value)!r}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True, help="JSON result path")
    parser.add_argument("--truth-e", type=float, default=0.0)
    parser.add_argument("--truth-omega", type=float, default=0.0)
    parser.add_argument("--truth-b", type=float, default=0.4)
    parser.add_argument("--observation-seed", type=int, default=23)
    parser.add_argument("--sample-count", type=int, default=100_000)
    parser.add_argument("--seed", type=int, default=20260925)
    parser.add_argument("--base-duration-days", type=float, default=DurationExperiment.base_duration_days)
    parser.add_argument("--exposure-days", type=float, default=DurationExperiment.exposure_days)
    parser.add_argument("--sigma-days", type=float, default=DurationExperiment.sigma_days)
    parser.add_argument("--duration-min-days", type=float, default=DurationExperiment.duration_min_days)
    parser.add_argument("--duration-max-days", type=float, default=DurationExperiment.duration_max_days)
    parser.add_argument("--eccentricity-max", type=float, default=DurationExperiment.eccentricity_max)
    parser.add_argument("--impact-max", type=float, default=DurationExperiment.impact_max)
    parser.add_argument("--duration-grid-size", type=int, default=DurationExperiment.duration_grid_size)
    args = parser.parse_args()

    experiment = DurationExperiment(
        base_duration_days=args.base_duration_days,
        exposure_days=args.exposure_days,
        sigma_days=args.sigma_days,
        duration_min_days=args.duration_min_days,
        duration_max_days=args.duration_max_days,
        eccentricity_max=args.eccentricity_max,
        impact_max=args.impact_max,
        duration_grid_size=args.duration_grid_size,
    )
    truth = synthetic_observation(
        experiment,
        eccentricity=args.truth_e,
        omega_radians=args.truth_omega,
        impact_parameter=args.truth_b,
        seed=args.observation_seed,
    )
    result = run_comparison(experiment, truth["observed_duration_days"], sample_count=args.sample_count, seed=args.seed)
    result["synthetic_truth"] = truth
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True, default=_json_default) + "\n", encoding="utf-8")
    print(f"wrote {args.output}")
    print(f"total variation (corrected): {result['equivalence']['total_variation_distance']:.3e}")
    print(f"total variation (omitted factor): {result['negative_control_omitted_factor']['total_variation_distance']:.3e}")


if __name__ == "__main__":
    main()
