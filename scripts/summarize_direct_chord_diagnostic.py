"""Summarize a bounded direct-chord nested-sampling diagnostic.

This program only describes the direct sampler output.  It does not compare
against the production eccentricity archive and does not make a population
claim.  A result with low effective posterior sample size is explicitly marked
as unsuitable for a parameterization comparison.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


def normalized_nested_weights(logwt: np.ndarray, logz_final: float) -> np.ndarray:
    """Recover normalized posterior weights from Dynesty's log-weight output."""
    raw = np.exp(np.asarray(logwt, dtype=float) - float(logz_final))
    if not np.all(np.isfinite(raw)) or raw.sum() <= 0:
        raise ValueError("nonfinite or zero nested posterior weights")
    return raw / raw.sum()


def weighted_quantiles(values: np.ndarray, weights: np.ndarray, quantiles: tuple[float, ...]) -> list[float]:
    """Return deterministic weighted quantiles for one posterior coordinate."""
    order = np.argsort(values)
    return [float(np.interp(q, np.cumsum(weights[order]), values[order])) for q in quantiles]


def summarize(path: Path) -> dict[str, object]:
    """Load the declared direct-run NPZ and record posterior quality metrics."""
    with np.load(path) as data:
        samples = np.asarray(data["samples"], dtype=float)
        weights = normalized_nested_weights(data["logwt"], float(data["logz"][-1]))
        logl = np.asarray(data["logl"], dtype=float)
    if samples.ndim != 2 or samples.shape[1] != 8 or len(samples) != len(weights):
        raise ValueError("unexpected direct diagnostic sample layout")
    ess = float(1.0 / np.sum(weights**2))
    return {
        "schema_version": 1,
        "scope": "Direct-chord diagnostic only; not a population result.",
        "input": str(path.resolve()),
        "nested_rows": int(len(weights)),
        "posterior_effective_sample_size": ess,
        "minimum_ess_for_route_comparison": 100.0,
        "suitable_for_route_comparison": bool(ess >= 100.0),
        "eccentricity_p16_p50_p84": weighted_quantiles(samples[:, 4], weights, (0.16, 0.50, 0.84)),
        "omega_radians_p16_p50_p84": weighted_quantiles(samples[:, 5], weights, (0.16, 0.50, 0.84)),
        "radius_ratio_p16_p50_p84": weighted_quantiles(samples[:, 2], weights, (0.16, 0.50, 0.84)),
        "impact_parameter_p16_p50_p84": weighted_quantiles(samples[:, 3], weights, (0.16, 0.50, 0.84)),
        "maximum_log_likelihood": float(np.max(logl)),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("samples", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = summarize(args.samples)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
