"""Audit every defensible interpretation of Sagear's underspecified disk GMM.

This script never promotes a count-matching model to canonical status. It holds
the planet sample fixed, varies only documented classifier decisions, and emits
both the method table and a transparent prior sweep showing how strongly the
reported counts depend on an unstated component prior.
"""

from __future__ import annotations

import argparse
from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import logsumexp

from common import load_config, output_dir, root_path
from diagnose_sample import classify_with_chemical_gmm, disk_counts


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default=None)
    parser.add_argument("--sample", default=None)
    args = parser.parse_args()

    cfg = load_config(args.config)
    sample_path = Path(args.sample) if args.sample else output_dir() / "canonical_sample_strict.csv"
    sample = pd.read_csv(sample_path)
    apogee_path = root_path(cfg, "apogee_dr17_chemical")
    if apogee_path is None or not apogee_path.exists():
        raise FileNotFoundError(f"APOGEE calibration table not found: {apogee_path}")
    apogee = pd.read_csv(apogee_path)

    rows: list[dict[str, object]] = []
    sweep_rows: list[dict[str, object]] = []
    target = cfg["sagear_targets"]
    target_vector = np.array(
        [
            target["thin_singles_planets"],
            target["thick_singles_planets"],
            target["thin_multi_planets"],
            target["thick_multi_planets"],
        ],
        dtype=float,
    )

    for model, velocity, scope, flags, prior in product(
        ["conditioned", "pooled"],
        ["old_astropy", "direct"],
        ["planet_hosts", "all_kepler"],
        ["finite", "all_zero"],
        ["equal", "empirical"],
    ):
        # The component-prior switch has no effect for an unsupervised pooled GMM.
        if model == "pooled" and prior == "empirical":
            continue
        try:
            classified, note = classify_with_chemical_gmm(
                sample,
                apogee,
                cfg,
                model=model,
                velocity=velocity,
                scope=scope,
                flag_policy=flags,
                component_prior=prior,
            )
        except Exception as exc:
            rows.append(
                {
                    "model": model,
                    "velocity": velocity,
                    "scope": scope,
                    "flag_policy": flags,
                    "component_prior": prior,
                    "status": f"error:{type(exc).__name__}:{exc}",
                }
            )
            continue
        counts = compact_counts(classified)
        vector = np.array(
            [counts["thin_singles"], counts["thick_singles"], counts["thin_multis"], counts["thick_multis"]]
        )
        row = {
            "model": model,
            "velocity": velocity,
            "scope": scope,
            "flag_policy": flags,
            "component_prior": prior,
            "status": "ok",
            **counts,
            "l1_planet_count_delta": int(np.abs(vector - target_vector).sum()),
            "method_note": note,
        }
        rows.append(row)

        # The sweep replaces the component prior explicitly. Run it once per
        # likelihood construction instead of duplicating it for both priors.
        if model != "conditioned" or prior != "equal":
            continue
        ll_low = classified["kinematic_loglike_thin"].to_numpy(float)
        ll_high = classified["kinematic_loglike_thick"].to_numpy(float)
        finite = np.isfinite(ll_low) & np.isfinite(ll_high)
        for p_high in np.linspace(0.02, 0.50, 97):
            logp = np.column_stack(
                [ll_low[finite] + np.log1p(-p_high), ll_high[finite] + np.log(p_high)]
            )
            p = np.exp(logp[:, 1] - logsumexp(logp, axis=1))
            swept = classified.copy()
            swept.loc[finite, "P_thick"] = p
            swept.loc[finite, "disk"] = np.where(p > cfg["cuts"]["p_thick_threshold"], "thick", "thin")
            counts_s = compact_counts(swept)
            vector_s = np.array(
                [counts_s["thin_singles"], counts_s["thick_singles"], counts_s["thin_multis"], counts_s["thick_multis"]]
            )
            sweep_rows.append(
                {
                    "velocity": velocity,
                    "scope": scope,
                    "flag_policy": flags,
                    "component_prior_high": float(p_high),
                    **counts_s,
                    "l1_planet_count_delta": int(np.abs(vector_s - target_vector).sum()),
                }
            )

    methods = pd.DataFrame(rows).sort_values("l1_planet_count_delta")
    sweep = pd.DataFrame(sweep_rows).sort_values("l1_planet_count_delta")
    methods_path = output_dir() / "classifier_reconstruction_methods.csv"
    sweep_path = output_dir() / "classifier_reconstruction_prior_sweep.csv"
    methods.to_csv(methods_path, index=False)
    sweep.to_csv(sweep_path, index=False)
    print(methods.head(20).to_string(index=False))
    print("\nBest prior-sweep rows (diagnostic, never a selection rule):")
    print(sweep.head(20).to_string(index=False))
    print(f"\nWrote: {methods_path}")
    print(f"Wrote: {sweep_path}")


def compact_counts(df: pd.DataFrame) -> dict[str, int]:
    tab = disk_counts(df)
    lookup = {(row.disk, row.system): int(row.planets) for row in tab.itertuples()}
    return {
        "thin_singles": lookup.get(("thin", "single"), 0),
        "thick_singles": lookup.get(("thick", "single"), 0),
        "thin_multis": lookup.get(("thin", "multi"), 0),
        "thick_multis": lookup.get(("thick", "multi"), 0),
    }


if __name__ == "__main__":
    main()
