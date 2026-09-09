"""Diagnostic fits for testing Sagear Table 3 population-row ordering.

This is deliberately separate from the canonical Rayleigh fitter.  It uses the
same posterior-grid contract and transit-selection convention, but reports
MAP parameter values and model means only.  It is intended to test whether a
single downstream array-order permutation can explain all Table 3 model
families.  It is not a claim of exact reproduction of the paper's NumPyro run.
"""
from __future__ import annotations

import argparse
import itertools
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import gammaln

from hierarchical_rayleigh import (
    POPULATIONS,
    apply_grid_support_exclusion,
    load_population_masses,
    recompute_grid_support_flags,
    validate_summary_contract,
)
from common import trapezoid


PUBLISHED_BY_MODEL = {
    "beta": {
        "thick_singles": 0.058,
        "thick_multis": 0.042,
        "thin_singles": 0.023,
        "thin_multis": 0.030,
    },
    "monotonic_beta": {
        "thick_singles": 0.041,
        "thick_multis": 0.025,
        "thin_singles": 0.022,
        "thin_multis": 0.026,
    },
    "half_gaussian": {
        "thick_singles": 0.066,
        "thick_multis": 0.037,
        "thin_singles": 0.025,
        "thin_multis": 0.033,
    },
}
PUBLISHED_ORDER = ("thick_singles", "thick_multis", "thin_singles", "thin_multis")


def normalized_density(e: np.ndarray, raw: np.ndarray) -> np.ndarray:
    raw = np.clip(np.asarray(raw, dtype=float), 0.0, np.inf)
    return raw / max(float(trapezoid(raw, e)), 1e-300)


def beta_density(e: np.ndarray, alpha: float, beta: float, outlier_floor: float = 0.0) -> np.ndarray:
    logf = (alpha - 1.0) * np.log(np.clip(e, 1e-12, None)) + (beta - 1.0) * np.log1p(-e)
    logf -= gammaln(alpha) + gammaln(beta) - gammaln(alpha + beta)
    return normalized_density(e, np.exp(np.clip(logf, -700, 700)) + outlier_floor)


def half_gaussian_density(e: np.ndarray, sigma: float, outlier_floor: float = 0.0) -> np.ndarray:
    return normalized_density(e, np.exp(-0.5 * (e / sigma) ** 2) + outlier_floor)


def model_density(
    model: str, e: np.ndarray, x: np.ndarray, outlier_floor: float = 0.0
) -> tuple[np.ndarray, float, dict]:
    if model == "beta":
        alpha, beta = np.exp(x)
        f = beta_density(e, alpha, beta, outlier_floor)
        return f, float(trapezoid(f * e, e)), {"alpha": alpha, "beta": beta}
    if model == "monotonic_beta":
        alpha, beta = np.exp(x)
        f = beta_density(e, alpha, beta, outlier_floor)
        return f, float(trapezoid(f * e, e)), {"alpha": alpha, "beta": beta}
    if model == "half_gaussian":
        sigma = float(np.exp(x[0]))
        f = half_gaussian_density(e, sigma, outlier_floor)
        return f, float(trapezoid(f * e, e)), {"sigma": sigma}
    raise ValueError(model)


def fit_model(
    masses: np.ndarray,
    e: np.ndarray,
    model: str,
    starts: list[np.ndarray],
    outlier_floor: float = 0.0,
) -> dict:
    # Each row is a planet's posterior mass over e after the same omega
    # selection correction used by the Rayleigh implementation.
    def objective(x: np.ndarray) -> float:
        if model == "beta":
            alpha, beta = np.exp(x)
            # Sagear samples tau=1/sqrt(alpha+beta) on (0,1].
            if alpha + beta < 1.0:
                return 1e100
        f, _, _ = model_density(model, e, x, outlier_floor)
        terms = np.clip(masses @ f, 1e-300, None)
        return float(-np.log(terms).sum())

    bounds = {
        "beta": [(-8.0, 8.0), (-8.0, 8.0)],
        # Sagear's monotonic Beta is restricted to alpha < 1 and beta > 1.
        "monotonic_beta": [(-8.0, -1e-8), (1e-8, 8.0)],
        "half_gaussian": [(-8.0, 0.0)],
    }[model]
    best = None
    for start in starts:
        result = minimize(objective, start, method="L-BFGS-B", bounds=bounds)
        if best is None or result.fun < best.fun:
            best = result
    assert best is not None
    _, mean_e, pars = model_density(model, e, best.x, outlier_floor)
    return {"status": "ok" if best.success else "optimizer_warning", "nll": best.fun, "mean_e": mean_e, **pars}


def starts_for(model: str) -> list[np.ndarray]:
    if model == "half_gaussian":
        return [np.array([np.log(s)]) for s in (0.02, 0.05, 0.1, 0.3)]
    return [np.log(np.asarray(v, dtype=float)) for v in ((0.5, 4), (1, 10), (2, 20), (0.3, 2))]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--summary", required=True)
    ap.add_argument("--tag", required=True)
    ap.add_argument("--selection-mode", choices=["none", "manuscript_reciprocal"], default="manuscript_reciprocal")
    ap.add_argument("--allow-non-dynesty-weights", action="store_true")
    ap.add_argument("--allow-nonpaired-impact", action="store_true")
    ap.add_argument("--allow-mixed-posterior-sources", action="store_true")
    ap.add_argument("--allow-missing-qc-manifest", action="store_true")
    ap.add_argument("--no-qc", action="store_true")
    ap.add_argument(
        "--outlier-floor",
        type=float,
        default=0.0,
        help="Add a constant to each analytic density before normalization (Gilbert diagnostic).",
    )
    args = ap.parse_args()
    if not np.isfinite(args.outlier_floor) or not 0.0 <= args.outlier_floor <= 1.0:
        ap.error("--outlier-floor must be finite and between 0 and 1")
    if args.outlier_floor > 0 and "FLOOR" not in args.tag.upper():
        ap.error("A nonzero --outlier-floor requires FLOOR in --tag")
    summary = pd.read_csv(args.summary)
    validate_summary_contract(
        summary,
        allow_mixed_sources=args.allow_mixed_posterior_sources,
        allow_nonpaired_impact=args.allow_nonpaired_impact,
        allow_non_dynesty_weights=args.allow_non_dynesty_weights,
        allow_missing_qc=args.allow_missing_qc_manifest,
    )
    summary = recompute_grid_support_flags(summary)
    summary = apply_grid_support_exclusion(summary)
    if not args.no_qc and "qc_primary_exclude" in summary:
        summary = summary[~summary.qc_primary_exclude.fillna(True).astype(bool)].copy()

    rows = []
    for disk, system, label in POPULATIONS:
        sub = summary[(summary.disk == disk) & (summary.system == system)].reset_index(drop=True)
        masses, e = load_population_masses(sub, args.selection_mode != "none", args.selection_mode)
        for model in ("beta", "monotonic_beta", "half_gaussian"):
            fit = fit_model(masses, e, model, starts_for(model), args.outlier_floor)
            fit["comparison_value"] = fit["sigma"] if model == "half_gaussian" else fit["mean_e"]
            rows.append({"population": label, "n": len(sub), "model": model, "selection_mode": args.selection_mode, "outlier_floor": args.outlier_floor, **fit})

    result = pd.DataFrame(rows)
    out = Path("outputs") / f"table3_model_order_diagnostic_{args.tag}.csv"
    result.to_csv(out, index=False)

    # Enumerate all assignments of the four locally fitted population rows to
    # the four published rows.  This is a diagnostic of bookkeeping only.
    # Sagear reports sigma_HG for the half-Gaussian, but expected eccentricity
    # for the Beta-family models.  Compare like with like.
    piv = result.pivot(index="model", columns="population", values="comparison_value")
    perm_rows = []
    labels = list(PUBLISHED_ORDER)
    for model in piv.index:
        vals = piv.loc[model]
        published = PUBLISHED_BY_MODEL[model]
        for source_order in itertools.permutations(labels):
            assigned = [float(vals[s]) for s in source_order]
            target = [published[k] for k in labels]
            rmse = float(np.sqrt(np.mean((np.asarray(assigned) - target) ** 2)))
            perm_rows.append(
                {
                    "model": model,
                    "source_order": "|".join(source_order),
                    "target_order": "|".join(labels),
                    "rmse": rmse,
                    "mean_abs_error": float(
                        np.mean(np.abs(np.asarray(assigned) - target))
                    ),
                }
            )
    perms = pd.DataFrame(perm_rows).sort_values(["model", "rmse"])
    p_out = Path("outputs") / f"table3_model_order_permutations_{args.tag}.csv"
    perms.to_csv(p_out, index=False)
    md = Path("outputs") / f"table3_model_order_diagnostic_{args.tag}.md"
    lines = [f"# Table 3 model-order diagnostic: {args.tag}", "", "This is a bookkeeping diagnostic, not a physical disk relabeling or exact NumPyro reproduction.", "", "## Local fits", "", result.to_markdown(index=False), "", "## Best assignment per model", ""]
    for model in piv.index:
        lines.append(perms[perms.model == model].head(1).to_markdown(index=False))
        lines.append("")
    md.write_text("\n".join(lines), encoding="utf-8")
    print(result.to_string(index=False))
    print(f"WROTE {out}")
    print(f"WROTE {md}")


if __name__ == "__main__":
    main()
