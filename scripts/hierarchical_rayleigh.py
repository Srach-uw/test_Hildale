from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import output_dir, trapezoid


POPULATIONS = [
    ("thick", "single", "thick_singles"),
    ("thin", "single", "thin_singles"),
    ("thick", "multi", "thick_multis"),
    ("thin", "multi", "thin_multis"),
]
SELECTION_MODES = [
    "legacy_forward_norm",
    "manuscript_reciprocal",
    "manuscript_reciprocal_with_norm",
    "none",
]

GRID_SUPPORT_FALLBACK_MIN = float(np.sqrt(1.0 - 0.95**2) / (1.0 + 0.95))
GRID_SUPPORT_FALLBACK_MAX = float(np.sqrt(1.0 - 0.95**2) / (1.0 - 0.95))


def main() -> None:
    parser = argparse.ArgumentParser(description="Fit Sagear-comparable Rayleigh population models to e,omega posterior grids.")
    parser.add_argument("--summary", default=None)
    parser.add_argument("--sigma-min", type=float, default=1e-4)
    parser.add_argument("--sigma-max", type=float, default=1.0)
    parser.add_argument("--n-sigma", type=int, default=2000)
    parser.add_argument("--log-sigma-grid", action="store_true", help="Diagnostic only; posterior still uses uniform sigma prior with cell widths.")
    parser.add_argument("--ignore-transit-selection", action="store_true", help="Diagnostic only: omit transit-probability weighting.")
    parser.add_argument(
        "--selection-mode",
        choices=SELECTION_MODES,
        default="manuscript_reciprocal",
        help=(
            "Transit-selection convention. legacy_forward_norm reproduces the older pipeline; "
            "manuscript_reciprocal follows the commented Sagear HBM equation using 1/p_transit; "
            "none omits transit selection."
        ),
    )
    parser.add_argument("--exclude-zeta-outside-grid-support", action="store_true", help="Drop rows whose zeta median/p16/p84 is outside e<=0.95 support.")
    parser.add_argument(
        "--exclude-qc-primary",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Drop rows marked qc_primary_exclude (canonical default); use --no-exclude-qc-primary only for sensitivity runs.",
    )
    parser.add_argument(
        "--allow-mixed-posterior-sources",
        action="store_true",
        help="Diagnostic only: allow archive and paired-ALDERAAN posterior products in one population fit.",
    )
    parser.add_argument(
        "--allow-nonpaired-impact",
        action="store_true",
        help="Diagnostic only: allow posterior products that did not preserve paired ALDERAAN impact samples.",
    )
    parser.add_argument(
        "--allow-missing-qc-manifest",
        action="store_true",
        help="Diagnostic only: allow a summary without deterministic qc_primary_exclude/qc_reasons fields.",
    )
    parser.add_argument("--exclude-kois", default=None, help="Comma-separated kepoi_name values to drop before fitting.")
    parser.add_argument("--out-tag", default=None, help="Suffix for output CSV names.")
    parser.add_argument("--diagnostics", action="store_true", help="Write leverage, leave-10%-out, and ECDF diagnostics.")
    parser.add_argument("--leaveout-trials", type=int, default=10)
    parser.add_argument("--random-seed", type=int, default=20260707)
    parser.add_argument("--min-berger-logg", type=float, default=None, help="Labeled dwarf-sample sensitivity cut.")
    parser.add_argument("--max-berger-radius", type=float, default=None, help="Labeled dwarf-sample sensitivity cut in solar radii.")
    parser.add_argument("--berger2018-evol", type=int, default=None, help="Labeled Berger+2018 evolutionary-state sensitivity.")
    parser.add_argument("--confirmed-only", action="store_true", help="Diagnostic KOI-disposition sensitivity.")
    parser.add_argument("--min-koi-snr", type=float, default=None, help="Diagnostic transit-S/N sensitivity.")
    parser.add_argument("--max-planet-radius", type=float, default=None, help="Diagnostic radius sensitivity in Earth radii.")
    args = parser.parse_args()

    if args.summary is None:
        parser.error("--summary is required so a canonical run cannot silently fit a stale posterior product")
    summary_path = Path(args.summary)
    if not summary_path.exists():
        raise FileNotFoundError(f"Posterior summary not found: {summary_path}")
    summary = pd.read_csv(summary_path)
    validate_summary_contract(
        summary,
        allow_mixed_sources=args.allow_mixed_posterior_sources,
        allow_nonpaired_impact=args.allow_nonpaired_impact,
        allow_missing_qc=args.allow_missing_qc_manifest,
    )
    summary = recompute_grid_support_flags(summary)
    for option, column, operator in [
        (args.min_berger_logg, "berger_logg", "min"),
        (args.max_berger_radius, "berger_rad", "max"),
        (args.berger2018_evol, "berger2018_evol", "equal"),
    ]:
        if option is None:
            continue
        if column not in summary.columns:
            raise ValueError(f"Requested stellar sensitivity cut requires summary column: {column}")
        values = pd.to_numeric(summary[column], errors="coerce")
        keep = values >= option if operator == "min" else (values <= option if operator == "max" else values == option)
        before = len(summary)
        summary = summary[keep].reset_index(drop=True)
        print(f"stellar sensitivity {column} {operator} {option}: dropped {before - len(summary)} rows")
    for option, column, label in [
        ("CONFIRMED" if args.confirmed_only else None, "koi_disposition", "confirmed-only"),
        (args.min_koi_snr, "koi_model_snr", "minimum KOI S/N"),
        (args.max_planet_radius, "koi_prad", "maximum planet radius"),
    ]:
        if option is None:
            continue
        if column not in summary.columns:
            raise ValueError(f"Requested {label} sensitivity requires summary column: {column}")
        before = len(summary)
        if column == "koi_disposition":
            keep = summary[column].astype(str).str.upper().eq(str(option))
        elif column == "koi_model_snr":
            keep = pd.to_numeric(summary[column], errors="coerce") >= float(option)
        else:
            keep = pd.to_numeric(summary[column], errors="coerce") < float(option)
        summary = summary[keep].reset_index(drop=True)
        print(f"diagnostic sensitivity {label} {option}: dropped {before - len(summary)} rows")
    if args.exclude_zeta_outside_grid_support:
        summary = apply_grid_support_exclusion(summary)
    if args.exclude_qc_primary:
        if "qc_primary_exclude" not in summary.columns:
            raise ValueError("--exclude-qc-primary requested, but summary lacks qc_primary_exclude")
        before = len(summary)
        summary = summary[~summary["qc_primary_exclude"].fillna(True).astype(bool)].reset_index(drop=True)
        print(f"--exclude-qc-primary dropped {before - len(summary)} rows")
    if args.exclude_kois:
        drop = {k.strip() for k in args.exclude_kois.split(",") if k.strip()}
        before = len(summary)
        summary = summary[~summary["kepoi_name"].astype(str).isin(drop)].reset_index(drop=True)
        print(f"--exclude-kois dropped {before - len(summary)} rows (requested {len(drop)})")

    sigmas = (np.geomspace if args.log_sigma_grid else np.linspace)(args.sigma_min, args.sigma_max, args.n_sigma)
    selection_mode = "none" if args.ignore_transit_selection else args.selection_mode
    apply_selection = selection_mode != "none"
    rows: list[dict[str, object]] = []
    leverage_rows: list[pd.DataFrame] = []
    leaveout_rows: list[dict[str, object]] = []
    topk_rows: list[dict[str, object]] = []
    rng = np.random.default_rng(args.random_seed)

    for disk, system, label in POPULATIONS:
        sub = summary[(summary["disk"] == disk) & (summary["system"] == system)].reset_index(drop=True)
        if len(sub) < 5:
            rows.append({"population": label, "n": len(sub), "status": "skipped_n_lt_5"})
            continue
        masses, e_grid = load_population_masses(
            sub, apply_transit_selection=apply_selection, selection_mode=selection_mode
        )
        fit = fit_from_mass_matrix(
            masses, e_grid, sigmas, apply_transit_selection=apply_selection, selection_mode=selection_mode
        )
        rows.append({"population": label, "n": len(sub), "status": "ok", "selection_mode": selection_mode, **fit})
        if args.diagnostics:
            leverage = per_planet_leverage(sub, masses, e_grid, sigmas, fit, apply_selection, label, selection_mode)
            leverage_rows.append(leverage)
            leaveout_rows.extend(
                leaveout_diagnostics(
                    sub, masses, e_grid, sigmas, apply_selection, label, rng, args.leaveout_trials, selection_mode
                )
            )
            topk_rows.extend(
                topk_diagnostics(sub, masses, e_grid, sigmas, apply_selection, label, leverage, selection_mode)
            )

    out = pd.DataFrame(rows)
    if selection_mode == "none":
        suffix = "no_transit_selection"
    elif selection_mode == "legacy_forward_norm":
        suffix = "transit_selection"
    else:
        suffix = f"transit_selection_{selection_mode}"
    tag = f"_{args.out_tag}" if args.out_tag else ""
    out_path = output_dir() / f"rayleigh_population_fit_{suffix}{tag}.csv"
    if not args.out_tag:
        out.to_csv(output_dir() / "rayleigh_population_fit.csv", index=False)
    out.to_csv(out_path, index=False)
    print(out.to_string(index=False))
    print(f"\nWrote: {out_path}")

    if args.diagnostics:
        if leverage_rows:
            leverage = pd.concat(leverage_rows, ignore_index=True)
            leverage_path = output_dir() / f"rayleigh_per_planet_leverage{tag}.csv"
            leverage.to_csv(leverage_path, index=False)
            print(f"Wrote: {leverage_path}")
        if leaveout_rows:
            leaveout = pd.DataFrame(leaveout_rows)
            leaveout_path = output_dir() / f"rayleigh_leave10out{tag}.csv"
            leaveout.to_csv(leaveout_path, index=False)
            print(f"Wrote: {leaveout_path}")
        if topk_rows:
            topk = pd.DataFrame(topk_rows)
            topk_path = output_dir() / f"rayleigh_top_leverage_out{tag}.csv"
            topk.to_csv(topk_path, index=False)
            print(f"Wrote: {topk_path}")
        ecdf_path = output_dir() / f"eccentricity_zeta_ecdf{tag}.png"
        make_ecdf_plot(summary, ecdf_path)
        print(f"Wrote: {ecdf_path}")


def validate_summary_contract(
    summary: pd.DataFrame,
    *,
    allow_mixed_sources: bool = False,
    allow_nonpaired_impact: bool = False,
    allow_missing_qc: bool = False,
) -> None:
    required = {"kepoi_name", "disk", "system", "posterior_file"}
    missing = sorted(required - set(summary.columns))
    if missing:
        raise ValueError(f"Posterior summary is missing required fields: {missing}")

    if not allow_missing_qc:
        qc_required = {"qc_primary_exclude", "qc_reasons"}
        qc_missing = sorted(qc_required - set(summary.columns))
        if qc_missing:
            raise ValueError(
                "Canonical population fit requires a deterministic QC manifest; "
                f"missing fields: {qc_missing}. Use --allow-missing-qc-manifest only for diagnostics."
            )

    if "posterior_source" not in summary.columns:
        if not allow_mixed_sources:
            raise ValueError(
                "Canonical population fit requires posterior_source provenance. "
                "Use --allow-mixed-posterior-sources only for legacy diagnostics."
            )
    else:
        sources = sorted(summary["posterior_source"].dropna().astype(str).unique())
        if len(sources) != 1 and not allow_mixed_sources:
            raise ValueError(
                "Canonical population fit cannot mix posterior constructions; "
                f"found posterior_source={sources}. Fit each source separately or rerun all targets uniformly."
            )

    if "impact_mode" not in summary.columns:
        if not allow_nonpaired_impact:
            raise ValueError(
                "Canonical population fit requires impact_mode provenance. "
                "Use --allow-nonpaired-impact only for legacy diagnostics."
            )
    else:
        impact_modes = sorted(summary["impact_mode"].dropna().astype(str).str.lower().unique())
        if impact_modes != ["alderaan"] and not allow_nonpaired_impact:
            raise ValueError(
                "Canonical population fit requires paired ALDERAAN impact samples; "
                f"found impact_mode={impact_modes}."
            )


def recompute_grid_support_flags(summary: pd.DataFrame) -> pd.DataFrame:
    out = summary.copy()
    zeta_summary_columns = {"zeta_median", "zeta_p16", "zeta_p84"}
    has_zeta_columns = zeta_summary_columns.issubset(out.columns)
    has_any_zeta = has_zeta_columns and bool(out[list(zeta_summary_columns)].notna().to_numpy().any())
    if not has_any_zeta:
        # Direct sample-level importance posteriors do not require a zeta KDE
        # or its finite-grid support gate.
        for col in zeta_summary_columns:
            if col not in out.columns:
                out[col] = np.nan
        out["zeta_support_applicable"] = False
        out["zeta_median_outside_grid_support"] = False
        out["zeta_p16_outside_grid_support"] = False
        out["zeta_p84_outside_grid_support"] = False
        out["zeta_any_summary_outside_grid_support"] = False
        return out
    out["zeta_support_applicable"] = True
    for col in ["zeta_median", "zeta_p16", "zeta_p84", "zeta_grid_min", "zeta_grid_max"]:
        if col in out.columns:
            out[col] = pd.to_numeric(out[col], errors="coerce")
    gmin = out["zeta_grid_min"] if "zeta_grid_min" in out.columns else pd.Series(np.nan, index=out.index)
    gmax = out["zeta_grid_max"] if "zeta_grid_max" in out.columns else pd.Series(np.nan, index=out.index)
    out["zeta_grid_min"] = gmin.fillna(GRID_SUPPORT_FALLBACK_MIN)
    out["zeta_grid_max"] = gmax.fillna(GRID_SUPPORT_FALLBACK_MAX)
    out["zeta_median_outside_grid_support"] = ~out["zeta_median"].between(out["zeta_grid_min"], out["zeta_grid_max"])
    out["zeta_p16_outside_grid_support"] = ~out["zeta_p16"].between(out["zeta_grid_min"], out["zeta_grid_max"])
    out["zeta_p84_outside_grid_support"] = ~out["zeta_p84"].between(out["zeta_grid_min"], out["zeta_grid_max"])
    out["zeta_any_summary_outside_grid_support"] = (
        out["zeta_median_outside_grid_support"]
        | out["zeta_p16_outside_grid_support"]
        | out["zeta_p84_outside_grid_support"]
    )
    return out


def apply_grid_support_exclusion(summary: pd.DataFrame) -> pd.DataFrame:
    summary = recompute_grid_support_flags(summary)
    bad = summary["zeta_any_summary_outside_grid_support"].fillna(True).astype(bool)
    dropped = summary[bad]
    if len(dropped):
        print(f"--exclude-zeta-outside-grid-support dropped {len(dropped)} rows:")
        print(dropped.groupby(["disk", "system"]).size().to_string())
        print("dropped kepoi_names: " + ", ".join(sorted(dropped["kepoi_name"].astype(str))))
    return summary[~bad].reset_index(drop=True)


def selection_mode_from_legacy(
    apply_transit_selection: bool | None = None,
    selection_mode: str | None = None,
) -> str:
    if selection_mode is not None:
        return selection_mode
    return "legacy_forward_norm" if apply_transit_selection else "none"


def fit_rayleigh(
    files: list[str],
    sigmas: np.ndarray,
    apply_transit_selection: bool = True,
    selection_mode: str | None = None,
) -> dict:
    mode = selection_mode_from_legacy(apply_transit_selection, selection_mode)
    masses, e_grid = load_population_masses_from_files(files, apply_transit_selection, mode)
    return fit_from_mass_matrix(masses, e_grid, sigmas, apply_transit_selection, mode)


def load_population_masses(
    summary: pd.DataFrame,
    apply_transit_selection: bool,
    selection_mode: str | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    mode = selection_mode_from_legacy(apply_transit_selection, selection_mode)
    return load_population_masses_from_files(summary["posterior_file"].tolist(), apply_transit_selection, mode)


def load_population_masses_from_files(
    files: list[str],
    apply_transit_selection: bool,
    selection_mode: str | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    mode = selection_mode_from_legacy(apply_transit_selection, selection_mode)
    e_masses = []
    e_grid = None
    omega_grid = None
    for file in files:
        with np.load(file, allow_pickle=False) as data:
            required = {"e_grid", "omega_grid", "posterior", "include_transit_prior"}
            missing = sorted(required - set(data.files))
            if missing:
                raise ValueError(f"Posterior product is missing required arrays {missing}: {file}")
            include_transit_prior = bool(np.asarray(data["include_transit_prior"]).item())
            if apply_transit_selection and include_transit_prior is not False:
                raise ValueError(
                    f"Canonical transit-selection fit requires include_transit_prior=False, got {include_transit_prior}: {file}"
                )
            this_e = np.asarray(data["e_grid"], dtype=float)
            this_omega = np.asarray(data["omega_grid"], dtype=float)
            posterior = np.asarray(data["posterior"], dtype=float)

        if this_e.ndim != 1 or this_omega.ndim != 1:
            raise ValueError(f"e_grid and omega_grid must be one-dimensional: {file}")
        if len(this_e) < 2 or len(this_omega) < 2:
            raise ValueError(f"Posterior grids are too short: {file}")
        if not np.all(np.isfinite(this_e)) or not np.all(np.diff(this_e) > 0):
            raise ValueError(f"e_grid must be finite and strictly increasing: {file}")
        if not np.all(np.isfinite(this_omega)) or not np.all(np.diff(this_omega) > 0):
            raise ValueError(f"omega_grid must be finite and strictly increasing: {file}")
        if posterior.shape != (len(this_e), len(this_omega)):
            raise ValueError(
                f"Posterior shape {posterior.shape} does not match grids "
                f"({len(this_e)}, {len(this_omega)}): {file}"
            )
        if not np.all(np.isfinite(posterior)):
            raise ValueError(f"Posterior contains NaN or infinite values: {file}")
        if np.any(posterior < 0):
            raise ValueError(f"Posterior contains negative probability mass: {file}")
        total = float(posterior.sum())
        if not np.isfinite(total) or total <= 0:
            raise ValueError(f"Posterior has zero or invalid total mass: {file}")
        posterior = posterior / total

        if e_grid is None:
            e_grid = this_e
            omega_grid = this_omega
        elif not np.array_equal(e_grid, this_e) or not np.array_equal(omega_grid, this_omega):
            raise ValueError(
                "Canonical population fit requires identical e/omega grids across planets; "
                f"grid mismatch at {file}"
            )
        if mode == "legacy_forward_norm":
            weight = transit_probability_weight(this_e, this_omega)
            e_mass = np.sum(posterior * weight, axis=1)
        elif mode in {"manuscript_reciprocal", "manuscript_reciprocal_with_norm"}:
            weight = transit_probability_weight(this_e, this_omega)
            e_mass = np.sum(posterior / weight, axis=1)
        elif mode == "none":
            e_mass = np.sum(posterior, axis=1)
        else:
            raise ValueError(f"Unknown selection mode: {mode}")
        e_masses.append(e_mass)
    if e_grid is None:
        raise ValueError("No posterior files supplied")
    return np.vstack(e_masses), e_grid


def fit_from_mass_matrix(
    mass_matrix: np.ndarray,
    e_grid: np.ndarray,
    sigmas: np.ndarray,
    apply_transit_selection: bool = True,
    selection_mode: str | None = None,
) -> dict:
    mode = selection_mode_from_legacy(apply_transit_selection, selection_mode)
    rays, normalizers = rayleigh_grid(e_grid, sigmas, apply_transit_selection, mode)
    terms = mass_matrix @ rays
    if mode in {"legacy_forward_norm", "manuscript_reciprocal_with_norm"}:
        terms = terms / normalizers
    log_terms = np.log(np.clip(terms, 1e-300, None))
    lls = np.sum(log_terms, axis=0)
    best = int(np.argmax(lls))
    weights = posterior_weights_from_ll(sigmas, lls)
    sigma_q16, sigma_q50, sigma_q84 = weighted_quantile(sigmas, weights, [0.16, 0.5, 0.84])
    expected_grid = sigmas * np.sqrt(np.pi / 2.0)
    rays, _ = rayleigh_grid(e_grid, sigmas, apply_transit_selection, mode)
    truncated_expected_grid = trapezoid(rays * e_grid[:, None], e_grid, axis=0)
    e_q16, e_q50, e_q84 = weighted_quantile(expected_grid, weights, [0.16, 0.5, 0.84])
    et_q16, et_q50, et_q84 = weighted_quantile(truncated_expected_grid, weights, [0.16, 0.5, 0.84])
    boundary = bool(best == 0 or best == len(sigmas) - 1 or weights[0] > 0.01 or weights[-1] > 0.01)
    return {
        "sigma_rayleigh": sigma_q50,
        "sigma_rayleigh_lo": sigma_q16,
        "sigma_rayleigh_hi": sigma_q84,
        "sigma_map": sigmas[best],
        "expected_e": e_q50,
        "expected_e_lo": e_q16,
        "expected_e_hi": e_q84,
        "expected_e_map": expected_grid[best],
        "expected_e_truncated": et_q50,
        "expected_e_truncated_lo": et_q16,
        "expected_e_truncated_hi": et_q84,
        "expected_e_truncated_map": truncated_expected_grid[best],
        "ll_max": lls[best],
        "transit_selection_applied": bool(mode != "none"),
        "selection_mode": mode,
        "sigma_at_grid_lower_edge": bool(best == 0),
        "sigma_at_grid_upper_edge": bool(best == len(sigmas) - 1),
        "boundary_flag": boundary,
        "posterior_mass_at_lower_edge": float(weights[0]),
        "posterior_mass_at_upper_edge": float(weights[-1]),
    }


def rayleigh_grid(
    e_grid: np.ndarray,
    sigmas: np.ndarray,
    apply_transit_selection: bool,
    selection_mode: str | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    mode = selection_mode_from_legacy(apply_transit_selection, selection_mode)
    rays = []
    normalizers = []
    for sigma in sigmas:
        ray = (e_grid / sigma**2) * np.exp(-(e_grid**2) / (2.0 * sigma**2))
        ray = ray / trapezoid(ray, e_grid)
        rays.append(ray)
        if mode == "legacy_forward_norm":
            normalizers.append(trapezoid(ray / (1.0 - e_grid**2), e_grid))
        elif mode == "manuscript_reciprocal_with_norm":
            # Mean over omega of 1 / p_transit(e, omega) for the dimensionless
            # transit-probability factor is sqrt(1-e^2). The manuscript's
            # commented equation does not include this population normalizer,
            # so this mode is diagnostic only.
            normalizers.append(trapezoid(ray * np.sqrt(np.clip(1.0 - e_grid**2, 0.0, np.inf)), e_grid))
        else:
            normalizers.append(1.0)
    return np.vstack(rays).T, np.asarray(normalizers, dtype=float)


def posterior_weights_from_ll(sigmas: np.ndarray, lls: np.ndarray) -> np.ndarray:
    widths = grid_cell_widths(sigmas)
    raw = np.exp(lls - np.nanmax(lls)) * widths
    total = raw.sum()
    if not np.isfinite(total) or total <= 0:
        return np.full(len(sigmas), 1.0 / len(sigmas))
    return raw / total


def grid_cell_widths(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    if len(x) == 1:
        return np.ones_like(x)
    edges = np.empty(len(x) + 1, dtype=float)
    edges[1:-1] = 0.5 * (x[:-1] + x[1:])
    edges[0] = max(0.0, x[0] - (edges[1] - x[0]))
    edges[-1] = x[-1] + (x[-1] - edges[-2])
    return np.diff(edges)


def weighted_quantile(values: np.ndarray, weights: np.ndarray, quantiles: list[float]) -> tuple[float, ...]:
    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)
    order = np.argsort(values)
    v = values[order]
    w = weights[order]
    cdf = np.cumsum(w)
    cdf = cdf / cdf[-1]
    return tuple(float(np.interp(q, cdf, v)) for q in quantiles)


def per_planet_leverage(
    summary: pd.DataFrame,
    mass_matrix: np.ndarray,
    e_grid: np.ndarray,
    sigmas: np.ndarray,
    fit: dict,
    apply_transit_selection: bool,
    population: str,
    selection_mode: str | None = None,
) -> pd.DataFrame:
    mode = selection_mode_from_legacy(apply_transit_selection, selection_mode)
    rays, normalizers = rayleigh_grid(e_grid, sigmas, apply_transit_selection, mode)
    terms = mass_matrix @ rays
    if mode in {"legacy_forward_norm", "manuscript_reciprocal_with_norm"}:
        terms = terms / normalizers
    log_terms = np.log(np.clip(terms, 1e-300, None))
    map_idx = int(np.argmin(np.abs(sigmas - float(fit["sigma_map"]))))
    low_idx = 0
    out = summary[["kepoi_name", "koi_target", "kepid", "disk", "system", "e50", "zeta_median"]].copy()
    out["population"] = population
    out["loglike_at_sigma_min"] = log_terms[:, low_idx]
    out["loglike_at_sigma_map"] = log_terms[:, map_idx]
    out["delta_loglike_map_minus_min"] = out["loglike_at_sigma_map"] - out["loglike_at_sigma_min"]
    return out.sort_values("delta_loglike_map_minus_min", ascending=False)


def leaveout_diagnostics(
    summary: pd.DataFrame,
    mass_matrix: np.ndarray,
    e_grid: np.ndarray,
    sigmas: np.ndarray,
    apply_transit_selection: bool,
    population: str,
    rng: np.random.Generator,
    trials: int,
    selection_mode: str | None = None,
) -> list[dict[str, object]]:
    rows = []
    n = len(summary)
    remove_n = max(1, int(round(0.10 * n)))
    mode = selection_mode_from_legacy(apply_transit_selection, selection_mode)
    full = fit_from_mass_matrix(mass_matrix, e_grid, sigmas, apply_transit_selection, mode)
    for trial in range(trials):
        drop = set(rng.choice(np.arange(n), size=remove_n, replace=False))
        keep = np.array([i not in drop for i in range(n)])
        fit = fit_from_mass_matrix(mass_matrix[keep], e_grid, sigmas, apply_transit_selection, mode)
        frac_shift = (fit["expected_e"] - full["expected_e"]) / full["expected_e"] if full["expected_e"] else np.nan
        rows.append(
            {
                "population": population,
                "trial": trial,
                "n_full": n,
                "n_removed": remove_n,
                "expected_e_full": full["expected_e"],
                "expected_e_leaveout": fit["expected_e"],
                "fractional_shift": frac_shift,
                "passes_sagear_5pct": bool(np.isfinite(frac_shift) and abs(frac_shift) <= 0.05),
            }
        )
    return rows


def topk_diagnostics(
    summary: pd.DataFrame,
    mass_matrix: np.ndarray,
    e_grid: np.ndarray,
    sigmas: np.ndarray,
    apply_transit_selection: bool,
    population: str,
    leverage: pd.DataFrame,
    selection_mode: str | None = None,
) -> list[dict[str, object]]:
    rows = []
    mode = selection_mode_from_legacy(apply_transit_selection, selection_mode)
    full = fit_from_mass_matrix(mass_matrix, e_grid, sigmas, apply_transit_selection, mode)
    order = leverage.reset_index().sort_values("delta_loglike_map_minus_min", ascending=False)["index"].to_numpy(int)
    for k in [0, 1, 3, 5, 10, 20, 50, 100]:
        if k >= len(summary):
            continue
        keep = np.ones(len(summary), dtype=bool)
        if k:
            keep[order[:k]] = False
        fit = fit_from_mass_matrix(mass_matrix[keep], e_grid, sigmas, apply_transit_selection, mode)
        frac_shift = (fit["expected_e"] - full["expected_e"]) / full["expected_e"] if full["expected_e"] else np.nan
        rows.append(
            {
                "population": population,
                "removed_top_leverage_n": k,
                "n_remaining": int(keep.sum()),
                "expected_e_full": full["expected_e"],
                "expected_e_after_removal": fit["expected_e"],
                "fractional_shift": frac_shift,
                "passes_sagear_5pct": bool(np.isfinite(frac_shift) and abs(frac_shift) <= 0.05),
            }
        )
    return rows


def transit_probability_weight(e_grid: np.ndarray, omega_grid: np.ndarray) -> np.ndarray:
    e, omega = np.meshgrid(e_grid, omega_grid, indexing="ij")
    weight = (1.0 + e * np.sin(omega)) / (1.0 - e**2)
    return np.clip(weight, 1e-12, np.inf)


def make_ecdf_plot(summary: pd.DataFrame, path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), constrained_layout=True)
    colors = {
        ("thin", "single"): "#008b8b",
        ("thick", "single"): "#9b0000",
        ("thin", "multi"): "#76b7b2",
        ("thick", "multi"): "#c76f6f",
    }
    for disk, system, label in POPULATIONS:
        sub = summary[(summary["disk"] == disk) & (summary["system"] == system)]
        color = colors[(disk, system)]
        for ax, col in [(axes[0], "e50"), (axes[1], "zeta_median")]:
            vals = pd.to_numeric(sub[col], errors="coerce").dropna().sort_values().to_numpy()
            if len(vals) == 0:
                continue
            y = np.arange(1, len(vals) + 1) / len(vals)
            ax.plot(vals, y, color=color, label=f"{label} (n={len(vals)})")
    axes[0].set_xlabel("posterior median e")
    axes[0].set_ylabel("cumulative fraction")
    axes[0].set_xlim(0, 0.95)
    axes[1].set_xlabel("median zeta")
    axes[1].set_xlim(0, 3)
    for ax in axes:
        ax.grid(alpha=0.2)
    axes[0].legend(fontsize=8)
    fig.savefig(path, dpi=220)
    plt.close(fig)


if __name__ == "__main__":
    main()
