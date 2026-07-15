from __future__ import annotations

"""Host-clustered robustness checks for hierarchical Rayleigh fits.

Planets in one transiting system share stellar-density information and fitted
system parameters.  This module complements planet-level leave-out checks by
resampling or removing complete KOI systems.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from common import output_dir
from hierarchical_rayleigh import (
    POPULATIONS,
    fit_from_mass_matrix,
    load_population_masses,
    posterior_weights_from_ll,
    rayleigh_grid,
    recompute_grid_support_flags,
    validate_summary_contract,
)


def expected_e_from_log_likelihood(sigmas: np.ndarray, log_likelihood: np.ndarray) -> float:
    weights = posterior_weights_from_ll(sigmas, log_likelihood)
    expected = sigmas * np.sqrt(np.pi / 2.0)
    order = np.argsort(expected)
    cdf = np.cumsum(weights[order])
    cdf /= cdf[-1]
    return float(np.interp(0.5, cdf, expected[order]))


def grouped_log_terms(
    mass_matrix: np.ndarray,
    e_grid: np.ndarray,
    sigmas: np.ndarray,
    groups: pd.Series,
    selection_mode: str,
) -> tuple[np.ndarray, np.ndarray]:
    rays, normalizers = rayleigh_grid(e_grid, sigmas, True, selection_mode)
    terms = mass_matrix @ rays
    if selection_mode in {"legacy_forward_norm", "manuscript_reciprocal_with_norm"}:
        terms = terms / normalizers
    log_terms = np.log(np.clip(terms, 1e-300, None))
    names = groups.astype(str).to_numpy()
    unique = pd.unique(names)
    grouped = np.vstack([log_terms[names == name].sum(axis=0) for name in unique])
    return np.asarray(unique, dtype=str), grouped


def analyze_population(
    summary: pd.DataFrame,
    mass_matrix: np.ndarray,
    e_grid: np.ndarray,
    sigmas: np.ndarray,
    *,
    population: str,
    selection_mode: str,
    trials: int,
    seed: int,
) -> tuple[dict[str, object], pd.DataFrame, pd.DataFrame]:
    if "koi_target" not in summary:
        raise ValueError("Host-clustered diagnostics require koi_target")
    host_names, host_ll = grouped_log_terms(
        mass_matrix, e_grid, sigmas, summary["koi_target"], selection_mode
    )
    full_ll = host_ll.sum(axis=0)
    full_e = expected_e_from_log_likelihood(sigmas, full_ll)
    rng = np.random.default_rng(seed)
    n_hosts = len(host_names)
    remove_n = max(1, int(round(0.10 * n_hosts)))

    bootstrap_e = np.empty(trials, dtype=float)
    leaveout_e = np.empty(trials, dtype=float)
    trial_rows: list[dict[str, object]] = []
    for trial in range(trials):
        sampled = rng.integers(0, n_hosts, size=n_hosts)
        bootstrap_e[trial] = expected_e_from_log_likelihood(sigmas, host_ll[sampled].sum(axis=0))
        dropped = np.sort(rng.choice(n_hosts, size=remove_n, replace=False))
        kept = np.ones(n_hosts, dtype=bool)
        kept[dropped] = False
        leaveout_e[trial] = expected_e_from_log_likelihood(sigmas, host_ll[kept].sum(axis=0))
        shift = (leaveout_e[trial] - full_e) / full_e if full_e else np.nan
        trial_rows.append(
            {
                "population": population,
                "trial": trial,
                "n_hosts": n_hosts,
                "n_hosts_removed": remove_n,
                "expected_e_full": full_e,
                "expected_e_leaveout": leaveout_e[trial],
                "fractional_shift": shift,
                "passes_5pct": bool(np.isfinite(shift) and abs(shift) <= 0.05),
            }
        )

    leverage_rows: list[dict[str, object]] = []
    for index, host in enumerate(host_names):
        without = expected_e_from_log_likelihood(sigmas, full_ll - host_ll[index])
        shift = (without - full_e) / full_e if full_e else np.nan
        leverage_rows.append(
            {
                "population": population,
                "koi_target": host,
                "n_planets": int((summary["koi_target"].astype(str) == host).sum()),
                "expected_e_full": full_e,
                "expected_e_without_host": without,
                "fractional_shift": shift,
                "absolute_fractional_shift": abs(shift),
            }
        )
    leverage = pd.DataFrame(leverage_rows).sort_values(
        "absolute_fractional_shift", ascending=False
    )

    q16, q50, q84 = np.quantile(bootstrap_e, [0.16, 0.5, 0.84])
    leaveout_shift = (leaveout_e - full_e) / full_e
    row = {
        "population": population,
        "n_planets": len(summary),
        "n_hosts": n_hosts,
        "expected_e_full": full_e,
        "host_bootstrap_e16": q16,
        "host_bootstrap_e50": q50,
        "host_bootstrap_e84": q84,
        "leave10_hosts_trials": trials,
        "leave10_hosts_pass_fraction_5pct": float(np.mean(np.abs(leaveout_shift) <= 0.05)),
        "leave10_hosts_p95_absolute_fractional_shift": float(np.quantile(np.abs(leaveout_shift), 0.95)),
        "leave10_hosts_max_absolute_fractional_shift": float(np.max(np.abs(leaveout_shift))),
        "largest_single_host_absolute_fractional_shift": float(
            leverage["absolute_fractional_shift"].max()
        ),
        "largest_single_host": str(leverage.iloc[0]["koi_target"]),
        "robust_under_5pct_contract": bool(np.all(np.abs(leaveout_shift) <= 0.05)),
    }
    return row, leverage, pd.DataFrame(trial_rows)


def run(
    summary: pd.DataFrame,
    *,
    selection_mode: str = "manuscript_reciprocal",
    trials: int = 1000,
    seed: int = 20260713,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    validate_summary_contract(summary)
    summary = recompute_grid_support_flags(summary)
    summary = summary[~summary["qc_primary_exclude"].fillna(True).astype(bool)].reset_index(drop=True)
    sigmas = np.linspace(1e-4, 1.0, 2000)
    rows: list[dict[str, object]] = []
    leverage: list[pd.DataFrame] = []
    trial_tables: list[pd.DataFrame] = []
    for offset, (disk, system, label) in enumerate(POPULATIONS):
        sub = summary[(summary["disk"] == disk) & (summary["system"] == system)].reset_index(drop=True)
        if len(sub) < 5:
            continue
        masses, e_grid = load_population_masses(sub, True, selection_mode)
        direct = fit_from_mass_matrix(masses, e_grid, sigmas, True, selection_mode)
        row, lev, trial_table = analyze_population(
            sub,
            masses,
            e_grid,
            sigmas,
            population=label,
            selection_mode=selection_mode,
            trials=trials,
            seed=seed + offset,
        )
        row["direct_fit_expected_e"] = direct["expected_e"]
        rows.append(row)
        leverage.append(lev)
        trial_tables.append(trial_table)
    return (
        pd.DataFrame(rows),
        pd.concat(leverage, ignore_index=True),
        pd.concat(trial_tables, ignore_index=True),
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--selection-mode", default="manuscript_reciprocal")
    parser.add_argument("--trials", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20260713)
    parser.add_argument("--out-tag", default="system_robustness")
    parser.add_argument("--output-dir", default=None)
    args = parser.parse_args()
    summary, leverage, trials = run(
        pd.read_csv(Path(args.summary)),
        selection_mode=args.selection_mode,
        trials=args.trials,
        seed=args.seed,
    )
    out = Path(args.output_dir).expanduser().resolve() if args.output_dir else output_dir()
    out.mkdir(parents=True, exist_ok=True)
    summary.to_csv(out / f"rayleigh_host_cluster_summary_{args.out_tag}.csv", index=False)
    leverage.to_csv(out / f"rayleigh_host_leverage_{args.out_tag}.csv", index=False)
    trials.to_csv(out / f"rayleigh_leave10hosts_{args.out_tag}.csv", index=False)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
