from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import logsumexp
from sklearn.mixture import GaussianMixture

from common import finite_columns, load_config, output_dir, root_path
from diagnose_sample import build_apogee_velocity_calibration
from toomre_diagnostics import ensure_velocity_aliases


COUNT_KEYS = [
    "thin_singles_planets",
    "thick_singles_planets",
    "thin_multi_planets",
    "thick_multi_planets",
    "thin_multi_hosts",
    "thick_multi_hosts",
]


def main() -> None:
    parser = argparse.ArgumentParser(description="Sweep disk-classifier thresholds against Sagear counts.")
    parser.add_argument("--config", default=None)
    parser.add_argument("--sample", default=None)
    parser.add_argument("--out-prefix", default="classifier_threshold_diagnostics")
    args = parser.parse_args()

    cfg = load_config(args.config)
    out_dir = output_dir()
    sample_path = Path(args.sample) if args.sample else out_dir / "canonical_sample_diagnostic.csv"
    sample = ensure_velocity_aliases(pd.read_csv(sample_path))

    probability_fields = build_probability_fields(sample, cfg)
    rows = []
    detail_rows = []
    thresholds = np.linspace(0.05, 0.95, 181)
    for name, probs, note in probability_fields:
        sweep = sweep_thresholds(sample, probs, thresholds, cfg)
        sweep["classifier"] = name
        sweep["note"] = note
        detail_rows.append(sweep)
        best = sweep.sort_values(["l1_planet_count_delta", "l1_all_count_delta", "threshold"]).iloc[0]
        rows.append(best.to_dict())

    summary = pd.DataFrame(rows)
    details = pd.concat(detail_rows, ignore_index=True) if detail_rows else pd.DataFrame()
    summary_path = out_dir / f"{args.out_prefix}_best.csv"
    detail_path = out_dir / f"{args.out_prefix}_grid.csv"
    summary.to_csv(summary_path, index=False)
    details.to_csv(detail_path, index=False)

    print("=== Best threshold per classifier ===")
    cols = [
        "classifier",
        "threshold",
        "l1_planet_count_delta",
        "l1_all_count_delta",
        "thin_singles_planets",
        "thick_singles_planets",
        "thin_multi_planets",
        "thick_multi_planets",
        "thin_multi_hosts",
        "thick_multi_hosts",
        "note",
    ]
    print(summary[[c for c in cols if c in summary.columns]].to_string(index=False))
    print(f"\nWrote: {summary_path}")
    print(f"Wrote: {detail_path}")


def build_probability_fields(sample: pd.DataFrame, cfg: dict) -> list[tuple[str, np.ndarray, str]]:
    fields: list[tuple[str, np.ndarray, str]] = []
    fields.extend(
        [
            gmm_probability_field(sample, sample, "planet_host_direct_gmm", "V_phi", "V_perp"),
            gmm_probability_field(sample, sample, "planet_host_geom_gmm", "V_phi_geom", "V_perp_geom"),
        ]
    )
    apogee_path = root_path(cfg, "apogee_dr17_chemical")
    if apogee_path and apogee_path.exists():
        apogee = pd.read_csv(apogee_path)
        calibration = build_apogee_velocity_calibration(apogee, "kepid", "feh", "mgfe", cfg)
        calibration = ensure_velocity_aliases(calibration)
        calibration = calibration[finite_columns(calibration, ["apogee_feh", "apogee_mgfe", "V_phi", "V_perp"])]
        high_alpha = calibration["apogee_mgfe"] > (
            cfg["cuts"]["high_alpha_slope"] * calibration["apogee_feh"] + cfg["cuts"]["high_alpha_intercept"]
        )
        fields.append(
            chemical_likelihood_probability_field(
                sample, calibration, high_alpha, "chem_1thin_1thick_prior_0p425", "V_phi", "V_perp", 1, 1, 0.425
            )
        )
        fields.append(
            chemical_likelihood_probability_field(
                sample, calibration, high_alpha, "chem_2thin_2thick_prior_0p350", "V_phi", "V_perp", 2, 2, 0.350
            )
        )
        fields.append(
            chemical_likelihood_probability_field(
                sample, calibration, high_alpha, "chem_2thin_1thick_prior_0p285", "V_phi", "V_perp", 2, 1, 0.285
            )
        )
    return fields


def gmm_probability_field(
    sample: pd.DataFrame,
    train: pd.DataFrame,
    name: str,
    x_col: str,
    y_col: str,
) -> tuple[str, np.ndarray, str]:
    use_train = finite_columns(train, [x_col, y_col])
    use_sample = finite_columns(sample, [x_col, y_col])
    probs = np.full(len(sample), np.nan)
    if use_train.sum() < 20 or use_sample.sum() < 20:
        return name, probs, "skipped_insufficient_finite"
    gmm = GaussianMixture(n_components=2, covariance_type="full", random_state=42, n_init=20)
    gmm.fit(train.loc[use_train, [x_col, y_col]].to_numpy())
    score = -gmm.means_[:, 0] + 0.25 * gmm.means_[:, 1]
    thick_i = int(np.argmax(score))
    probs[use_sample.to_numpy()] = gmm.predict_proba(sample.loc[use_sample, [x_col, y_col]].to_numpy())[:, thick_i]
    note = f"{x_col},{y_col}; means={gmm.means_.round(3).tolist()}; thick_component={thick_i}"
    return name, probs, note


def chemical_likelihood_probability_field(
    sample: pd.DataFrame,
    calibration: pd.DataFrame,
    high_alpha: pd.Series,
    name: str,
    x_col: str,
    y_col: str,
    n_thin: int,
    n_thick: int,
    prior: float,
) -> tuple[str, np.ndarray, str]:
    use_cal = finite_columns(calibration, [x_col, y_col])
    use_sample = finite_columns(sample, [x_col, y_col])
    probs = np.full(len(sample), np.nan)
    high = high_alpha.loc[use_cal]
    cal = calibration.loc[use_cal]
    if high.sum() < 10 or (~high).sum() < 10 or use_sample.sum() < 20:
        return name, probs, "skipped_insufficient_finite"

    thin_g = GaussianMixture(n_components=n_thin, covariance_type="full", random_state=42, n_init=10)
    thick_g = GaussianMixture(n_components=n_thick, covariance_type="full", random_state=42, n_init=10)
    thin_g.fit(cal.loc[~high, [x_col, y_col]].to_numpy())
    thick_g.fit(cal.loc[high, [x_col, y_col]].to_numpy())
    x_sample = sample.loc[use_sample, [x_col, y_col]].to_numpy()
    log_l_thin = thin_g.score_samples(x_sample)
    log_l_thick = thick_g.score_samples(x_sample)
    log_num = log_l_thick + np.log(prior)
    log_den = logsumexp(np.vstack([log_l_thin + np.log(1.0 - prior), log_num]), axis=0)
    probs[use_sample.to_numpy()] = np.exp(log_num - log_den)
    note = (
        f"n_thin={n_thin}; n_thick={n_thick}; prior={prior}; "
        f"thin_means={thin_g.means_.round(3).tolist()}; thick_means={thick_g.means_.round(3).tolist()}"
    )
    return name, probs, note


def sweep_thresholds(sample: pd.DataFrame, probs: np.ndarray, thresholds: np.ndarray, cfg: dict) -> pd.DataFrame:
    rows = []
    for threshold in thresholds:
        rows.append(count_row(sample, probs, threshold, cfg))
    return pd.DataFrame(rows)


def count_row(sample: pd.DataFrame, probs: np.ndarray, threshold: float, cfg: dict) -> dict[str, object]:
    disk = np.where(probs > threshold, "thick", "thin").astype(object)
    disk[~np.isfinite(probs)] = "unclassified"

    def planets(d: str, s: str) -> int:
        return int(((disk == d) & (sample["system"] == s)).sum())

    def hosts(d: str, s: str) -> int:
        return int(sample.loc[(disk == d) & (sample["system"] == s), "kepid"].nunique())

    counts = {
        "thin_singles_planets": planets("thin", "single"),
        "thick_singles_planets": planets("thick", "single"),
        "thin_multi_planets": planets("thin", "multi"),
        "thick_multi_planets": planets("thick", "multi"),
        "thin_multi_hosts": hosts("thin", "multi"),
        "thick_multi_hosts": hosts("thick", "multi"),
    }
    targets = cfg["sagear_targets"]
    deltas = {f"delta_{k}": counts[k] - int(targets[k]) for k in COUNT_KEYS}
    l1_planets = sum(abs(deltas[f"delta_{k}"]) for k in COUNT_KEYS[:4])
    l1_all = sum(abs(deltas[f"delta_{k}"]) for k in COUNT_KEYS)
    return {
        "threshold": float(threshold),
        **counts,
        **deltas,
        "l1_planet_count_delta": int(l1_planets),
        "l1_all_count_delta": int(l1_all),
    }


if __name__ == "__main__":
    main()
