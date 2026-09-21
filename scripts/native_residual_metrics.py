"""Summarize native standardized photometric residuals without refitting."""

import csv
import json
import math
from pathlib import Path

import numpy as np


METRIC_COLUMNS = [
    "quarter", "n", "mean", "median", "rms", "centered_robust_mad",
    "fraction_absolute_residual_gt_3", "fraction_absolute_residual_gt_5",
    "max_absolute_residual", "largest_one_percent_count",
    "largest_one_percent_squared_residual_fraction", "lag1_pair_count",
    "lag1_pearson_correlation", "lag1_correlation_status",
]


def _as_finite_vector(values, name):
    vector = np.asarray(values, dtype=float)
    if vector.ndim != 1 or not len(vector):
        raise ValueError(f"{name} must be a nonempty one-dimensional array")
    if not np.isfinite(vector).all():
        raise ValueError(f"{name} must be finite")
    return vector


def _point_metrics(residuals):
    residuals = _as_finite_vector(residuals, "residuals")
    n = len(residuals)
    centered_median = float(np.median(residuals))
    squared = np.square(residuals)
    largest_count = int(math.ceil(n * 0.01))
    squared_total = float(squared.sum())
    largest_fraction = None if squared_total == 0.0 else float(
        np.sort(squared)[-largest_count:].sum() / squared_total
    )
    return {
        "n": int(n),
        "mean": float(np.mean(residuals)),
        "median": centered_median,
        "rms": float(np.sqrt(np.mean(squared))),
        "centered_robust_mad": float(1.4826 * np.median(np.abs(residuals - centered_median))),
        "fraction_absolute_residual_gt_3": float(np.mean(np.abs(residuals) > 3.0)),
        "fraction_absolute_residual_gt_5": float(np.mean(np.abs(residuals) > 5.0)),
        "max_absolute_residual": float(np.max(np.abs(residuals))),
        "largest_one_percent_count": largest_count,
        "largest_one_percent_squared_residual_fraction": largest_fraction,
    }


def residual_metrics(residuals, times, exposure):
    """Return point diagnostics and gap-aware adjacent-pair correlation."""
    residuals = _as_finite_vector(residuals, "residuals")
    times = _as_finite_vector(times, "times")
    if len(times) != len(residuals):
        raise ValueError("times and residuals must have the same length")
    if not math.isfinite(exposure) or exposure <= 0:
        raise ValueError("exposure must be positive and finite")

    order = np.argsort(times, kind="stable")
    times = times[order]
    residuals = residuals[order]
    if np.any(np.diff(times) == 0.0):
        raise ValueError("times must not contain duplicates within a quarter")
    metrics = _point_metrics(residuals)

    adjacent = (np.diff(times) > 0.0) & (np.diff(times) <= 1.5 * exposure)
    left = residuals[:-1][adjacent]
    right = residuals[1:][adjacent]
    pair_count = int(len(left))
    if pair_count < 2:
        lag1_correlation = None
        lag1_status = "undefined_inadequate_pairs"
    elif np.std(left) == 0.0 or np.std(right) == 0.0:
        lag1_correlation = None
        lag1_status = "undefined_zero_variance"
    else:
        lag1_correlation = float(np.corrcoef(left, right)[0, 1])
        lag1_status = "ok"

    return {
        **metrics,
        "lag1_pair_count": pair_count,
        "lag1_pearson_correlation": lag1_correlation,
        "lag1_correlation_status": lag1_status,
    }


def write_residual_report(groups, output_dir, *, exposure, sample_index, stored_lnlike):
    """Write quarter and global metrics for (quarter, time, residual) groups."""
    if not groups:
        raise ValueError("residual report requires at least one photometric quarter")
    rows = []
    all_residuals = []
    for quarter, times, residuals in groups:
        if not isinstance(quarter, (int, np.integer)):
            raise ValueError("quarter must be an integer")
        metrics = residual_metrics(residuals, times, exposure)
        rows.append({"quarter": int(quarter), **metrics})
        all_residuals.append(_as_finite_vector(residuals, "residuals"))
    rows.sort(key=lambda row: row["quarter"])

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    quarter_path = output_dir / "quarter_residual_metrics.csv"
    with quarter_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=METRIC_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)

    # Aggregate adjacent pairs quarter-by-quarter so concatenation cannot bridge gaps.
    global_pairs_left = []
    global_pairs_right = []
    for _, times, residuals in groups:
        times = _as_finite_vector(times, "times")
        residuals = _as_finite_vector(residuals, "residuals")
        order = np.argsort(times, kind="stable")
        times, residuals = times[order], residuals[order]
        if np.any(np.diff(times) == 0.0):
            raise ValueError("times must not contain duplicates within a quarter")
        adjacent = (np.diff(times) > 0.0) & (np.diff(times) <= 1.5 * exposure)
        global_pairs_left.append(residuals[:-1][adjacent])
        global_pairs_right.append(residuals[1:][adjacent])
    all_metrics = _point_metrics(np.concatenate(all_residuals))
    left = np.concatenate(global_pairs_left)
    right = np.concatenate(global_pairs_right)
    all_metrics["lag1_pair_count"] = int(len(left))
    if len(left) < 2:
        all_metrics["lag1_pearson_correlation"] = None
        all_metrics["lag1_correlation_status"] = "undefined_inadequate_pairs"
    elif np.std(left) == 0.0 or np.std(right) == 0.0:
        all_metrics["lag1_pearson_correlation"] = None
        all_metrics["lag1_correlation_status"] = "undefined_zero_variance"
    else:
        all_metrics["lag1_pearson_correlation"] = float(np.corrcoef(left, right)[0, 1])
        all_metrics["lag1_correlation_status"] = "ok"
    report = {
        "scope": "Photometric standardized residuals at the single best stored LN_LIKE sample; a single best-fit point, not posterior predictive.",
        "sample_index": int(sample_index),
        "stored_lnlike": float(stored_lnlike),
        "quarters": [row["quarter"] for row in rows],
        "global_metrics": all_metrics,
        "lag1_caveat": "Global lag-1 correlation pools only adjacent pairs within native quarters and within 1.5 native LC exposures; it does not cross quarter or transit gaps.",
    }
    global_path = output_dir / "global_residual_metrics.json"
    global_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return quarter_path, global_path, report
