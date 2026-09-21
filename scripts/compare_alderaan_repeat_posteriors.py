"""Compare overlapping original and recovery ALDERAAN shape posteriors.

This is a repeatability diagnostic for saved posterior FITS. It does not imply
that the two fits used identical light curves, settings, code, or priors.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits


PARAMETERS = ("ROR", "IMPACT", "DUR14")
QUANTILES = (0.16, 0.50, 0.84)
PARAMETER_RE = re.compile(r"^(ROR|IMPACT|DUR14)_(\d+)$")
PLANET_COLUMNS = [
    "target",
    "planet_index",
    "npl",
    "qc_selected",
    "qc_selected_planet_id",
    "fixed_summary_population",
    "fixed_summary_disk",
    "fixed_summary_system",
    "fixed_summary_posterior_source",
    "fixed_summary_match_status",
    "source_stratum",
    "original_sample_count",
    "recovery_sample_count",
    "original_effective_sample_size",
    "recovery_effective_sample_size",
]
for _parameter in PARAMETERS:
    _lower = _parameter.lower()
    PLANET_COLUMNS.extend(
        f"{prefix}_{_lower}_{suffix}"
        for prefix in ("original", "recovery")
        for suffix in ("q16", "median", "q84")
    )
    PLANET_COLUMNS.append(f"delta_{_lower}_median_recovery_minus_original")


class ComparisonFailure(ValueError):
    """Expected validation failure with a stable machine-readable code."""

    def __init__(self, code: str, detail: str):
        super().__init__(detail)
        self.code = code
        self.detail = detail


@dataclass
class PosteriorFit:
    path: Path
    target: str
    npl: int
    sample_count: int
    parameter_columns: tuple[str, ...]
    weights: np.ndarray
    values: dict[str, np.ndarray]

    @property
    def effective_sample_size(self) -> float:
        return float(1.0 / np.sum(np.square(self.weights)))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def normalize_log_weights(log_weights: np.ndarray) -> np.ndarray:
    """Normalize one fit's finite Dynesty log weights without resampling."""
    values = np.asarray(log_weights, dtype=float)
    if values.ndim != 1 or values.size == 0:
        raise ComparisonFailure("invalid_ln_wt", "LN_WT must be a nonempty 1-D array")
    if not np.isfinite(values).all():
        raise ComparisonFailure("invalid_ln_wt", "LN_WT contains non-finite values")
    shifted = values - np.max(values)
    weights = np.exp(shifted)
    total = float(weights.sum())
    if not np.isfinite(total) or total <= 0:
        raise ComparisonFailure("invalid_ln_wt", "LN_WT cannot be normalized")
    weights /= total
    return weights


def weighted_quantiles(
    values: np.ndarray, weights: np.ndarray, quantiles: tuple[float, ...] = QUANTILES
) -> np.ndarray:
    """Calculate weighted quantiles while retaining each joint sample row."""
    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)
    if values.ndim != 1 or values.shape != weights.shape:
        raise ComparisonFailure("invalid_parameter_values", "values and weights must be paired 1-D arrays")
    if not np.isfinite(values).all():
        raise ComparisonFailure("invalid_parameter_values", "parameter column contains non-finite values")
    requested = np.asarray(quantiles, dtype=float)
    if np.any((requested < 0) | (requested > 1)):
        raise ValueError("quantiles must lie in [0, 1]")
    order = np.argsort(values, kind="stable")
    ordered_values = values[order]
    cumulative = np.cumsum(weights[order])
    cumulative[-1] = 1.0
    return np.interp(requested, cumulative, ordered_values)


def expected_parameter_columns(npl: int) -> tuple[str, ...]:
    return tuple(f"{parameter}_{index}" for index in range(npl) for parameter in PARAMETERS)


def read_fit(path: Path) -> PosteriorFit:
    path = Path(path)
    try:
        with fits.open(path, memmap=False) as hdul:
            target = str(hdul[0].header.get("TARGET", "")).strip()
            try:
                npl = int(hdul[0].header.get("NPL", -1))
            except (TypeError, ValueError) as exc:
                raise ComparisonFailure("invalid_npl", "NPL is not an integer") from exc
            if npl <= 0:
                raise ComparisonFailure("invalid_npl", f"NPL must be positive, got {npl}")
            if "SAMPLES" not in hdul:
                raise ComparisonFailure("missing_samples_hdu", "missing SAMPLES HDU")
            table = hdul["SAMPLES"].data
            names = tuple(hdul["SAMPLES"].columns.names)
            if "LN_WT" not in names:
                raise ComparisonFailure("missing_ln_wt", "SAMPLES is missing LN_WT")
            parameter_columns = tuple(name for name in names if PARAMETER_RE.match(name))
            expected = expected_parameter_columns(npl)
            if set(parameter_columns) != set(expected):
                missing = sorted(set(expected) - set(parameter_columns))
                extra = sorted(set(parameter_columns) - set(expected))
                raise ComparisonFailure(
                    "invalid_parameter_column_set", f"missing={missing}; extra={extra}"
                )
            weights = normalize_log_weights(np.asarray(table["LN_WT"], dtype=float))
            values = {name: np.asarray(table[name], dtype=float) for name in expected}
    except ComparisonFailure:
        raise
    except Exception as exc:
        raise ComparisonFailure("fits_read_failure", f"{type(exc).__name__}: {exc}") from exc
    return PosteriorFit(
        path=path,
        target=target,
        npl=npl,
        sample_count=len(weights),
        parameter_columns=expected,
        weights=weights,
        values=values,
    )


def validate_pair(expected_target: str, original: PosteriorFit, recovery: PosteriorFit) -> None:
    if original.target != expected_target or recovery.target != expected_target:
        raise ComparisonFailure(
            "target_mismatch",
            f"filename target={expected_target}; original TARGET={original.target}; recovery TARGET={recovery.target}",
        )
    if original.npl != recovery.npl:
        raise ComparisonFailure(
            "npl_mismatch", f"original NPL={original.npl}; recovery NPL={recovery.npl}"
        )
    if set(original.parameter_columns) != set(recovery.parameter_columns):
        raise ComparisonFailure(
            "parameter_column_set_mismatch",
            "original and recovery indexed ROR/IMPACT/DUR14 columns differ",
        )


def _bool_mask(values: pd.Series, label: str) -> pd.Series:
    normalized = values.astype(str).str.strip().str.lower()
    if not normalized.isin(("true", "false")).all():
        raise ValueError(f"fixed summary has invalid {label} values")
    return normalized.eq("true")


def load_fixed_summary(path: Path) -> pd.DataFrame:
    summary = pd.read_csv(path)
    required = {
        "koi_target",
        "kepoi_name",
        "alderaan_planet_index",
        "disk",
        "system",
        "posterior_source",
        "qc_primary_exclude",
    }
    missing = sorted(required - set(summary.columns))
    if missing:
        raise ValueError(f"fixed summary is missing columns: {missing}")
    summary = summary.copy()
    summary["_qc_selected"] = ~_bool_mask(summary["qc_primary_exclude"], "qc_primary_exclude")
    return summary


def summary_context(summary: pd.DataFrame, target: str, npl: int) -> dict:
    rows = summary.loc[summary["koi_target"].astype(str).eq(target)].copy()
    unavailable = {
        "status": "target_absent",
        "selected_by_index": {},
        "population": "",
        "disk": "",
        "system": "",
        "posterior_source": "",
    }
    if rows.empty:
        return unavailable
    selected = rows.loc[rows["_qc_selected"]].copy()
    if selected.empty:
        return {**unavailable, "status": "no_qc_selected_planets"}
    try:
        selected["_index"] = selected["alderaan_planet_index"].astype(int)
    except (TypeError, ValueError):
        return {**unavailable, "status": "ambiguous_invalid_planet_index"}
    if (
        selected["_index"].duplicated().any()
        or selected["kepoi_name"].astype(str).duplicated().any()
        or ((selected["_index"] < 0) | (selected["_index"] >= npl)).any()
    ):
        return {**unavailable, "status": "ambiguous_planet_mapping"}
    fields = {}
    for column in ("disk", "system", "posterior_source"):
        unique = sorted(set(selected[column].dropna().astype(str)))
        if len(unique) != 1:
            return {**unavailable, "status": f"ambiguous_{column}"}
        fields[column] = unique[0]
    return {
        "status": "unambiguous",
        "selected_by_index": dict(zip(selected["_index"], selected["kepoi_name"].astype(str))),
        "population": f"{fields['disk']}_{fields['system']}",
        **fields,
    }


def compare_pair(
    target: str,
    original: PosteriorFit,
    recovery: PosteriorFit,
    context: dict,
    source_stratum: str,
) -> list[dict]:
    validate_pair(target, original, recovery)
    rows = []
    selected = context["selected_by_index"]
    for index in range(original.npl):
        row = {
            "target": target,
            "planet_index": index,
            "npl": original.npl,
            "qc_selected": index in selected if context["status"] == "unambiguous" else pd.NA,
            "qc_selected_planet_id": selected.get(index, ""),
            "fixed_summary_population": context["population"],
            "fixed_summary_disk": context["disk"],
            "fixed_summary_system": context["system"],
            "fixed_summary_posterior_source": context["posterior_source"],
            "fixed_summary_match_status": context["status"],
            "source_stratum": source_stratum,
            "original_sample_count": original.sample_count,
            "recovery_sample_count": recovery.sample_count,
            "original_effective_sample_size": original.effective_sample_size,
            "recovery_effective_sample_size": recovery.effective_sample_size,
        }
        for parameter in PARAMETERS:
            column = f"{parameter}_{index}"
            original_q = weighted_quantiles(original.values[column], original.weights)
            recovery_q = weighted_quantiles(recovery.values[column], recovery.weights)
            lower = parameter.lower()
            for prefix, values in (("original", original_q), ("recovery", recovery_q)):
                row[f"{prefix}_{lower}_q16"] = values[0]
                row[f"{prefix}_{lower}_median"] = values[1]
                row[f"{prefix}_{lower}_q84"] = values[2]
            row[f"delta_{lower}_median_recovery_minus_original"] = recovery_q[1] - original_q[1]
        rows.append(row)
    return rows


def discover_results(root: Path) -> dict[str, list[Path]]:
    found: dict[str, list[Path]] = {}
    for path in sorted(Path(root).rglob("*-results.fits")):
        target = path.name.removesuffix("-results.fits")
        found.setdefault(target, []).append(path.resolve())
    return found


def aggregate_rows(comparison: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "stratification",
        "stratum",
        "target_count",
        "planet_count",
        "median_original_effective_sample_size",
        "median_recovery_effective_sample_size",
        "median_delta_ror_recovery_minus_original",
        "median_delta_impact_recovery_minus_original",
        "median_delta_dur14_recovery_minus_original",
    ]
    if comparison.empty:
        return pd.DataFrame(columns=columns)
    groups: list[tuple[str, str, pd.DataFrame]] = [("overall", "all", comparison)]
    for npl, group in comparison.groupby("npl", sort=True):
        groups.append(("npl", str(npl), group))
    available_selected = comparison.loc[comparison["qc_selected"].notna()]
    for selected, group in available_selected.groupby("qc_selected", sort=True):
        groups.append(("qc_selected", "selected" if bool(selected) else "not_selected", group))
    available_population = comparison.loc[comparison["fixed_summary_population"].astype(str).ne("")]
    for population, group in available_population.groupby("fixed_summary_population", sort=True):
        groups.append(("population", str(population), group))
    for source, group in comparison.groupby("source_stratum", sort=True):
        groups.append(("source", str(source), group))
    output = []
    for stratification, stratum, group in groups:
        output.append(
            {
                "stratification": stratification,
                "stratum": stratum,
                "target_count": group["target"].nunique(),
                "planet_count": len(group),
                "median_original_effective_sample_size": group["original_effective_sample_size"].median(),
                "median_recovery_effective_sample_size": group["recovery_effective_sample_size"].median(),
                "median_delta_ror_recovery_minus_original": group[
                    "delta_ror_median_recovery_minus_original"
                ].median(),
                "median_delta_impact_recovery_minus_original": group[
                    "delta_impact_median_recovery_minus_original"
                ].median(),
                "median_delta_dur14_recovery_minus_original": group[
                    "delta_dur14_median_recovery_minus_original"
                ].median(),
            }
        )
    return pd.DataFrame(output, columns=columns)


def run_comparison(
    original_dir: Path,
    recovery_dir: Path,
    fixed_summary_path: Path,
    output_dir: Path,
    *,
    original_source_label: str = "source_original_archive",
    recovery_source_label: str = "source_recovery_local",
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict]:
    original_dir = Path(original_dir).resolve()
    recovery_dir = Path(recovery_dir).resolve()
    fixed_summary_path = Path(fixed_summary_path).resolve()
    output_dir = Path(output_dir).resolve()
    original_files = discover_results(original_dir)
    recovery_files = discover_results(recovery_dir)
    overlap = sorted(set(original_files) & set(recovery_files))
    summary = load_fixed_summary(fixed_summary_path)
    source_stratum = f"{original_source_label}_vs_{recovery_source_label}"
    audit_rows = []
    planet_rows = []
    for target in overlap:
        original_paths = original_files[target]
        recovery_paths = recovery_files[target]
        audit = {
            "target": target,
            "status": "pending",
            "failure_code": "",
            "failure_detail": "",
            "original_path": "|".join(map(str, original_paths)),
            "recovery_path": "|".join(map(str, recovery_paths)),
            "original_sha256": "",
            "recovery_sha256": "",
            "original_header_target": "",
            "recovery_header_target": "",
            "original_npl": pd.NA,
            "recovery_npl": pd.NA,
            "original_sample_count": pd.NA,
            "recovery_sample_count": pd.NA,
            "original_parameter_columns": "",
            "recovery_parameter_columns": "",
            "fixed_summary_match_status": "not_checked",
        }
        try:
            if len(original_paths) != 1 or len(recovery_paths) != 1:
                raise ComparisonFailure(
                    "duplicate_target_files",
                    f"original files={len(original_paths)}; recovery files={len(recovery_paths)}",
                )
            original_path, recovery_path = original_paths[0], recovery_paths[0]
            audit["original_sha256"] = sha256(original_path)
            audit["recovery_sha256"] = sha256(recovery_path)
            original = read_fit(original_path)
            recovery = read_fit(recovery_path)
            audit.update(
                {
                    "original_header_target": original.target,
                    "recovery_header_target": recovery.target,
                    "original_npl": original.npl,
                    "recovery_npl": recovery.npl,
                    "original_sample_count": original.sample_count,
                    "recovery_sample_count": recovery.sample_count,
                    "original_parameter_columns": "|".join(original.parameter_columns),
                    "recovery_parameter_columns": "|".join(recovery.parameter_columns),
                }
            )
            validate_pair(target, original, recovery)
            context = summary_context(summary, target, original.npl)
            audit["fixed_summary_match_status"] = context["status"]
            rows = compare_pair(target, original, recovery, context, source_stratum)
            planet_rows.extend(rows)
            audit["status"] = "compared"
        except ComparisonFailure as exc:
            audit["status"] = "not_compared"
            audit["failure_code"] = exc.code
            audit["failure_detail"] = exc.detail
        audit_rows.append(audit)

    comparison = pd.DataFrame(planet_rows, columns=PLANET_COLUMNS)
    if not comparison.empty:
        comparison = comparison.sort_values(["target", "planet_index"]).reset_index(drop=True)
    audit_frame = pd.DataFrame(audit_rows).sort_values("target").reset_index(drop=True)
    aggregate = aggregate_rows(comparison)
    output_dir.mkdir(parents=True, exist_ok=True)
    comparison_path = output_dir / "per_planet_comparison.csv"
    aggregate_path = output_dir / "aggregate_shifts.csv"
    audit_path = output_dir / "target_audit.csv"
    comparison.to_csv(comparison_path, index=False)
    aggregate.to_csv(aggregate_path, index=False)
    audit_frame.to_csv(audit_path, index=False)
    failure_counts = (
        audit_frame.loc[audit_frame["failure_code"].ne(""), "failure_code"].value_counts().sort_index().to_dict()
    )
    scope = {
        "diagnostic": "bounded local ALDERAAN saved-posterior repeatability comparison",
        "interpretation_boundary": (
            "This compares saved posterior outputs only. It is not evidence that the runs used "
            "the same light curves, settings, code, priors, or stochastic configuration."
        ),
        "eccentricity_analysis": "not_performed",
        "joint_sample_treatment": "row-paired parameters retained; no resampling; LN_WT normalized separately per fit",
        "quantiles": list(QUANTILES),
        "dur14_unit": "native FITS value (days)",
        "standardized_shift": "not_calculated",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "original_dir": str(original_dir),
            "recovery_dir": str(recovery_dir),
            "fixed_summary": str(fixed_summary_path),
            "fixed_summary_sha256": sha256(fixed_summary_path),
            "original_source_label": original_source_label,
            "recovery_source_label": recovery_source_label,
        },
        "counts": {
            "original_targets": len(original_files),
            "recovery_targets": len(recovery_files),
            "overlapping_targets": len(overlap),
            "compared_targets": int(audit_frame["status"].eq("compared").sum()),
            "not_compared_targets": int(audit_frame["status"].eq("not_compared").sum()),
            "compared_planets": len(comparison),
            "fixed_summary_unambiguous_compared_targets": int(
                audit_frame["fixed_summary_match_status"].eq("unambiguous").sum()
            ),
        },
        "failure_counts": {key: int(value) for key, value in failure_counts.items()},
        "status": "completed_with_recorded_mismatches" if failure_counts else "completed",
        "outputs": {
            path.name: {"sha256": sha256(path), "bytes": path.stat().st_size}
            for path in (comparison_path, aggregate_path, audit_path)
        },
    }
    (output_dir / "scope.json").write_text(json.dumps(scope, indent=2) + "\n", encoding="utf-8")
    return comparison, aggregate, audit_frame, scope


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--original-dir", type=Path, required=True)
    parser.add_argument("--recovery-dir", type=Path, required=True)
    parser.add_argument("--fixed-summary", type=Path, required=True)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "metadata/alderaan_repeatability_20260920",
    )
    parser.add_argument("--original-source-label", default="source_original_archive")
    parser.add_argument("--recovery-source-label", default="source_recovery_local")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    _, _, _, scope = run_comparison(
        args.original_dir,
        args.recovery_dir,
        args.fixed_summary,
        args.output_dir,
        original_source_label=args.original_source_label,
        recovery_source_label=args.recovery_source_label,
    )
    print(json.dumps(scope["counts"], sort_keys=True))
    print(f"status={scope['status']} output_dir={Path(args.output_dir).resolve()}")


if __name__ == "__main__":
    main()
