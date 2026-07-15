"""Analyze returned ALDERAAN factorial-validation fits as paired experiments.

The validation bundle contains six arms.  This script discovers their FITS
files, regenerates eccentricity posteriors with the direct MacDougall-style
extractor, compares only the intended paired targets, and calibrates practical
effect flags against the observed repeat-run distribution.

The evidence labels are deliberately descriptive.  Exceeding the observed
repeatability distribution is an exploratory robustness signal, not a formal
hypothesis test or evidence that an arm caused an astrophysical difference.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
import sys
from typing import Callable, Iterable

import numpy as np
import pandas as pd
from astropy.io import fits

from common import load_config, normalize_dynesty_weights
from extract_eccentricity_posteriors import stable_seed
from extract_eccentricity_posteriors_direct import (
    match_planets_by_period,
    process_target as extract_direct_target,
    read_alderaan_planets,
    weighted_quantile,
)


REPO_ROOT = Path(__file__).resolve().parent


@dataclass(frozen=True)
class ArmSpec:
    run_id: str
    target_set: str


@dataclass(frozen=True)
class PairSpec:
    comparison_id: str
    effect: str
    comparison_arm: str
    baseline_arm: str
    target_set: str


ARM_SPECS = {
    "original_lc": ArmSpec("sagear_validation_original_lc", "all"),
    "reference_lc": ArmSpec("sagear_validation_reference_lc", "all"),
    "original_lc_repeat": ArmSpec("sagear_validation_original_lc_repeat", "repeat"),
    "original_lcsc": ArmSpec("sagear_validation_original_lcsc", "sc"),
    "reference_lcsc": ArmSpec("sagear_validation_reference_lcsc", "sc"),
    "paper_priors_original_lc": ArmSpec("sagear_validation_paper_priors_original_lc", "repeat"),
}

PAIR_SPECS = (
    PairSpec("ld_reference_lc", "limb_darkening", "reference_lc", "original_lc", "all"),
    PairSpec(
        "sampler_repeatability",
        "sampler_repeatability",
        "original_lc_repeat",
        "original_lc",
        "repeat",
    ),
    PairSpec("cadence_original_ld", "cadence", "original_lcsc", "original_lc", "sc"),
    PairSpec("cadence_reference_ld", "cadence", "reference_lcsc", "reference_lc", "sc"),
    PairSpec(
        "paper_prior_ambiguity",
        "prior_ambiguity",
        "paper_priors_original_lc",
        "original_lc",
        "repeat",
    ),
)

TARGET_MANIFESTS = {
    "all": "targets_ld_reference_validation.csv",
    "repeat": "targets_repeatability_validation.csv",
    "sc": "targets_short_cadence_validation.csv",
}

PARAMETERS = {
    "t14_hr": ("t14_hr_p16", "t14_hr_p50", "t14_hr_p84", True),
    "impact": ("impact_p16", "impact_p50", "impact_p84", False),
    "rp_over_rs": ("rp_over_rs_p16", "rp_over_rs_p50", "rp_over_rs_p84", True),
    "e": ("e16", "e50", "e84", False),
    "zeta": ("zeta_p16", "zeta_median", "zeta_p84", True),
}


class ValidationAnalysisError(RuntimeError):
    """An input or extraction failure that makes paired analysis invalid."""


class DirectExtractionError(ValidationAnalysisError):
    """Direct extraction failed with an explicit per-planet exclusion manifest."""

    def __init__(self, message: str, exclusions: pd.DataFrame):
        super().__init__(message)
        self.exclusions = exclusions


def resolve_path(path: str | Path) -> Path:
    candidate = Path(path)
    return candidate if candidate.is_absolute() else REPO_ROOT / candidate


def _require_columns(frame: pd.DataFrame, columns: Iterable[str], label: str) -> None:
    missing = sorted(set(columns) - set(frame.columns))
    if missing:
        raise ValidationAnalysisError(f"{label} is missing required columns: {', '.join(missing)}")


def read_target_sets(metadata_root: Path) -> dict[str, set[str]]:
    target_sets: dict[str, set[str]] = {}
    for name, filename in TARGET_MANIFESTS.items():
        path = metadata_root / filename
        if not path.is_file():
            raise ValidationAnalysisError(f"Required target manifest not found: {path}")
        frame = pd.read_csv(path)
        _require_columns(frame, ["koi_target"], str(path))
        targets = set(frame["koi_target"].dropna().astype(str))
        if not targets:
            raise ValidationAnalysisError(f"Target manifest contains no targets: {path}")
        target_sets[name] = targets
    if not target_sets["repeat"].issubset(target_sets["all"]):
        raise ValidationAnalysisError("Repeatability targets are not a subset of the all-target manifest")
    if not target_sets["sc"].issubset(target_sets["all"]):
        raise ValidationAnalysisError("Short-cadence targets are not a subset of the all-target manifest")
    return target_sets


def discover_validation_fits(
    validation_root: Path,
    target_sets: dict[str, set[str]],
) -> pd.DataFrame:
    """Return one audit row for every expected arm/target FITS file."""
    if not validation_root.is_dir():
        raise ValidationAnalysisError(f"Validation root is not a directory: {validation_root}")

    found: dict[tuple[str, str], list[Path]] = {}
    arm_names = set(ARM_SPECS)
    for path in sorted(validation_root.rglob("*-results.fits")):
        matching_arms = arm_names.intersection(path.parts)
        if len(matching_arms) != 1:
            continue
        arm = next(iter(matching_arms))
        target = path.name.removesuffix("-results.fits")
        if target in target_sets[ARM_SPECS[arm].target_set]:
            found.setdefault((arm, target), []).append(path.resolve())

    rows: list[dict[str, object]] = []
    duplicate_messages: list[str] = []
    for arm, spec in ARM_SPECS.items():
        for target in sorted(target_sets[spec.target_set]):
            paths = found.get((arm, target), [])
            if len(paths) > 1:
                duplicate_messages.append(f"{arm}/{target}: " + ", ".join(map(str, paths)))
            path = paths[0] if len(paths) == 1 else None
            rows.append(
                {
                    "arm": arm,
                    "run_id": spec.run_id,
                    "target_set": spec.target_set,
                    "koi_target": target,
                    "status": "present" if path is not None else ("duplicate" if paths else "missing"),
                    "results_file": str(path) if path is not None else "",
                }
            )
    if duplicate_messages:
        raise ValidationAnalysisError(
            "Multiple FITS files matched an expected arm/target:\n  " + "\n  ".join(duplicate_messages)
        )
    return pd.DataFrame(rows)


def build_planet_sample(inventory_path: Path, sample_path: Path, all_targets: set[str]) -> pd.DataFrame:
    inventory = pd.read_csv(inventory_path)
    sample = pd.read_csv(sample_path)
    _require_columns(
        inventory,
        ["koi_target", "kepoi_name", "kepid", "period_days", "included_in_alderaan_system"],
        str(inventory_path),
    )
    _require_columns(sample, ["kepid", "rho_log", "rho_log_upper", "rho_log_lower"], str(sample_path))

    included = inventory["included_in_alderaan_system"]
    if included.dtype != bool:
        included = included.astype(str).str.strip().str.lower().isin({"true", "1", "yes"})
    inventory = inventory[included & inventory["koi_target"].astype(str).isin(all_targets)].copy()
    if inventory.empty:
        raise ValidationAnalysisError(f"No included validation planets found in inventory: {inventory_path}")

    context_columns = ["kepid", "rho_log", "rho_log_upper", "rho_log_lower"]
    context_columns += [column for column in ("disk", "system") if column in sample]
    context = sample[context_columns].drop_duplicates("kepid")
    planets = inventory.merge(context, on="kepid", how="left", validate="many_to_one")
    planets = planets.rename(columns={"period_days": "koi_period"})
    planets["koi_target"] = planets["koi_target"].astype(str)
    planets["kepoi_name"] = planets["kepoi_name"].astype(str)
    for column in ("kepid", "koi_period", "rho_log", "rho_log_upper", "rho_log_lower"):
        planets[column] = pd.to_numeric(planets[column], errors="coerce")

    bad = planets[
        planets[["kepid", "koi_period", "rho_log"]].isna().any(axis=1)
        | (planets["koi_period"] <= 0)
    ]
    if not bad.empty:
        labels = ", ".join(bad["kepoi_name"].head(10))
        raise ValidationAnalysisError(
            f"Inventory planets lack a finite KIC, period, or stellar density in {sample_path}: {labels}"
        )
    return planets.sort_values(["koi_target", "koi_period"]).reset_index(drop=True)


def _native_sample_frame(data: fits.FITS_rec) -> pd.DataFrame:
    array = np.array(data)
    if not array.dtype.isnative:
        array = array.byteswap().view(array.dtype.newbyteorder("="))
    return pd.DataFrame(array)


def extract_shape_rows(
    results_file: Path,
    planets: pd.DataFrame,
    arm: str,
    period_tol: float,
) -> list[dict[str, object]]:
    with fits.open(results_file, memmap=False) as hdul:
        if "SAMPLES" not in hdul:
            raise ValidationAnalysisError(f"Missing SAMPLES HDU: {results_file}")
        samples = _native_sample_frame(hdul["SAMPLES"].data)
        fit_planets = read_alderaan_planets(hdul)
    if "LN_WT" not in samples:
        raise ValidationAnalysisError(f"Missing LN_WT sample column: {results_file}")

    matches = match_planets_by_period(planets, fit_planets, period_tol)
    if len(matches) != len(planets):
        matched = {row_index for row_index, *_ in matches}
        missing = [str(planets.iloc[index]["kepoi_name"]) for index in range(len(planets)) if index not in matched]
        raise ValidationAnalysisError(
            f"Shape period matching failed for {arm} {results_file.name}: {', '.join(missing)}"
        )

    weights = normalize_dynesty_weights(samples["LN_WT"].to_numpy(float))
    rows: list[dict[str, object]] = []
    for planet_index, fit_index, fit_period, relative_difference in matches:
        planet = planets.iloc[planet_index]
        columns = {
            "t14_hr": f"DUR14_{fit_index}",
            "rp_over_rs": f"ROR_{fit_index}",
            "impact": f"IMPACT_{fit_index}",
        }
        missing_columns = [column for column in columns.values() if column not in samples]
        if missing_columns:
            raise ValidationAnalysisError(
                f"Missing paired shape columns in {results_file}: {', '.join(missing_columns)}"
            )
        values = {
            "t14_hr": samples[columns["t14_hr"]].to_numpy(float) * 24.0,
            "rp_over_rs": samples[columns["rp_over_rs"]].to_numpy(float),
            "impact": samples[columns["impact"]].to_numpy(float),
        }
        row: dict[str, object] = {
            "arm": arm,
            "koi_target": planet["koi_target"],
            "kepoi_name": planet["kepoi_name"],
            "kepid": int(planet["kepid"]),
            "koi_period": float(planet["koi_period"]),
            "disk": planet.get("disk", np.nan),
            "system": planet.get("system", np.nan),
            "fit_planet_index": fit_index,
            "fit_period_days": fit_period,
            "period_relative_difference": relative_difference,
            "nested_sample_count": len(samples),
            "results_file": str(results_file),
        }
        for parameter, array in values.items():
            q16, q50, q84 = weighted_quantile(array, weights, [0.16, 0.5, 0.84])
            row[f"{parameter}_p16"] = q16
            row[f"{parameter}_p50"] = q50
            row[f"{parameter}_p84"] = q84
        rows.append(row)
    return rows


def zeta_quantiles(posterior_file: Path) -> tuple[float, float, float]:
    with np.load(posterior_file, allow_pickle=False) as posterior_data:
        required = {"e_grid", "omega_grid", "posterior"}
        missing = required - set(posterior_data.files)
        if missing:
            raise ValidationAnalysisError(
                f"Direct posterior is missing {', '.join(sorted(missing))}: {posterior_file}"
            )
        e_grid = np.asarray(posterior_data["e_grid"], dtype=float)
        omega_grid = np.asarray(posterior_data["omega_grid"], dtype=float)
        weights = np.asarray(posterior_data["posterior"], dtype=float)
    if weights.shape != (len(e_grid), len(omega_grid)) or not np.isfinite(weights).all() or weights.sum() <= 0:
        raise ValidationAnalysisError(f"Direct posterior grid is invalid: {posterior_file}")
    eccentricity, omega = np.meshgrid(e_grid, omega_grid, indexing="ij")
    with np.errstate(invalid="ignore", divide="ignore"):
        zeta = np.sqrt(1.0 - eccentricity**2) / (1.0 + eccentricity * np.sin(omega))
    quantiles = weighted_quantile(zeta.ravel(), weights.ravel(), [0.16, 0.5, 0.84])
    if not np.isfinite(quantiles).all():
        raise ValidationAnalysisError(f"Finite zeta summaries could not be produced: {posterior_file}")
    return tuple(float(value) for value in quantiles)


def extract_all_arms(
    discovery: pd.DataFrame,
    planets: pd.DataFrame,
    direct_root: Path,
    *,
    e_grid: np.ndarray,
    omega_grid: np.ndarray,
    n_proposals: int,
    e_max: float,
    density_error_mode: str,
    period_tol: float,
    min_importance_ess: float,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    shape_rows: list[dict[str, object]] = []
    direct_rows: list[dict[str, object]] = []
    exclusions: list[dict[str, object]] = []

    for audit in discovery.sort_values(["arm", "koi_target"]).itertuples(index=False):
        arm = str(audit.arm)
        target = str(audit.koi_target)
        target_planets = planets[planets["koi_target"].eq(target)].sort_values("koi_period").reset_index(drop=True)
        if target_planets.empty:
            raise ValidationAnalysisError(f"No inventory planets found for expected target {target}")
        results_file = Path(str(audit.results_file))
        shape_rows.extend(extract_shape_rows(results_file, target_planets, arm, period_tol))

        arm_out = direct_root / arm
        arm_out.mkdir(parents=True, exist_ok=True)
        try:
            summaries, target_exclusions = extract_direct_target(
                results_file,
                target_planets,
                arm_out,
                e_grid,
                omega_grid,
                n_proposals=n_proposals,
                e_max=e_max,
                density_error_mode=density_error_mode,
                period_tol=period_tol,
                min_importance_ess=min_importance_ess,
                allow_density_error_fallback=False,
            )
        except Exception as exc:
            raise ValidationAnalysisError(
                f"Direct eccentricity extraction crashed for {arm}/{target} ({results_file}): {exc}"
            ) from exc
        for exclusion in target_exclusions:
            exclusions.append({"arm": arm, **exclusion})
        for summary in summaries:
            p16, median, p84 = zeta_quantiles(Path(str(summary["posterior_file"])))
            direct_rows.append(
                {
                    "arm": arm,
                    **summary,
                    "zeta_p16": p16,
                    "zeta_median": median,
                    "zeta_p84": p84,
                }
            )

    exclusion_frame = pd.DataFrame(
        exclusions,
        columns=["arm", "koi_target", "kepid", "kepoi_name", "koi_period", "results_file", "stage", "reason"],
    )
    if exclusions:
        examples = "; ".join(
            f"{row['arm']}/{row.get('kepoi_name')}: {row.get('stage')} - {row.get('reason')}"
            for row in exclusions[:8]
        )
        raise DirectExtractionError(
            "Direct eccentricity outputs could not be produced for every expected planet. " + examples,
            exclusion_frame,
        )
    if not direct_rows:
        raise ValidationAnalysisError("Direct eccentricity extraction produced no planet summaries")

    shape = pd.DataFrame(shape_rows)
    direct = pd.DataFrame(direct_rows)
    keys = ["arm", "koi_target", "kepoi_name"]
    direct_keep = keys + [
        "e16",
        "e50",
        "e84",
        "zeta_p16",
        "zeta_median",
        "zeta_p84",
        "importance_ess",
        "valid_proposal_fraction",
        "minimum_importance_ess",
        "qc_importance_ess_low",
        "qc_primary_exclude",
        "qc_reasons",
        "posterior_file",
    ]
    planets_out = shape.merge(direct[direct_keep], on=keys, how="outer", validate="one_to_one", indicator=True)
    unmatched = planets_out[planets_out["_merge"].ne("both")]
    if not unmatched.empty:
        labels = ", ".join(
            f"{row['arm']}/{row['kepoi_name']}:{row['_merge']}" for _, row in unmatched.iterrows()
        )
        raise ValidationAnalysisError(f"Shape/direct outputs did not match one-to-one: {labels}")
    return planets_out.drop(columns="_merge"), exclusion_frame


def build_paired_outputs(
    arm_planets: pd.DataFrame,
    target_sets: dict[str, set[str]],
    *,
    allow_incomplete: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    identifiers = ["koi_target", "kepoi_name", "kepid", "koi_period", "disk", "system"]
    summary_columns = sorted({column for spec in PARAMETERS.values() for column in spec[:3]})
    quality_columns = [
        "importance_ess",
        "valid_proposal_fraction",
        "qc_importance_ess_low",
        "qc_primary_exclude",
        "qc_reasons",
    ]
    wide_rows: list[pd.DataFrame] = []
    long_rows: list[dict[str, object]] = []

    for pair in PAIR_SPECS:
        targets = target_sets[pair.target_set]
        baseline = arm_planets[
            arm_planets["arm"].eq(pair.baseline_arm) & arm_planets["koi_target"].isin(targets)
        ][identifiers + summary_columns + quality_columns]
        comparison = arm_planets[
            arm_planets["arm"].eq(pair.comparison_arm) & arm_planets["koi_target"].isin(targets)
        ][identifiers + summary_columns + quality_columns]
        paired = baseline.merge(
            comparison,
            on=["koi_target", "kepoi_name"],
            how="outer" if not allow_incomplete else "inner",
            suffixes=("_baseline", "_comparison"),
            validate="one_to_one",
            indicator=True,
        )
        unmatched = paired[paired["_merge"].ne("both")]
        if not allow_incomplete and not unmatched.empty:
            labels = ", ".join(
                f"{row['koi_target']}/{row['kepoi_name']}:{row['_merge']}" for _, row in unmatched.iterrows()
            )
            raise ValidationAnalysisError(f"Incomplete paired planets for {pair.comparison_id}: {labels}")
        paired = paired.drop(columns="_merge")
        if paired.empty:
            continue
        paired.insert(0, "comparison_id", pair.comparison_id)
        paired.insert(1, "effect", pair.effect)
        paired.insert(2, "comparison_arm", pair.comparison_arm)
        paired.insert(3, "baseline_arm", pair.baseline_arm)

        for parameter, (low, center, high, fractional) in PARAMETERS.items():
            baseline_center = pd.to_numeric(paired[f"{center}_baseline"], errors="coerce")
            comparison_center = pd.to_numeric(paired[f"{center}_comparison"], errors="coerce")
            delta = comparison_center - baseline_center
            paired[f"delta_{parameter}"] = delta
            if fractional:
                paired[f"fractional_delta_{parameter}"] = delta / baseline_center.replace(0.0, np.nan)
            for index, row in paired.iterrows():
                long_rows.append(
                    {
                        "comparison_id": pair.comparison_id,
                        "effect": pair.effect,
                        "comparison_arm": pair.comparison_arm,
                        "baseline_arm": pair.baseline_arm,
                        "koi_target": row["koi_target"],
                        "kepoi_name": row["kepoi_name"],
                        "kepid": row["kepid_baseline"],
                        "koi_period": row["koi_period_baseline"],
                        "disk": row["disk_baseline"],
                        "system": row["system_baseline"],
                        "parameter": parameter,
                        "baseline_p16": row[f"{low}_baseline"],
                        "baseline_p50": row[f"{center}_baseline"],
                        "baseline_p84": row[f"{high}_baseline"],
                        "comparison_p16": row[f"{low}_comparison"],
                        "comparison_p50": row[f"{center}_comparison"],
                        "comparison_p84": row[f"{high}_comparison"],
                        "delta": delta.loc[index],
                        "abs_delta": abs(delta.loc[index]),
                        "fractional_delta": (
                            delta.loc[index] / baseline_center.loc[index]
                            if fractional and baseline_center.loc[index] != 0
                            else np.nan
                        ),
                        "baseline_direct_qc_exclude": bool(row["qc_primary_exclude_baseline"]),
                        "comparison_direct_qc_exclude": bool(row["qc_primary_exclude_comparison"]),
                    }
                )
        wide_rows.append(paired)
    if not wide_rows:
        raise ValidationAnalysisError("No complete arm/planet pairs were available for comparison")
    return pd.concat(wide_rows, ignore_index=True), pd.DataFrame(long_rows)


def system_cluster_bootstrap(
    frame: pd.DataFrame,
    value_column: str,
    statistic: Callable[[np.ndarray], float],
    *,
    n_bootstrap: int,
    seed: int,
) -> np.ndarray:
    """Resample systems, retaining all planet rows within each sampled system."""
    clean = frame[["koi_target", value_column]].copy()
    clean[value_column] = pd.to_numeric(clean[value_column], errors="coerce")
    clean = clean[np.isfinite(clean[value_column])]
    systems = sorted(clean["koi_target"].astype(str).unique())
    if not systems:
        return np.full(n_bootstrap, np.nan)
    grouped = {
        system: clean.loc[clean["koi_target"].astype(str).eq(system), value_column].to_numpy(float)
        for system in systems
    }
    rng = np.random.default_rng(seed)
    output = np.empty(n_bootstrap, dtype=float)
    for index in range(n_bootstrap):
        sampled = rng.choice(systems, size=len(systems), replace=True)
        values = np.concatenate([grouped[system] for system in sampled])
        output[index] = statistic(values)
    return output


def _direct_qc_mask(frame: pd.DataFrame) -> pd.Series:
    baseline = frame.get("baseline_direct_qc_exclude", pd.Series(False, index=frame.index)).fillna(True).astype(bool)
    comparison = frame.get("comparison_direct_qc_exclude", pd.Series(False, index=frame.index)).fillna(True).astype(bool)
    return baseline | comparison


def _finite_quantile(values: np.ndarray, quantile: float) -> float:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    return float(np.quantile(finite, quantile)) if len(finite) else np.nan


def attach_repeatability_evidence(metrics: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    direct_qc = _direct_qc_mask(metrics)
    repeat = metrics[metrics["comparison_id"].eq("sampler_repeatability") & ~direct_qc]
    threshold_rows: list[dict[str, object]] = []
    for parameter in PARAMETERS:
        values = pd.to_numeric(repeat.loc[repeat["parameter"].eq(parameter), "abs_delta"], errors="coerce")
        values = values[np.isfinite(values)]
        threshold_rows.append(
            {
                "parameter": parameter,
                "repeatability_n_planets": len(values),
                "repeatability_abs_delta_p50": float(np.quantile(values, 0.5)) if len(values) else np.nan,
                "repeatability_abs_delta_p95": float(np.quantile(values, 0.95)) if len(values) else np.nan,
                "repeatability_abs_delta_max": float(np.max(values)) if len(values) else np.nan,
            }
        )
    thresholds = pd.DataFrame(threshold_rows)
    out = metrics.merge(thresholds, on="parameter", how="left", validate="many_to_one")
    threshold = out["repeatability_abs_delta_p95"]
    out["repeatability_ratio"] = np.where(
        threshold > 0,
        out["abs_delta"] / threshold,
        np.where(out["abs_delta"].eq(0) & threshold.eq(0), 0.0, np.nan),
    )
    out["practical_evidence"] = "within_observed_repeatability_p95"
    out.loc[threshold.isna(), "practical_evidence"] = "not_assessable_without_repeatability"
    out.loc[out["abs_delta"] > threshold, "practical_evidence"] = (
        "exceeds_observed_repeatability_p95_exploratory"
    )
    out.loc[out["comparison_id"].eq("sampler_repeatability"), "practical_evidence"] = (
        "repeatability_reference_distribution"
    )
    output_direct_qc = _direct_qc_mask(out)
    out.loc[output_direct_qc, "practical_evidence"] = "not_assessable_direct_qc"
    return out, thresholds


def summarize_arms(metrics: pd.DataFrame, *, n_bootstrap: int, seed: int) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for (comparison_id, parameter), group in metrics.groupby(["comparison_id", "parameter"], sort=True):
        direct_qc = _direct_qc_mask(group)
        finite = group[
            np.isfinite(pd.to_numeric(group["delta"], errors="coerce")) & ~direct_qc
        ].copy()
        pair = next(spec for spec in PAIR_SPECS if spec.comparison_id == comparison_id)
        signed_bootstrap = system_cluster_bootstrap(
            finite,
            "delta",
            np.median,
            n_bootstrap=n_bootstrap,
            seed=stable_seed("factorial-validation", seed, comparison_id, parameter, "signed"),
        )
        absolute_bootstrap = system_cluster_bootstrap(
            finite,
            "abs_delta",
            np.median,
            n_bootstrap=n_bootstrap,
            seed=stable_seed("factorial-validation", seed, comparison_id, parameter, "absolute"),
        )
        threshold = finite["repeatability_abs_delta_p95"].iloc[0] if len(finite) else np.nan
        n_exceeding = int((finite["abs_delta"] > threshold).sum()) if np.isfinite(threshold) else 0
        if comparison_id == "sampler_repeatability":
            evidence = "repeatability_reference_distribution"
        elif not np.isfinite(threshold):
            evidence = "not_assessable_without_repeatability"
        elif n_exceeding:
            evidence = "some_shifts_exceed_repeatability_p95_exploratory"
        else:
            evidence = "no_shift_exceeds_repeatability_p95"
        n_systems = int(finite["koi_target"].nunique())
        sample_adequacy = (
            "underpowered_fewer_than_5_systems"
            if n_systems < 5
            else "descriptive_targeted_validation_sample"
        )
        rows.append(
            {
                "comparison_id": comparison_id,
                "effect": pair.effect,
                "comparison_arm": pair.comparison_arm,
                "baseline_arm": pair.baseline_arm,
                "parameter": parameter,
                "n_planets_total": len(group),
                "n_planets_direct_qc_excluded": int(direct_qc.sum()),
                "n_planets": len(finite),
                "n_systems": n_systems,
                "sample_adequacy": sample_adequacy,
                "median_delta": finite["delta"].median(),
                "median_delta_ci_low": _finite_quantile(signed_bootstrap, 0.025),
                "median_delta_ci_high": _finite_quantile(signed_bootstrap, 0.975),
                "median_abs_delta": finite["abs_delta"].median(),
                "median_abs_delta_ci_low": _finite_quantile(absolute_bootstrap, 0.025),
                "median_abs_delta_ci_high": _finite_quantile(absolute_bootstrap, 0.975),
                "p95_abs_delta": finite["abs_delta"].quantile(0.95),
                "repeatability_abs_delta_p95": threshold,
                "n_exceeding_repeatability_p95": n_exceeding,
                "fraction_exceeding_repeatability_p95": n_exceeding / len(finite) if len(finite) else np.nan,
                "practical_evidence": evidence,
                "bootstrap_unit": "koi_target_system",
                "bootstrap_statistic": "median_paired_planet_delta",
                "bootstrap_replicates": n_bootstrap,
                "bootstrap_seed": seed,
            }
        )
    return pd.DataFrame(rows)


def write_report(path: Path, discovery: pd.DataFrame, summaries: pd.DataFrame, thresholds: pd.DataFrame) -> None:
    incomplete = int(discovery["status"].ne("present").sum())
    lines = [
        "# Factorial validation analysis",
        "",
        "All comparisons are paired within planet and restricted to the validation bundle's intended target set.",
        "Eccentricity summaries were regenerated with `extract_eccentricity_posteriors_direct.py`; zeta summaries",
        "were derived from those direct e-omega posterior grids.",
        "",
        "## Interpretation boundary",
        "",
        "`exceeds_repeatability` labels mean only that an observed absolute paired shift is larger than the",
        "95th percentile of the available repeat-run shifts for that parameter. The repeat sample is small and",
        "targeted, so these labels are exploratory robustness signals, not p-values, causal claims, or proof that",
        "an analysis choice changes the population-level scientific conclusion.",
        "Any summary with fewer than five target systems is labeled underpowered and must not be used for subgroup inference.",
        "",
        "## Coverage",
        "",
        f"- Analysis mode: {'partial snapshot' if incomplete else 'complete matrix'}",
        f"- Expected and discovered FITS: {int(discovery['status'].eq('present').sum())}/{len(discovery)}",
        f"- Paired comparison summaries: {len(summaries)}",
        "- Bootstrap unit: target system (all planets in a sampled system remain clustered)",
        "",
        "## Repeatability thresholds",
        "",
        thresholds.to_markdown(index=False),
        "",
        "## Arm summaries",
        "",
        summaries.to_markdown(index=False),
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def run_analysis(
    *,
    validation_root: Path,
    metadata_root: Path,
    sample_path: Path,
    inventory_path: Path,
    output_dir: Path,
    n_proposals: int,
    e_max: float,
    density_error_mode: str,
    period_tol: float,
    min_importance_ess: float,
    e_grid_size: int,
    omega_grid_size: int,
    n_bootstrap: int,
    seed: int,
    allow_incomplete: bool = False,
) -> dict[str, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    target_sets = read_target_sets(metadata_root)
    discovery = discover_validation_fits(validation_root, target_sets)
    discovery_path = output_dir / "factorial_validation_discovery.csv"
    discovery.to_csv(discovery_path, index=False)
    missing = discovery[discovery["status"].ne("present")]
    if not allow_incomplete and not missing.empty:
        examples = ", ".join(f"{row.arm}/{row.koi_target}" for row in missing.head(12).itertuples())
        raise ValidationAnalysisError(
            f"Validation FITS are incomplete ({len(missing)} missing of {len(discovery)} expected): {examples}. "
            f"Discovery audit: {discovery_path}"
        )

    present_discovery = discovery[discovery["status"].eq("present")].copy()
    if present_discovery.empty:
        raise ValidationAnalysisError(
            f"No validation FITS were discovered. Discovery audit: {discovery_path}"
        )

    planets = build_planet_sample(inventory_path, sample_path, target_sets["all"])
    e_grid = np.linspace(0.0, e_max, e_grid_size)
    omega_grid = np.linspace(0.0, 2.0 * np.pi, omega_grid_size, endpoint=False)
    exclusion_path = output_dir / "factorial_validation_direct_exclusions.csv"
    try:
        arm_planets, exclusions = extract_all_arms(
            present_discovery,
            planets,
            output_dir / "direct_posteriors",
            e_grid=e_grid,
            omega_grid=omega_grid,
            n_proposals=n_proposals,
            e_max=e_max,
            density_error_mode=density_error_mode,
            period_tol=period_tol,
            min_importance_ess=min_importance_ess,
        )
    except DirectExtractionError as exc:
        exc.exclusions.to_csv(exclusion_path, index=False)
        raise
    except ValidationAnalysisError:
        pd.DataFrame(columns=["arm", "koi_target", "kepoi_name", "stage", "reason"]).to_csv(
            exclusion_path, index=False
        )
        raise
    exclusions.to_csv(exclusion_path, index=False)

    paired_planets, paired_metrics = build_paired_outputs(
        arm_planets,
        target_sets,
        allow_incomplete=allow_incomplete,
    )
    paired_metrics, thresholds = attach_repeatability_evidence(paired_metrics)
    summaries = summarize_arms(paired_metrics, n_bootstrap=n_bootstrap, seed=seed)

    paths = {
        "discovery": discovery_path,
        "direct_exclusions": exclusion_path,
        "arm_planets": output_dir / "factorial_validation_arm_planets.csv",
        "paired_planets": output_dir / "factorial_validation_paired_planets.csv",
        "paired_metrics": output_dir / "factorial_validation_paired_metrics.csv",
        "repeatability": output_dir / "factorial_validation_repeatability_thresholds.csv",
        "arm_summaries": output_dir / "factorial_validation_arm_summaries.csv",
        "report": output_dir / "factorial_validation_report.md",
    }
    arm_planets.to_csv(paths["arm_planets"], index=False)
    paired_planets.to_csv(paths["paired_planets"], index=False)
    paired_metrics.to_csv(paths["paired_metrics"], index=False)
    thresholds.to_csv(paths["repeatability"], index=False)
    summaries.to_csv(paths["arm_summaries"], index=False)
    write_report(paths["report"], discovery, summaries, thresholds)
    return paths


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Paired analysis of returned cloud factorial-validation FITS.",
        epilog=(
            "This command only reads returned FITS and writes local analysis products; it never launches cloud jobs. "
            "Practical-evidence labels are calibrated to observed repeat runs and are exploratory."
        ),
    )
    parser.add_argument(
        "--validation-root",
        default="cloud_ld_validation_batch",
        help="Unpacked validation result root; FITS are discovered recursively below it.",
    )
    parser.add_argument(
        "--metadata-root",
        default=None,
        help="Directory containing the three target manifests (default: validation root).",
    )
    parser.add_argument(
        "--inventory",
        default=None,
        help="Full-system inventory CSV (default: <metadata-root>/full_system_inventory.csv).",
    )
    parser.add_argument("--sample", default="outputs/canonical_sample_old_astropy_rawcc.csv")
    parser.add_argument("--config", default="config.json", help="Config used for default direct grid sizes.")
    parser.add_argument("--output-dir", default="outputs/factorial_validation")
    parser.add_argument("--n-proposals", type=int, default=150_000)
    parser.add_argument("--e-max", type=float, default=0.95)
    parser.add_argument(
        "--density-error-mode",
        choices=["symmetric-average", "split"],
        default="symmetric-average",
    )
    parser.add_argument("--period-tol", type=float, default=0.01)
    parser.add_argument(
        "--min-importance-ess",
        type=float,
        default=100.0,
        help="Direct-extractor ESS threshold retained as QC flags in paired outputs.",
    )
    parser.add_argument("--e-grid-size", type=int, default=None)
    parser.add_argument("--omega-grid-size", type=int, default=None)
    parser.add_argument("--bootstrap-replicates", type=int, default=10_000)
    parser.add_argument("--seed", type=int, default=20260710)
    parser.add_argument(
        "--allow-incomplete",
        action="store_true",
        help=(
            "Analyze only discovered FITS and complete within-planet arm pairs. "
            "The default remains strict and requires the full 82-fit matrix."
        ),
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.n_proposals <= 0:
        parser.error("--n-proposals must be positive")
    if not 0.0 < args.e_max < 1.0:
        parser.error("--e-max must be strictly between 0 and 1")
    if args.period_tol <= 0:
        parser.error("--period-tol must be positive")
    if args.min_importance_ess <= 0:
        parser.error("--min-importance-ess must be positive")
    if args.bootstrap_replicates <= 0:
        parser.error("--bootstrap-replicates must be positive")

    validation_root = resolve_path(args.validation_root)
    metadata_root = resolve_path(args.metadata_root) if args.metadata_root else validation_root
    inventory_path = resolve_path(args.inventory) if args.inventory else metadata_root / "full_system_inventory.csv"
    cfg = load_config(resolve_path(args.config))
    e_grid_size = args.e_grid_size or int(cfg["alderaan"]["eccentricity_grid_size"])
    omega_grid_size = args.omega_grid_size or int(cfg["alderaan"]["omega_grid_size"])
    if e_grid_size < 2 or omega_grid_size < 2:
        parser.error("direct grid sizes must each be at least 2")

    try:
        paths = run_analysis(
            validation_root=validation_root,
            metadata_root=metadata_root,
            sample_path=resolve_path(args.sample),
            inventory_path=inventory_path,
            output_dir=resolve_path(args.output_dir),
            n_proposals=args.n_proposals,
            e_max=args.e_max,
            density_error_mode=args.density_error_mode,
            period_tol=args.period_tol,
            min_importance_ess=args.min_importance_ess,
            e_grid_size=e_grid_size,
            omega_grid_size=omega_grid_size,
            n_bootstrap=args.bootstrap_replicates,
            seed=args.seed,
            allow_incomplete=args.allow_incomplete,
        )
    except (ValidationAnalysisError, FileNotFoundError, pd.errors.ParserError) as exc:
        print(f"ERROR: factorial validation analysis failed: {exc}", file=sys.stderr)
        return 2

    print("Wrote factorial validation analysis:")
    for label, path in paths.items():
        print(f"  {label}: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
