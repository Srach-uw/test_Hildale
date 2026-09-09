from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits

from common import load_config, normalize_dynesty_weights, root_path
from extract_eccentricity_posteriors_direct import (
    macdougall_rho_star_samp,
    weighted_quantile,
)


def read_berger_radius_errors(cfg: dict) -> pd.DataFrame:
    path = root_path(cfg, "berger_table2")
    if path is None:
        raise FileNotFoundError("Berger table2 is not configured")
    rows = []
    with gzip.open(path, "rt") as handle:
        for line in handle:
            try:
                rows.append(
                    {
                        "kepid": int(line[0:8].strip()),
                        "berger_rad": float(line[91:98]),
                        "berger_rad_err_hi": float(line[99:106]),
                        "berger_rad_err_lo": abs(float(line[107:115])),
                    }
                )
            except (TypeError, ValueError):
                continue
    return pd.DataFrame(rows).drop_duplicates("kepid")


def weighted_fraction(mask: np.ndarray, weights: np.ndarray) -> float:
    mask = np.asarray(mask, dtype=bool)
    weights = np.asarray(weights, dtype=float)
    if mask.shape != weights.shape:
        raise ValueError("mask and weights must have the same shape")
    total = weights.sum()
    if not np.isfinite(total) or total <= 0:
        return np.nan
    return float(weights[mask].sum() / total)


def grazing_exceeds_limit(fraction: float, limit: float = 0.05) -> bool:
    return bool(np.isfinite(fraction) and fraction > limit)


def catalog_density_solar(planet: pd.Series) -> float:
    direct = float(planet.get("rho_true_solar", np.nan))
    if np.isfinite(direct) and direct > 0:
        return direct
    rho_log = float(planet.get("rho_log", np.nan))
    if np.isfinite(rho_log):
        return 10.0**rho_log
    return np.nan


def radius_fractional_uncertainty(
    ror: np.ndarray,
    weights: np.ndarray,
    stellar_radius: float,
    radius_err_hi: float,
    radius_err_lo: float,
) -> tuple[float, float, float, float]:
    ror16, ror50, ror84 = weighted_quantile(ror, weights, [0.16, 0.5, 0.84])
    ror_sigma = 0.5 * (ror84 - ror16)
    stellar_sigma = 0.5 * (radius_err_hi + radius_err_lo)
    if ror50 <= 0 or stellar_radius <= 0:
        return ror16, ror50, ror84, np.nan
    frac = np.hypot(ror_sigma / ror50, stellar_sigma / stellar_radius)
    return ror16, ror50, ror84, float(frac)


def build_result_index(*roots: Path) -> dict[str, list[Path]]:
    index: dict[str, list[Path]] = {}
    seen: set[Path] = set()
    for root in roots:
        if not root.exists():
            continue
        for path in root.rglob("*-results.fits"):
            resolved = path.resolve()
            if resolved in seen or path.stat().st_size <= 0:
                continue
            seen.add(resolved)
            target = path.name.removesuffix("-results.fits")
            index.setdefault(target, []).append(resolved)
    return index


def result_file(
    target: str,
    result_index: dict[str, list[Path]],
) -> tuple[Path | None, bool, int]:
    matches = result_index.get(target, [])
    if not matches:
        return None, False, 0
    if len(matches) == 1:
        return matches[0], False, 1
    sizes = {path.stat().st_size for path in matches}
    # Prefer the largest product and record the decision in the audit ledger.
    selected = sorted(matches, key=lambda path: (-path.stat().st_size, str(path)))[0]
    return selected, len(sizes) > 1, len(matches)


def audit(
    summary: pd.DataFrame,
    archive_dir: Path,
    cloud_dir: Path,
    radius_errors: pd.DataFrame,
    *,
    grazing_fraction_limit: float = 0.05,
) -> pd.DataFrame:
    radius_lookup = radius_errors.set_index("kepid")
    result_index = build_result_index(archive_dir, cloud_dir)
    rows: list[dict] = []
    group_columns = ["koi_target"]
    if "transit_fit_source" in summary.columns:
        group_columns.append("transit_fit_source")
    grouped = summary.groupby(group_columns, sort=True)
    for group_key, planets in grouped:
        target = group_key[0] if isinstance(group_key, tuple) else group_key
        path, conflicting_duplicates, candidate_count = result_file(
            str(target), result_index
        )
        if path is None:
            for _, planet in planets.iterrows():
                rows.append(
                    {
                        "kepoi_name": planet["kepoi_name"],
                        "raw_result_file": "",
                        "raw_result_candidate_count": 0,
                        "gilbert_qc_available": False,
                        "gilbert_qc_reason": "raw_result_missing",
                    }
                )
            continue
        with fits.open(path, memmap=True) as hdul:
            samples = hdul["SAMPLES"].data
            names = set(samples.names)
            weights = normalize_dynesty_weights(np.asarray(samples["LN_WT"], dtype=float))
            for _, planet in planets.iterrows():
                idx = int(planet["alderaan_planet_index"])
                ror_col = f"ROR_{idx}"
                impact_col = f"IMPACT_{idx}"
                duration_col = f"DUR14_{idx}"
                base = {
                    "kepoi_name": planet["kepoi_name"],
                    "raw_result_file": str(path.resolve()),
                    "raw_result_candidate_count": candidate_count,
                    "raw_result_selection_rule": "largest_nonempty_then_path",
                    "raw_result_conflicting_duplicates": conflicting_duplicates,
                }
                if ror_col not in names or impact_col not in names or duration_col not in names:
                    rows.append(
                        {
                            **base,
                            "gilbert_qc_available": False,
                            "gilbert_qc_reason": "transit_columns_missing",
                        }
                    )
                    continue
                ror = np.asarray(samples[ror_col], dtype=float)
                impact = np.asarray(samples[impact_col], dtype=float)
                duration_days = np.asarray(samples[duration_col], dtype=float)
                valid = (
                    np.isfinite(ror)
                    & np.isfinite(impact)
                    & np.isfinite(duration_days)
                    & np.isfinite(weights)
                    & (weights >= 0)
                    & (ror > 0)
                    & (duration_days > 0)
                )
                if not valid.any() or weights[valid].sum() <= 0:
                    rows.append(
                        {
                            **base,
                            "gilbert_qc_available": False,
                            "gilbert_qc_reason": "no_valid_nested_samples",
                        }
                    )
                    continue
                w = weights[valid]
                w = w / w.sum()
                r = ror[valid]
                b = impact[valid]
                duration = duration_days[valid]
                grazing_fraction = weighted_fraction(b > (1.0 - r), w)
                rho_circular = macdougall_rho_star_samp(
                    float(planet["koi_period"]) * 86400.0,
                    duration * 86400.0,
                    r,
                    b,
                    np.zeros_like(r),
                    np.zeros_like(r),
                )
                rho_valid = np.isfinite(rho_circular) & (rho_circular > 0)
                if rho_valid.any() and w[rho_valid].sum() > 0:
                    rho_weights = w[rho_valid] / w[rho_valid].sum()
                    rho16, rho50, rho84 = weighted_quantile(
                        rho_circular[rho_valid],
                        rho_weights,
                        [0.16, 0.5, 0.84],
                    )
                    rho_catalog = catalog_density_solar(planet)
                    rho_ratio = rho50 / rho_catalog
                else:
                    rho16 = rho50 = rho84 = rho_ratio = np.nan
                if int(planet["kepid"]) in radius_lookup.index:
                    star = radius_lookup.loc[int(planet["kepid"])]
                    r16, r50, r84, rp_frac = radius_fractional_uncertainty(
                        r,
                        w,
                        float(star["berger_rad"]),
                        float(star["berger_rad_err_hi"]),
                        float(star["berger_rad_err_lo"]),
                    )
                else:
                    r16 = r50 = r84 = rp_frac = np.nan
                rows.append(
                    {
                        **base,
                        "gilbert_qc_available": True,
                        "gilbert_qc_reason": "",
                        "nested_grazing_fraction": grazing_fraction,
                        "ror16": r16,
                        "ror50": r50,
                        "ror84": r84,
                        "planet_radius_fractional_uncertainty_approx": rp_frac,
                        "rho_circular_p16_solar": rho16,
                        "rho_circular_median_solar": rho50,
                        "rho_circular_p84_solar": rho84,
                        "rho_catalog_solar": catalog_density_solar(planet),
                        "rho_circular_to_catalog_ratio": rho_ratio,
                        "gilbert_grazing_exclude": grazing_exceeds_limit(
                            grazing_fraction, grazing_fraction_limit
                        ),
                        "gilbert_grazing_fraction_limit": grazing_fraction_limit,
                        "gilbert_radius_precision_exclude": bool(rp_frac > 0.20),
                    }
                )
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Audit the published Gilbert et al. post-fit grazing and radius-precision cuts."
    )
    parser.add_argument(
        "--summary",
        default="outputs/eccentricity_posterior_summary_uniform_paired_full.csv",
    )
    parser.add_argument(
        "--archive-dir",
        default="inputs/alderaan_archive_raw_20260612/ALDERAAN_posteriors",
    )
    parser.add_argument(
        "--cloud-dir",
        default="alderaan_project/Results",
    )
    parser.add_argument(
        "--output",
        default="outputs/gilbert_postfit_qc_audit_uniform_paired_full.csv",
    )
    parser.add_argument(
        "--filtered-summary",
        default="outputs/eccentricity_posterior_summary_uniform_paired_gilbert_qc.csv",
        help="Summary with the published Gilbert grazing and radius-precision exclusions applied.",
    )
    parser.add_argument(
        "--grazing-fraction-limit",
        type=float,
        default=0.05,
        help=(
            "Maximum weighted grazing-sample fraction. Gilbert's analysis notebook "
            "sets 0.05; its command-line parser separately defaults to 0.01."
        ),
    )
    args = parser.parse_args()
    if not 0.0 <= args.grazing_fraction_limit <= 1.0:
        parser.error("--grazing-fraction-limit must be between 0 and 1")

    summary = pd.read_csv(args.summary)
    required = {
        "kepoi_name",
        "kepid",
        "koi_target",
        "alderaan_planet_index",
    }
    missing = required - set(summary.columns)
    if missing:
        raise ValueError(f"Summary is missing required columns: {sorted(missing)}")

    cfg = load_config()
    result = audit(
        summary,
        Path(args.archive_dir),
        Path(args.cloud_dir),
        read_berger_radius_errors(cfg),
        grazing_fraction_limit=args.grazing_fraction_limit,
    )
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(out, index=False)
    merged = summary.merge(result, on="kepoi_name", how="left", validate="one_to_one")
    available = merged["gilbert_qc_available"].fillna(False).astype(bool)
    grazing = merged["gilbert_grazing_exclude"].fillna(False).astype(bool)
    radius_precision = merged["gilbert_radius_precision_exclude"].fillna(False).astype(bool)
    merged["gilbert_postfit_qc_exclude"] = ~available | grazing | radius_precision
    merged["gilbert_postfit_qc_reasons"] = np.select(
        [
            ~available,
            grazing & radius_precision,
            grazing,
            radius_precision,
        ],
        [
            "qc_unavailable",
            "grazing_fraction_and_radius_precision",
            "grazing_fraction",
            "radius_precision",
        ],
        default="",
    )
    filtered = merged.loc[~merged["gilbert_postfit_qc_exclude"]].copy()
    filtered_out = Path(args.filtered_summary)
    filtered.to_csv(filtered_out, index=False)
    print(f"Wrote {len(result)} rows: {out}")
    print(f"Wrote {len(filtered)} Gilbert-QC rows: {filtered_out}")
    print(result["gilbert_qc_available"].value_counts(dropna=False).to_string())
    for col in ["gilbert_grazing_exclude", "gilbert_radius_precision_exclude"]:
        print(result[col].fillna(False).astype(bool).value_counts().to_string())


if __name__ == "__main__":
    main()
