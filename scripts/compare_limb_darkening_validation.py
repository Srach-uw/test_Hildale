"""Compare baseline ALDERAAN fits with limb-darkening-reference validation reruns.

Run this after the validation FITS are copied back into
alderaan_project/Results/sagear_ld_reference_validation.

It can be run before eccentricity extraction: in that case it reports direct
shape changes only. After running extract_eccentricity_posteriors.py for the
validation run, it also reports e/zeta changes.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits

from common import load_config, normalize_dynesty_weights, output_dir, root_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-run-id", default="sagear_missing")
    parser.add_argument("--validation-run-id", default="sagear_ld_reference_validation")
    parser.add_argument("--targets", default="outputs/ld_reference_validation_targets.csv")
    parser.add_argument("--sample", default="outputs/canonical_sample_old_astropy_rawcc.csv")
    parser.add_argument("--validation-catalog", default="cloud_ld_validation_batch/sagear_ld_reference_catalog.csv")
    parser.add_argument("--baseline-summary", default="outputs/eccentricity_posterior_summary_merged_paired_exact_qcprimary.csv")
    parser.add_argument("--validation-summary", default="outputs/eccentricity_posterior_summary_ld_reference_validation.csv")
    parser.add_argument("--out-tag", default="LD_REFERENCE_VALIDATION")
    return parser.parse_args()


def repo_root() -> Path:
    return Path(__file__).resolve().parent


def resolve(path: str | Path) -> Path:
    p = Path(path)
    return p if p.is_absolute() else repo_root() / p


def result_file(project: Path, run_id: str, target: str) -> Path:
    return project / "Results" / run_id / target / f"{target}-results.fits"


def weighted_quantile(values: np.ndarray, weights: np.ndarray, quantiles: list[float]) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)
    mask = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    if mask.sum() == 0:
        return np.full(len(quantiles), np.nan)
    values = values[mask]
    weights = weights[mask]
    order = np.argsort(values)
    values = values[order]
    weights = weights[order]
    cdf = np.cumsum(weights)
    cdf = cdf / cdf[-1]
    return np.interp(quantiles, cdf, values)


def shape_rows_for_run(project: Path, run_id: str, sample: pd.DataFrame, targets: list[str]) -> pd.DataFrame:
    rows = []
    for target in targets:
        path = result_file(project, run_id, target)
        target_sample = sample[sample["koi_target"].eq(target)].sort_values("koi_period").reset_index(drop=True)
        if not path.exists():
            for _, row in target_sample.iterrows():
                rows.append(
                    {
                        "run_id": run_id,
                        "koi_target": target,
                        "kepoi_name": row["kepoi_name"],
                        "kepid": row["kepid"],
                        "status": "missing_results_fits",
                    }
                )
            continue
        with fits.open(path, memmap=False) as hdul:
            samples = pd.DataFrame(np.array(hdul["SAMPLES"].data).byteswap().newbyteorder())
        weights = normalize_dynesty_weights(np.asarray(samples["LN_WT"], dtype=float))
        n_planets = sum(col.startswith("DUR14") for col in samples.columns)
        for idx in range(min(n_planets, len(target_sample))):
            suffix = f"_{idx}"
            if f"DUR14{suffix}" not in samples.columns:
                suffix = "" if idx == 0 and "DUR14" in samples.columns else suffix
            row = target_sample.iloc[idx]
            dur = np.asarray(samples[f"DUR14{suffix}"], dtype=float) * 24.0
            ror = np.asarray(samples[f"ROR{suffix}"], dtype=float)
            impact = np.asarray(samples[f"IMPACT{suffix}"], dtype=float)
            dur_q = weighted_quantile(dur, weights, [0.16, 0.5, 0.84])
            ror_q = weighted_quantile(ror, weights, [0.16, 0.5, 0.84])
            impact_q = weighted_quantile(impact, weights, [0.16, 0.5, 0.84])
            rows.append(
                {
                    "run_id": run_id,
                    "koi_target": target,
                    "kepoi_name": row["kepoi_name"],
                    "kepid": row["kepid"],
                    "disk": row["disk"],
                    "system": row["system"],
                    "koi_period": row["koi_period"],
                    "status": "matched",
                    "fit_planet_index": idx,
                    "shape_samples": len(samples),
                    "dur14_hr_p16": dur_q[0],
                    "dur14_hr_p50": dur_q[1],
                    "dur14_hr_p84": dur_q[2],
                    "ror_p16": ror_q[0],
                    "ror_p50": ror_q[1],
                    "ror_p84": ror_q[2],
                    "impact_p16": impact_q[0],
                    "impact_p50": impact_q[1],
                    "impact_p84": impact_q[2],
                }
            )
    return pd.DataFrame(rows)


def compare_shapes(shape: pd.DataFrame) -> pd.DataFrame:
    base = shape[shape["run_id"].eq("sagear_missing")].copy()
    val = shape[~shape["run_id"].eq("sagear_missing")].copy()
    if base.empty or val.empty:
        return pd.DataFrame()
    keep = [
        "koi_target",
        "kepoi_name",
        "status",
        "dur14_hr_p50",
        "ror_p50",
        "impact_p50",
        "shape_samples",
    ]
    out = base[keep].merge(val[keep], on=["koi_target", "kepoi_name"], suffixes=("_baseline", "_validation"), how="outer")
    out["delta_dur14_hr"] = out["dur14_hr_p50_validation"] - out["dur14_hr_p50_baseline"]
    out["frac_delta_dur14"] = out["delta_dur14_hr"] / out["dur14_hr_p50_baseline"]
    out["delta_ror"] = out["ror_p50_validation"] - out["ror_p50_baseline"]
    out["frac_delta_ror"] = out["delta_ror"] / out["ror_p50_baseline"]
    out["delta_impact"] = out["impact_p50_validation"] - out["impact_p50_baseline"]
    return out


def attach_eccentricity_comparison(comp: pd.DataFrame, baseline_path: Path, validation_path: Path) -> pd.DataFrame:
    if comp.empty or not baseline_path.exists() or not validation_path.exists():
        return comp
    base = pd.read_csv(baseline_path)
    val = pd.read_csv(validation_path)
    keep = ["kepoi_name", "e16", "e50", "e84", "zeta_median", "zeta_p16", "zeta_p84"]
    out = comp.merge(base[keep], on="kepoi_name", how="left", suffixes=("", "_baseline_e"))
    out = out.rename(columns={c: f"{c}_baseline" for c in keep if c != "kepoi_name"})
    out = out.merge(val[keep], on="kepoi_name", how="left")
    out = out.rename(columns={c: f"{c}_validation" for c in keep if c != "kepoi_name"})
    out["delta_e50"] = out["e50_validation"] - out["e50_baseline"]
    out["delta_zeta_median"] = out["zeta_median_validation"] - out["zeta_median_baseline"]
    return out


def summarize(comp: pd.DataFrame) -> pd.DataFrame:
    if comp.empty:
        return pd.DataFrame()
    rows = []
    comp = comp.copy()
    comp["population"] = comp.get("population", np.nan)
    for label, grp in comp.groupby("population", dropna=False):
        rows.append(
            {
                "population": label,
                "n": len(grp),
                "n_validation_shape_present": int(grp["status_validation"].eq("matched").sum())
                if "status_validation" in grp
                else 0,
                "median_frac_delta_dur14": grp["frac_delta_dur14"].median(),
                "median_delta_impact": grp["delta_impact"].median(),
                "median_delta_ror": grp["delta_ror"].median(),
                "median_delta_e50": grp["delta_e50"].median() if "delta_e50" in grp else np.nan,
                "median_delta_zeta": grp["delta_zeta_median"].median() if "delta_zeta_median" in grp else np.nan,
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    args = parse_args()
    cfg = load_config(repo_root() / "config.json")
    project = root_path(cfg, "alderaan_project")
    if project is None:
        raise ValueError("Missing alderaan_project in config")

    targets = pd.read_csv(resolve(args.targets))
    target_names = sorted(targets["koi_target"].dropna().unique())
    sample = pd.read_csv(resolve(args.sample))
    sample = sample[sample["koi_target"].isin(target_names)].copy()
    validation_catalog = pd.read_csv(resolve(args.validation_catalog), index_col=0)

    baseline_shape = shape_rows_for_run(project, args.baseline_run_id, sample, target_names)
    validation_shape = shape_rows_for_run(project, args.validation_run_id, sample, target_names)
    shape = pd.concat([baseline_shape, validation_shape], ignore_index=True)

    comp = compare_shapes(shape)
    if not comp.empty:
        context = targets[["koi_target", "selection_reason", "validation_priority"]].drop_duplicates("koi_target")
        pop = sample[["koi_target", "kepoi_name", "disk", "system"]].copy()
        pop["population"] = pop["disk"].astype(str) + "_" + pop["system"].astype(str) + "s"
        ld = validation_catalog.groupby("koi_id")[["limbdark_1", "limbdark_2"]].first().reset_index()
        ld = ld.rename(columns={"koi_id": "koi_target", "limbdark_1": "validation_limbdark_1", "limbdark_2": "validation_limbdark_2"})
        comp = comp.merge(context, on="koi_target", how="left")
        comp = comp.merge(pop[["koi_target", "kepoi_name", "population"]], on=["koi_target", "kepoi_name"], how="left")
        comp = comp.merge(ld, on="koi_target", how="left")
        comp = attach_eccentricity_comparison(comp, resolve(args.baseline_summary), resolve(args.validation_summary))

    out = output_dir()
    shape_path = out / f"ld_reference_validation_shape_rows_{args.out_tag}.csv"
    comp_path = out / f"ld_reference_validation_comparison_{args.out_tag}.csv"
    summary_path = out / f"ld_reference_validation_summary_{args.out_tag}.csv"
    shape.to_csv(shape_path, index=False)
    comp.to_csv(comp_path, index=False)
    summary = summarize(comp)
    summary.to_csv(summary_path, index=False)

    print("Wrote:")
    for path in [shape_path, comp_path, summary_path]:
        print(f"  {path}")
    print()
    print(f"Baseline shape rows: {len(baseline_shape)}")
    print(f"Validation shape rows: {len(validation_shape)}")
    if validation_shape["status"].eq("missing_results_fits").all():
        print()
        print("Validation FITS are not present yet. After the cloud run returns, run:")
        print(
            "  python extract_eccentricity_posteriors.py "
            "--sample outputs/canonical_sample_old_astropy_rawcc.csv "
            f"--run-id {args.validation_run_id} "
            "--impact-mode alderaan "
            "--posterior-subdir eccentricity_posteriors_ld_reference_validation "
            "--summary-out outputs/eccentricity_posterior_summary_ld_reference_validation.csv "
            "--coverage-out outputs/eccentricity_posterior_coverage_ld_reference_validation.csv "
            "--excluded-out outputs/eccentricity_posterior_excluded_ld_reference_validation.csv"
        )
        print("Then rerun this comparison script.")
    elif not summary.empty:
        print()
        print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
