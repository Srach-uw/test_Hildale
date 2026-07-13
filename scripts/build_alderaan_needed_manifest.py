from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from alderaan_batch import build_alderaan_catalog
from common import load_config, output_dir, root_path


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build final pre-ALDERAAN planet and system manifests using the best current sample labels."
    )
    parser.add_argument("--config", default=None)
    parser.add_argument("--sample", default=None)
    parser.add_argument("--summary", default=None)
    parser.add_argument("--coverage", default=None)
    parser.add_argument("--shape", default=None)
    parser.add_argument("--copy-to-external-output", default=None)
    args = parser.parse_args()

    cfg = load_config(args.config)
    out = output_dir()
    sample_path = Path(args.sample) if args.sample else out / "canonical_sample_old_astropy_rawcc.csv"
    summary_path = Path(args.summary) if args.summary else out / "eccentricity_posterior_summary_old_astropy_rawcc.csv"
    coverage_path = Path(args.coverage) if args.coverage else out / "alderaan_existing_coverage_detail.csv"
    shape_path = Path(args.shape) if args.shape else out / "alderaan_shape_diagnostics.csv"

    sample = pd.read_csv(sample_path)
    summary = pd.read_csv(summary_path)
    coverage = pd.read_csv(coverage_path)
    shape = pd.read_csv(shape_path)

    all_status = build_all_planet_status(sample, summary, coverage, shape)
    needed_planets = all_status[all_status["needs_alderaan"]].copy()
    runnable_needed_planets = needed_planets[needed_planets["can_seed_alderaan"]].copy()
    unseeded_needed_planets = needed_planets[~needed_planets["can_seed_alderaan"]].copy()
    all_needed_targets = build_target_manifest(sample, needed_planets)
    needed_targets = all_needed_targets[all_needed_targets["runnable_needed_planets"] > 0].copy()
    expanded_catalog_rows = sample[sample["koi_target"].isin(needed_targets["koi_target"])].copy()
    catalog = build_alderaan_catalog(expanded_catalog_rows, cfg)
    summary_table = build_summary(
        all_status,
        needed_planets,
        runnable_needed_planets,
        unseeded_needed_planets,
        all_needed_targets,
        needed_targets,
        catalog,
    )

    all_status_path = out / "alderaan_all_planets_status_best.csv"
    planet_path = out / "alderaan_needed_planets_best.csv"
    runnable_planet_path = out / "alderaan_runnable_needed_planets_best.csv"
    unseeded_planet_path = out / "alderaan_unseeded_needed_planets_best.csv"
    all_target_path = out / "alderaan_needed_targets_all_best.csv"
    target_path = out / "alderaan_needed_targets_best.csv"
    catalog_rows_path = out / "alderaan_needed_catalog_rows_best.csv"
    catalog_path = out / "alderaan_needed_catalog_best.csv"
    summary_path_out = out / "alderaan_needed_summary_best.csv"
    md_path = out / "alderaan_needed_manifest.md"

    all_status.to_csv(all_status_path, index=False)
    needed_planets.to_csv(planet_path, index=False)
    runnable_needed_planets.to_csv(runnable_planet_path, index=False)
    unseeded_needed_planets.to_csv(unseeded_planet_path, index=False)
    all_needed_targets.to_csv(all_target_path, index=False)
    needed_targets.to_csv(target_path, index=False)
    expanded_catalog_rows.to_csv(catalog_rows_path, index=False)
    catalog.to_csv(catalog_path)
    summary_table.to_csv(summary_path_out, index=False)
    write_markdown(
        md_path,
        summary_table,
        needed_planets,
        runnable_needed_planets,
        unseeded_needed_planets,
        all_needed_targets,
        needed_targets,
        catalog_path,
    )

    project = root_path(cfg, "alderaan_project")
    if project is not None:
        for rel in ["Catalogs", "Scripts", "Data", "Results", "Figures"]:
            (project / rel).mkdir(parents=True, exist_ok=True)
        catalog.to_csv(project / "Catalogs" / "sagear_needed_catalog.csv")
        needed_targets.to_csv(project / "sagear_needed_targets.csv", index=False)
        write_needed_scripts(cfg, needed_targets, "sagear_needed_catalog.csv", "sagear_needed_targets.csv")

    copy_root = Path(args.copy_to_external_output) if args.copy_to_external_output else None
    if copy_root:
        copy_root.mkdir(parents=True, exist_ok=True)
        for path in [
            all_status_path,
            planet_path,
            runnable_planet_path,
            unseeded_planet_path,
            all_target_path,
            target_path,
            catalog_rows_path,
            catalog_path,
            summary_path_out,
            md_path,
        ]:
            (copy_root / path.name).write_bytes(path.read_bytes())

    print("=== Needed ALDERAAN summary ===")
    print(summary_table.to_string(index=False))
    print(f"\nWrote: {all_status_path}")
    print(f"Wrote: {planet_path}")
    print(f"Wrote: {runnable_planet_path}")
    print(f"Wrote: {unseeded_planet_path}")
    print(f"Wrote: {all_target_path}")
    print(f"Wrote: {target_path}")
    print(f"Wrote: {catalog_rows_path}")
    print(f"Wrote: {catalog_path}")
    print(f"Wrote: {summary_path_out}")
    print(f"Wrote: {md_path}")
    if project is not None:
        print(f"Wrote project catalog: {project / 'Catalogs' / 'sagear_needed_catalog.csv'}")
        print(f"Wrote project targets: {project / 'sagear_needed_targets.csv'}")


def build_all_planet_status(
    sample: pd.DataFrame,
    summary: pd.DataFrame,
    coverage: pd.DataFrame,
    shape: pd.DataFrame,
) -> pd.DataFrame:
    cols = [
        "kepid",
        "kepoi_name",
        "koi_target",
        "disk",
        "system",
        "system_after_all_cuts",
        "system_definition",
        "P_thick",
        "koi_period",
        "koi_time0bk",
        "koi_depth",
        "koi_duration",
        "koi_impact",
        "koi_ror",
        "koi_model_snr",
        "koi_prad",
        "koi_count",
        "berger_teff",
        "rho_log",
        "berger2018_bin",
    ]
    df = sample[[c for c in cols if c in sample.columns]].copy()

    posterior_cols = [
        "kepoi_name",
        "e16",
        "e50",
        "e84",
        "zeta_median",
        "zeta_p16",
        "zeta_p84",
        "posterior_file",
    ]
    df = df.merge(summary[[c for c in posterior_cols if c in summary.columns]], on="kepoi_name", how="left")
    df["has_extracted_posterior_best"] = df["posterior_file"].notna()

    coverage_keep = [
        "kepoi_name",
        "has_alderaan_results_file",
        "has_extracted_posterior",
        "coverage_state",
    ]
    df = df.merge(
        coverage[[c for c in coverage_keep if c in coverage.columns]].drop_duplicates("kepoi_name"),
        on="kepoi_name",
        how="left",
    )
    missing_coverage_state = df["coverage_state"].isna()
    df.loc[missing_coverage_state, "coverage_state"] = np.where(
        df.loc[missing_coverage_state, "has_extracted_posterior_best"],
        "posterior_extracted",
        "unknown_missing_in_coverage_table",
    )
    df["has_alderaan_results_file"] = df["has_alderaan_results_file"].fillna(False)
    df["has_extracted_posterior"] = df["has_extracted_posterior_best"]

    shape_cols = [
        "kepoi_name",
        "alderaan_to_koi_duration_ratio",
        "alderaan_to_koi_ror_ratio",
        "alderaan_impact_med",
        "alderaan_impact_p84",
        "rho_frac_err",
    ]
    df = df.merge(shape[[c for c in shape_cols if c in shape.columns]], on="kepoi_name", how="left")

    duration_ratio = numeric(df, "alderaan_to_koi_duration_ratio")
    ror_ratio = numeric(df, "alderaan_to_koi_ror_ratio")
    impact_med = numeric(df, "alderaan_impact_med")
    zeta = numeric(df, "zeta_median")
    e50 = numeric(df, "e50")
    snr = numeric(df, "koi_model_snr")
    seed_period = numeric(df, "koi_period")
    seed_epoch = numeric(df, "koi_time0bk")
    seed_depth = numeric(df, "koi_depth")
    seed_duration = numeric(df, "koi_duration")

    df["flag_missing_posterior"] = ~df["has_extracted_posterior_best"]
    df["flag_no_results_file"] = df["coverage_state"].eq("no_results_file")
    df["flag_results_file_bad_match"] = df["coverage_state"].eq("results_file_no_match_or_bad_zeta")
    df["flag_missing_seed_period"] = ~np.isfinite(seed_period) | (seed_period <= 0)
    df["flag_missing_seed_epoch"] = ~np.isfinite(seed_epoch)
    df["flag_missing_seed_depth"] = ~np.isfinite(seed_depth) | (seed_depth <= 0)
    df["flag_missing_seed_duration"] = ~np.isfinite(seed_duration) | (seed_duration <= 0)
    df["can_seed_alderaan"] = ~(
        df["flag_missing_seed_period"]
        | df["flag_missing_seed_epoch"]
        | df["flag_missing_seed_depth"]
        | df["flag_missing_seed_duration"]
    )
    df["flag_duration_shift_25pct"] = np.abs(np.log(duration_ratio)) > np.log(1.25)
    df["flag_duration_shift_50pct"] = np.abs(np.log(duration_ratio)) > np.log(1.50)
    df["flag_ror_shift_25pct"] = np.abs(np.log(ror_ratio)) > np.log(1.25)
    df["flag_ror_shift_50pct"] = np.abs(np.log(ror_ratio)) > np.log(1.50)
    df["flag_impact_gt_0p85"] = impact_med > 0.85
    df["flag_impact_gt_1"] = impact_med > 1.0
    df["flag_zeta_lt_0p7"] = zeta < 0.7
    df["flag_zeta_gt_1p3"] = zeta > 1.3
    df["flag_e50_gt_0p5"] = e50 > 0.5
    df["flag_low_snr_lt_20"] = snr < 20

    df["hard_shape_problem"] = (
        df["flag_duration_shift_25pct"]
        | df["flag_ror_shift_25pct"]
        | df["flag_impact_gt_1"]
        | df["flag_results_file_bad_match"]
    )
    df["high_e_shape_tail"] = df["flag_e50_gt_0p5"] & (df["flag_zeta_lt_0p7"] | df["flag_zeta_gt_1p3"])
    df["low_snr_high_e_tail"] = df["flag_low_snr_lt_20"] & (df["flag_e50_gt_0p5"] | df["flag_zeta_lt_0p7"])

    df["priority_tier"] = np.select(
        [
            df["flag_missing_posterior"],
            df["hard_shape_problem"],
            df["high_e_shape_tail"],
            df["low_snr_high_e_tail"],
        ],
        [1, 2, 3, 4],
        default=0,
    )
    df["needs_alderaan"] = df["priority_tier"] > 0
    df["recommended_action"] = np.select(
        [
            df["flag_missing_posterior"] & df["flag_no_results_file"],
            df["flag_missing_posterior"] & df["flag_results_file_bad_match"],
            df["flag_missing_posterior"],
            df["hard_shape_problem"],
            df["high_e_shape_tail"],
            df["low_snr_high_e_tail"],
        ],
        [
            "run_missing_no_results_file",
            "rerun_results_file_bad_match_or_bad_zeta",
            "run_missing_unknown_status",
            "refit_hard_shape_pathology",
            "refit_high_e_zeta_tail",
            "refit_low_snr_high_e_tail",
        ],
        default="no_alderaan_needed_preliminary",
    )
    df["run_reason_detail"] = df.apply(reason_detail, axis=1)
    df["seed_blocker_detail"] = df.apply(seed_blocker_detail, axis=1)

    sort_cols = ["priority_tier", "disk", "system", "koi_target", "koi_period"]
    return df.sort_values(sort_cols, ascending=[True, True, True, True, True]).reset_index(drop=True)


def numeric(df: pd.DataFrame, col: str) -> pd.Series:
    if col not in df.columns:
        return pd.Series(np.nan, index=df.index)
    return pd.to_numeric(df[col], errors="coerce")


def reason_detail(row: pd.Series) -> str:
    reasons = []
    for flag, label in [
        ("flag_no_results_file", "no_results_file"),
        ("flag_results_file_bad_match", "results_file_no_match_or_bad_zeta"),
        ("flag_duration_shift_25pct", "duration_shift_gt_25pct"),
        ("flag_ror_shift_25pct", "ror_shift_gt_25pct"),
        ("flag_impact_gt_1", "impact_gt_1"),
        ("flag_impact_gt_0p85", "impact_gt_0p85"),
        ("flag_zeta_lt_0p7", "zeta_lt_0p7"),
        ("flag_zeta_gt_1p3", "zeta_gt_1p3"),
        ("flag_e50_gt_0p5", "e50_gt_0p5"),
        ("flag_low_snr_lt_20", "snr_lt_20"),
    ]:
        if bool(row.get(flag, False)):
            reasons.append(label)
    return ";".join(reasons) if reasons else "none"


def seed_blocker_detail(row: pd.Series) -> str:
    blockers = []
    for flag, label in [
        ("flag_missing_seed_period", "missing_or_invalid_period"),
        ("flag_missing_seed_epoch", "missing_epoch"),
        ("flag_missing_seed_depth", "missing_or_invalid_depth"),
        ("flag_missing_seed_duration", "missing_or_invalid_duration"),
    ]:
        if bool(row.get(flag, False)):
            blockers.append(label)
    return ";".join(blockers) if blockers else "seed_ready"


def build_target_manifest(sample: pd.DataFrame, needed_planets: pd.DataFrame) -> pd.DataFrame:
    planet_reasons = (
        needed_planets.groupby("koi_target")
        .agg(
            needed_planets=("kepoi_name", "count"),
            runnable_needed_planets=("can_seed_alderaan", "sum"),
            unseeded_needed_planets=("can_seed_alderaan", lambda x: int((~x.astype(bool)).sum())),
            needed_kepoi_names=("kepoi_name", lambda x: ",".join(map(str, x))),
            min_priority_tier=("priority_tier", "min"),
            recommended_actions=("recommended_action", lambda x: ";".join(sorted(set(map(str, x))))),
            run_reason_detail=("run_reason_detail", lambda x: "|".join(sorted(set(map(str, x))))),
            seed_blocker_detail=("seed_blocker_detail", lambda x: "|".join(sorted(set(map(str, x))))),
        )
        .reset_index()
    )
    runnable_names = (
        needed_planets[needed_planets["can_seed_alderaan"]]
        .groupby("koi_target")["kepoi_name"]
        .agg(lambda x: ",".join(map(str, x)))
        .rename("runnable_kepoi_names")
        .reset_index()
    )
    unseeded_names = (
        needed_planets[~needed_planets["can_seed_alderaan"]]
        .groupby("koi_target")["kepoi_name"]
        .agg(lambda x: ",".join(map(str, x)))
        .rename("unseeded_kepoi_names")
        .reset_index()
    )
    planet_reasons = planet_reasons.merge(runnable_names, on="koi_target", how="left")
    planet_reasons = planet_reasons.merge(unseeded_names, on="koi_target", how="left")
    planet_reasons["runnable_kepoi_names"] = planet_reasons["runnable_kepoi_names"].fillna("")
    planet_reasons["unseeded_kepoi_names"] = planet_reasons["unseeded_kepoi_names"].fillna("")
    target_meta = (
        sample.groupby("koi_target")
        .agg(
            kepid=("kepid", "first"),
            disk=("disk", lambda x: ",".join(sorted(set(map(str, x))))),
            system=("system", lambda x: ",".join(sorted(set(map(str, x))))),
            P_thick=("P_thick", "median"),
            target_sample_planets=("kepoi_name", "count"),
            all_kepoi_names=("kepoi_name", lambda x: ",".join(map(str, x))),
            min_period=("koi_period", "min"),
            max_period=("koi_period", "max"),
            max_snr=("koi_model_snr", "max"),
            max_prad=("koi_prad", "max"),
        )
        .reset_index()
    )
    out = planet_reasons.merge(target_meta, on="koi_target", how="left")
    out["target_run_type"] = np.where(
        out["runnable_needed_planets"].eq(0),
        "blocked_missing_transit_seed",
        np.where(
        out["min_priority_tier"].eq(1),
        "required_missing_or_bad_extraction",
        np.where(out["min_priority_tier"].eq(2), "required_refit_hard_shape", "recommended_refit_tail"),
        ),
    )
    return out.sort_values(["min_priority_tier", "disk", "system", "koi_target"]).reset_index(drop=True)


def build_summary(
    all_status: pd.DataFrame,
    needed_planets: pd.DataFrame,
    runnable_needed_planets: pd.DataFrame,
    unseeded_needed_planets: pd.DataFrame,
    all_needed_targets: pd.DataFrame,
    needed_targets: pd.DataFrame,
    catalog: pd.DataFrame,
) -> pd.DataFrame:
    rows = []
    rows.append(
        {
            "metric": "all_current_best_sample_planets",
            "count": int(len(all_status)),
            "note": "Rows in canonical_sample_old_astropy_rawcc.csv.",
        }
    )
    rows.append(
        {
            "metric": "existing_extracted_posteriors",
            "count": int(all_status["has_extracted_posterior_best"].sum()),
            "note": "Planets with extracted e,omega grids from local ALDERAAN FITS.",
        }
    )
    rows.append(
        {
            "metric": "needed_planet_rows_total",
            "count": int(len(needed_planets)),
            "note": "Tier 1 missing plus Tier 2/3/4 refit candidates.",
        }
    )
    rows.append(
        {
            "metric": "runnable_needed_planets",
            "count": int(len(runnable_needed_planets)),
            "note": "Needed planets with finite period, epoch, depth, and duration seeds for ALDERAAN.",
        }
    )
    rows.append(
        {
            "metric": "unseeded_needed_planets",
            "count": int(len(unseeded_needed_planets)),
            "note": "Needed planets lacking ALDERAAN launch seeds, usually missing KOI depth.",
        }
    )
    for tier, label in [
        (1, "tier1_missing_or_bad_extraction_planets"),
        (2, "tier2_hard_shape_refit_planets"),
        (3, "tier3_high_e_zeta_tail_refit_planets"),
        (4, "tier4_low_snr_tail_refit_planets"),
    ]:
        rows.append(
            {
                "metric": label,
                "count": int((needed_planets["priority_tier"] == tier).sum()),
                "note": "Planet-level count.",
            }
        )
    rows.append(
        {
            "metric": "all_needed_unique_koi_targets",
            "count": int(len(all_needed_targets)),
            "note": "Unique KOI systems with at least one planet needing ALDERAAN.",
        }
    )
    rows.append(
        {
            "metric": "runnable_unique_koi_targets",
            "count": int(len(needed_targets)),
            "note": "Unique KOI systems with at least one runnable needed planet.",
        }
    )
    rows.append(
        {
            "metric": "alderaan_catalog_rows_after_system_expansion",
            "count": int(len(catalog)),
            "note": "Runnable-system catalog rows after expanding to all seedable planets in those systems.",
        }
    )
    for (disk, system), sub in needed_planets.groupby(["disk", "system"]):
        rows.append(
            {
                "metric": f"needed_planets_{disk}_{system}",
                "count": int(len(sub)),
                "note": "Planet-level needed rows by current best labels.",
            }
        )
    return pd.DataFrame(rows)


def write_needed_scripts(cfg: dict, targets: pd.DataFrame, catalog_name: str, target_name: str) -> None:
    project = root_path(cfg, "alderaan_project")
    repo = root_path(cfg, "alderaan_repo")
    assert project is not None
    run_id = cfg["alderaan"].get("needed_run_id", "sagear_needed")
    mission = cfg["alderaan"]["mission"]
    cadence = cfg["alderaan"].get("cadence", "long")
    scripts = project / "Scripts"
    data_dir = project / "Data"
    logs_dir = project / "Logs"
    data_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)
    target_path = project / target_name

    download = scripts / "download_needed_lightcurves.ps1"
    lines = [
        "param(",
        f"    [string]$TargetCsv = '{target_path}',",
        "    [int]$Start = 0,",
        "    [int]$Limit = 0,",
        f"    [string]$Cadence = '{cadence}'",
        ")",
        "$ErrorActionPreference = 'Stop'",
        "$env:OMP_NUM_THREADS = '1'",
        "$env:MKL_NUM_THREADS = '1'",
        "$env:OPENBLAS_NUM_THREADS = '1'",
        "$env:NUMEXPR_NUM_THREADS = '1'",
        f"$DataDir = '{data_dir}'",
        "New-Item -ItemType Directory -Force -Path $DataDir | Out-Null",
    ]
    if repo and (repo / "bin" / "get_kepler_data.py").exists():
        lines += [
            f"$AlderaanGet = '{repo / 'bin' / 'get_kepler_data.py'}'",
            "$targets = Import-Csv $TargetCsv | Sort-Object koi_target",
            "if ($Limit -gt 0) {",
            "    $targets = $targets | Select-Object -Skip $Start -First $Limit",
            "} elseif ($Start -gt 0) {",
            "    $targets = $targets | Select-Object -Skip $Start",
            "}",
            "$kepids = $targets | Select-Object -ExpandProperty kepid -Unique",
            "Push-Location $DataDir",
            "try {",
            "    foreach ($kepid in $kepids) {",
            "        Write-Host \"Downloading KIC $kepid\"",
            "        python $AlderaanGet $kepid -c $Cadence -t lightcurve -o \"get_${kepid}_lc.sh\" --cmdtype wget",
            "        bash \"get_${kepid}_lc.sh\"",
            "    }",
            "} finally {",
            "    Pop-Location",
            "}",
        ]
    else:
        lines.append("# ALDERAAN repo not found; cannot generate data download commands.")
    download.write_text("\n".join(lines) + "\n", encoding="utf-8")

    run = scripts / "run_needed_alderaan.ps1"
    lines = [
        "param(",
        f"    [string]$TargetCsv = '{target_path}',",
        f"    [string]$CatalogName = '{catalog_name}',",
        f"    [string]$RunId = '{run_id}',",
        "    [int]$MaxPriority = 4,",
        "    [int]$Start = 0,",
        "    [int]$Limit = 0,",
        "    [ValidateSet('all','detrend','noise','fit')] [string]$Phase = 'all'",
        ")",
        "$ErrorActionPreference = 'Stop'",
        "# Run this from an activated ALDERAAN conda environment.",
        "$env:OMP_NUM_THREADS = '1'",
        "$env:MKL_NUM_THREADS = '1'",
        "$env:OPENBLAS_NUM_THREADS = '1'",
        "$env:NUMEXPR_NUM_THREADS = '1'",
        f"$ProjectDir = '{project}'",
        "$LogDir = Join-Path $ProjectDir \"Logs\"",
        "New-Item -ItemType Directory -Force -Path $LogDir | Out-Null",
    ]
    if repo and repo.exists():
        lines += [
            f"$AlderaanRepo = '{repo}'",
            "$targets = Import-Csv $TargetCsv | Where-Object { [int]$_.min_priority_tier -le $MaxPriority } | Sort-Object min_priority_tier,koi_target",
            "if ($Limit -gt 0) {",
            "    $targets = $targets | Select-Object -Skip $Start -First $Limit",
            "} elseif ($Start -gt 0) {",
            "    $targets = $targets | Select-Object -Skip $Start",
            "}",
            "$failed = @()",
            "Push-Location $AlderaanRepo",
            "try {",
            "    foreach ($row in $targets) {",
            "        $target = $row.koi_target",
            "        $log = Join-Path $LogDir \"$RunId-$target.log\"",
            "        Write-Host \"Running $target priority=$($row.min_priority_tier) action=$($row.target_run_type)\"",
            "        try {",
            "            if ($Phase -in @('all','detrend')) {",
            f"                python bin/detrend_and_estimate_ttvs.py --mission {mission} --target $target --run_id $RunId --project_dir $ProjectDir --data_dir \"$ProjectDir\\Data\\\" --catalog $CatalogName *>> $log",
            "            }",
            "            if ($Phase -in @('all','noise')) {",
            f"                python bin/analyze_autocorrelated_noise.py --mission {mission} --target $target --run_id $RunId --project_dir $ProjectDir --data_dir \"$ProjectDir\\Data\\\" --catalog $CatalogName *>> $log",
            "            }",
            "            if ($Phase -in @('all','fit')) {",
            f"                python bin/fit_transit_shape_simultaneous_nested.py --mission {mission} --target $target --run_id $RunId --project_dir $ProjectDir *>> $log",
            "            }",
            "        } catch {",
            "            $failed += $target",
            "            \"FAILED ${target}: $($_.Exception.Message)\" | Tee-Object -Append -FilePath (Join-Path $LogDir \"$RunId-failures.log\")",
            "        }",
            "    }",
            "} finally {",
            "    Pop-Location",
            "}",
            "if ($failed.Count -gt 0) {",
            "    throw \"ALDERAAN failed for $($failed.Count) target(s). See $LogDir.\"",
            "}",
        ]
    else:
        lines.append("# ALDERAAN repo not found; cannot generate run commands.")
    run.write_text("\n".join(lines) + "\n", encoding="utf-8")

    download_sh = scripts / "download_needed_lightcurves.sh"
    download_sh.write_text(
        "\n".join(
            [
                "#!/usr/bin/env bash",
                "set -euo pipefail",
                "export OMP_NUM_THREADS=1",
                "export MKL_NUM_THREADS=1",
                "export OPENBLAS_NUM_THREADS=1",
                "export NUMEXPR_NUM_THREADS=1",
                'SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"',
                'PROJECT_DIR="${PROJECT_DIR:-$(cd "$SCRIPT_DIR/.." && pwd)}"',
                'RESEARCH_ROOT="${RESEARCH_ROOT:-$(cd "$PROJECT_DIR/../.." && pwd)}"',
                'ALDERAAN_REPO="${ALDERAAN_REPO:-$RESEARCH_ROOT/external/alderaan}"',
                'TARGET_CSV="${1:-$PROJECT_DIR/' + target_name + '}"',
                'DATA_DIR="${DATA_DIR:-$PROJECT_DIR/Data}"',
                'CADENCE="${CADENCE:-' + cadence + '}"',
                'START="${START:-0}"',
                'LIMIT="${LIMIT:-0}"',
                'mkdir -p "$DATA_DIR"',
                'if [[ ! -f "$ALDERAAN_REPO/bin/get_kepler_data.py" ]]; then',
                '  echo "Cannot find ALDERAAN repo. Set ALDERAAN_REPO=/path/to/alderaan" >&2',
                "  exit 2",
                "fi",
                'TMP_KEPIDS="$(mktemp)"',
                "python - \"$TARGET_CSV\" \"$START\" \"$LIMIT\" > \"$TMP_KEPIDS\" <<'PY'",
                "import sys",
                "import pandas as pd",
                "target_csv, start, limit = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])",
                "df = pd.read_csv(target_csv).sort_values('koi_target')",
                "if start:",
                "    df = df.iloc[start:]",
                "if limit:",
                "    df = df.iloc[:limit]",
                "for kepid in df['kepid'].dropna().astype(int).drop_duplicates():",
                "    print(kepid)",
                "PY",
                'cd "$DATA_DIR"',
                'while read -r kepid; do',
                '  [[ -z "$kepid" ]] && continue',
                '  echo "Downloading KIC $kepid"',
                '  python "$ALDERAAN_REPO/bin/get_kepler_data.py" "$kepid" -c "$CADENCE" -t lightcurve -o "get_${kepid}_lc.sh" --cmdtype wget',
                '  bash "get_${kepid}_lc.sh"',
                'done < "$TMP_KEPIDS"',
                'rm -f "$TMP_KEPIDS"',
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    run_sh = scripts / "run_needed_alderaan_parallel.sh"
    run_sh.write_text(
        "\n".join(
            [
                "#!/usr/bin/env bash",
                "set -euo pipefail",
                "export OMP_NUM_THREADS=1",
                "export MKL_NUM_THREADS=1",
                "export OPENBLAS_NUM_THREADS=1",
                "export NUMEXPR_NUM_THREADS=1",
                'SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"',
                'PROJECT_DIR="${PROJECT_DIR:-$(cd "$SCRIPT_DIR/.." && pwd)}"',
                'RESEARCH_ROOT="${RESEARCH_ROOT:-$(cd "$PROJECT_DIR/../.." && pwd)}"',
                'ALDERAAN_REPO="${ALDERAAN_REPO:-$RESEARCH_ROOT/external/alderaan}"',
                'TARGET_CSV="${1:-$PROJECT_DIR/' + target_name + '}"',
                'CATALOG_NAME="${CATALOG_NAME:-' + catalog_name + '}"',
                'RUN_ID="${RUN_ID:-' + run_id + '}"',
                'MAX_PRIORITY="${MAX_PRIORITY:-4}"',
                'START="${START:-0}"',
                'LIMIT="${LIMIT:-0}"',
                'JOBS="${JOBS:-1}"',
                'PHASE="${PHASE:-all}"',
                f'MISSION="{mission}"',
                'LOG_DIR="$PROJECT_DIR/Logs/$RUN_ID"',
                'mkdir -p "$LOG_DIR"',
                'if [[ ! -d "$ALDERAAN_REPO" ]]; then',
                '  echo "Cannot find ALDERAAN repo. Set ALDERAAN_REPO=/path/to/alderaan" >&2',
                "  exit 2",
                "fi",
                'TMP_TARGETS="$(mktemp)"',
                "python - \"$TARGET_CSV\" \"$MAX_PRIORITY\" \"$START\" \"$LIMIT\" > \"$TMP_TARGETS\" <<'PY'",
                "import sys",
                "import pandas as pd",
                "target_csv, max_priority, start, limit = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])",
                "df = pd.read_csv(target_csv)",
                "df = df[pd.to_numeric(df['min_priority_tier'], errors='coerce') <= max_priority]",
                "df = df.sort_values(['min_priority_tier', 'koi_target'])",
                "if start:",
                "    df = df.iloc[start:]",
                "if limit:",
                "    df = df.iloc[:limit]",
                "for target in df['koi_target'].drop_duplicates():",
                "    print(target)",
                "PY",
                "run_one() {",
                '  local target="$1"',
                '  local log="$LOG_DIR/${target}.log"',
                '  echo "Running $target" | tee "$log"',
                '  cd "$ALDERAAN_REPO"',
                '  if [[ "$PHASE" == "all" || "$PHASE" == "detrend" ]]; then',
                '    python bin/detrend_and_estimate_ttvs.py --mission "$MISSION" --target "$target" --run_id "$RUN_ID" --project_dir "$PROJECT_DIR" --data_dir "$PROJECT_DIR/Data/" --catalog "$CATALOG_NAME" >> "$log" 2>&1',
                "  fi",
                '  if [[ "$PHASE" == "all" || "$PHASE" == "noise" ]]; then',
                '    python bin/analyze_autocorrelated_noise.py --mission "$MISSION" --target "$target" --run_id "$RUN_ID" --project_dir "$PROJECT_DIR" --data_dir "$PROJECT_DIR/Data/" --catalog "$CATALOG_NAME" >> "$log" 2>&1',
                "  fi",
                '  if [[ "$PHASE" == "all" || "$PHASE" == "fit" ]]; then',
                '    python bin/fit_transit_shape_simultaneous_nested.py --mission "$MISSION" --target "$target" --run_id "$RUN_ID" --project_dir "$PROJECT_DIR" >> "$log" 2>&1',
                "  fi",
                "}",
                "export -f run_one",
                "export ALDERAAN_REPO CATALOG_NAME LOG_DIR MISSION PHASE PROJECT_DIR RUN_ID",
                'xargs -a "$TMP_TARGETS" -n 1 -P "$JOBS" bash -c \'run_one "$@"\' _',
                'rm -f "$TMP_TARGETS"',
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def write_markdown(
    path: Path,
    summary_table: pd.DataFrame,
    needed_planets: pd.DataFrame,
    runnable_needed_planets: pd.DataFrame,
    unseeded_needed_planets: pd.DataFrame,
    all_needed_targets: pd.DataFrame,
    needed_targets: pd.DataFrame,
    catalog_path: Path,
) -> None:
    by_tier = (
        needed_planets.groupby(["priority_tier", "recommended_action"])
        .size()
        .rename("planets")
        .reset_index()
        .sort_values(["priority_tier", "recommended_action"])
    )
    by_pop = (
        needed_planets.groupby(["disk", "system", "priority_tier"])
        .size()
        .rename("planets")
        .reset_index()
        .sort_values(["disk", "system", "priority_tier"])
    )
    target_by_tier = (
        needed_targets.groupby(["min_priority_tier", "target_run_type"])
        .size()
        .rename("targets")
        .reset_index()
        .sort_values(["min_priority_tier"])
    )
    all_target_by_type = (
        all_needed_targets.groupby(["target_run_type"])
        .size()
        .rename("targets")
        .reset_index()
        .sort_values(["target_run_type"])
    )
    unseeded_by_pop = (
        unseeded_needed_planets.groupby(["disk", "system", "seed_blocker_detail"])
        .size()
        .rename("planets")
        .reset_index()
        .sort_values(["disk", "system", "seed_blocker_detail"])
    )
    lines = [
        "# ALDERAAN Needed Manifest",
        "",
        "This manifest uses the current best pre-ALDERAAN sample labels:",
        "",
        "- old-Astropy disk classifier;",
        "- raw confirmed/candidate KOI multiplicity for single/multi labels.",
        "",
        "Priority tiers:",
        "",
        "- Tier 1: missing extracted posterior or existing results file failed matching/zeta extraction.",
        "- Tier 2: existing posterior but hard shape pathology, such as >25% duration/RpRs shift or impact >1.",
        "- Tier 3: existing posterior in high-e and extreme-zeta tail.",
        "- Tier 4: lower-priority low-SNR high-e/zeta-tail stress tests.",
        "",
        "Launchability:",
        "",
        "- `alderaan_needed_planets_best.csv` is the full scientific need list.",
        "- `alderaan_runnable_needed_planets_best.csv` is the immediate ALDERAAN run list at planet level.",
        "- `alderaan_unseeded_needed_planets_best.csv` contains needed planets lacking launch seeds, usually missing KOI depth.",
        "- `alderaan_needed_targets_best.csv` is the runnable system-level execution list.",
        "",
        "## Summary",
        "",
        summary_table.to_markdown(index=False),
        "",
        "## Needed Planets By Tier",
        "",
        by_tier.to_markdown(index=False),
        "",
        "## Needed Planets By Population",
        "",
        by_pop.to_markdown(index=False),
        "",
        "## Needed KOI Targets By Tier",
        "",
        target_by_tier.to_markdown(index=False),
        "",
        "## All Needed Target Types",
        "",
        all_target_by_type.to_markdown(index=False),
        "",
        "## Unseeded Needed Planets",
        "",
        unseeded_by_pop.to_markdown(index=False) if len(unseeded_by_pop) else "None.",
        "",
        unseeded_needed_planets[
            [
                "kepoi_name",
                "kepid",
                "disk",
                "system",
                "koi_period",
                "koi_duration",
                "recommended_action",
                "seed_blocker_detail",
            ]
        ]
        .head(50)
        .to_markdown(index=False, floatfmt=".4f")
        if len(unseeded_needed_planets)
        else "",
        "",
        "## First 30 Target Systems",
        "",
        needed_targets.head(30).to_markdown(index=False, floatfmt=".4f"),
        "",
        f"ALDERAAN catalog written to `{catalog_path}`.",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    main()
