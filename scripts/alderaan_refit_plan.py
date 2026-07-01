from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from alderaan_batch import build_alderaan_catalog, write_download_script, write_run_script
from common import load_config, output_dir, root_path


def main() -> None:
    parser = argparse.ArgumentParser(description="Build ALDERAAN missing/refit/control manifests from current diagnostics.")
    parser.add_argument("--config", default=None)
    parser.add_argument("--sample", default=None)
    parser.add_argument("--coverage", default=None)
    parser.add_argument("--shape", default=None)
    parser.add_argument("--flagged", default=None)
    parser.add_argument("--n-suspicious", type=int, default=40)
    parser.add_argument("--n-controls-per-bin", type=int, default=5)
    parser.add_argument("--n-missing-validation-per-bin", type=int, default=8)
    args = parser.parse_args()

    cfg = load_config(args.config)
    out = output_dir()
    sample = pd.read_csv(Path(args.sample) if args.sample else out / "canonical_sample_diagnostic.csv")
    coverage = pd.read_csv(Path(args.coverage) if args.coverage else out / "alderaan_existing_coverage_detail.csv")
    shape = pd.read_csv(Path(args.shape) if args.shape else out / "alderaan_shape_diagnostics.csv")
    flagged = pd.read_csv(Path(args.flagged) if args.flagged else out / "alderaan_shape_diagnostics_flagged.csv")

    missing_systems = build_missing_manifest(coverage)
    suspicious_systems = build_suspicious_manifest(flagged, args.n_suspicious)
    controls = build_control_manifest(shape, args.n_controls_per_bin)
    missing_validation = build_missing_validation_manifest(coverage, sample, args.n_missing_validation_per_bin)

    full_manifest = pd.concat(
        [
            missing_systems.assign(batch="all_missing"),
            suspicious_systems.assign(batch="suspicious_existing_refit"),
            controls.assign(batch="clean_controls"),
            missing_validation.assign(batch="missing_validation"),
        ],
        ignore_index=True,
    )
    full_manifest = fill_target_metadata(full_manifest, sample)
    full_manifest = full_manifest.sort_values(["batch", "disk", "system", "koi_target"]).drop_duplicates(
        ["batch", "koi_target", "reason"]
    )

    full_path = out / "alderaan_refit_full_manifest.csv"
    full_manifest.to_csv(full_path, index=False)

    validation = pd.concat([suspicious_systems, controls, missing_validation], ignore_index=True)
    validation = fill_target_metadata(validation, sample)
    validation = validation.drop_duplicates("koi_target").sort_values(["priority", "disk", "system", "koi_target"])
    validation_path = out / "alderaan_refit_validation_targets.csv"
    validation.to_csv(validation_path, index=False)

    missing_path = out / "alderaan_missing_targets_all.csv"
    missing_systems.to_csv(missing_path, index=False)

    project = root_path(cfg, "alderaan_project")
    if project is not None:
        for rel in ["Catalogs", "Data", "Results", "Figures", "Scripts"]:
            (project / rel).mkdir(parents=True, exist_ok=True)
        selected = sample[sample["koi_target"].isin(validation["koi_target"])]
        catalog = build_alderaan_catalog(selected, cfg)
        catalog_path = project / "Catalogs" / "sagear_refit_validation_catalog.csv"
        catalog.to_csv(catalog_path)
        # alderaan_batch.write_run_script currently emits the historical
        # validation catalog name, so keep this alias in sync.
        catalog.to_csv(project / "Catalogs" / "sagear_validation_catalog.csv")
        project_targets = project / "sagear_refit_validation_targets.csv"
        validation.to_csv(project_targets, index=False)
        write_download_script(cfg, str(project_targets))
        write_run_script(cfg, str(project_targets))

    write_markdown(full_manifest, validation, missing_systems)

    print(f"Wrote: {full_path}")
    print(f"Wrote: {validation_path}")
    print(f"Wrote: {missing_path}")
    print("\nValidation target counts:")
    print(validation.groupby(["reason", "disk", "system"]).size().rename("systems").reset_index().to_string(index=False))


def fill_target_metadata(manifest: pd.DataFrame, sample: pd.DataFrame) -> pd.DataFrame:
    target_meta = (
        sample.groupby("koi_target")
        .agg(
            kepid_fill=("kepid", "first"),
            disk_fill=("disk", "first"),
            system_fill=("system", "first"),
            P_thick_fill=("P_thick", "median"),
            n_planets_fill=("kepoi_name", "count"),
            min_period_fill=("koi_period", "min"),
            max_period_fill=("koi_period", "max"),
        )
        .reset_index()
    )
    out = manifest.merge(target_meta, on="koi_target", how="left")
    for col in ["kepid", "disk", "system", "P_thick", "n_planets", "min_period", "max_period"]:
        fill = f"{col}_fill"
        if col not in out:
            out[col] = np.nan
        out[col] = out[col].where(out[col].notna(), out[fill])
        out = out.drop(columns=[fill])
    out["kepid"] = pd.to_numeric(out["kepid"], errors="coerce").astype("Int64")
    out["n_planets"] = pd.to_numeric(out["n_planets"], errors="coerce").astype("Int64")
    return out


def build_missing_manifest(coverage: pd.DataFrame) -> pd.DataFrame:
    missing = coverage[coverage["coverage_state"] == "no_results_file"].copy()
    rows = (
        missing.groupby("koi_target")
        .agg(
            kepid=("kepid", "first"),
            disk=("disk", "first"),
            system=("system", "first"),
            P_thick=("P_thick", "median"),
            n_planets=("kepoi_name", "count"),
            min_period=("koi_period", "min"),
            max_period=("koi_period", "max"),
        )
        .reset_index()
    )
    rows["reason"] = "missing_alderaan_results_file"
    rows["priority"] = 50
    return rows


def build_suspicious_manifest(flagged: pd.DataFrame, n_suspicious: int) -> pd.DataFrame:
    flagged = flagged.sort_values(["n_flags", "e50"], ascending=False).copy()
    top = flagged.drop_duplicates("koi_target").head(n_suspicious)
    rows = top[
        [
            "koi_target",
            "kepoi_name",
            "disk",
            "system",
            "P_thick",
            "e50",
            "zeta_median",
            "n_flags",
            "alderaan_to_koi_duration_ratio",
            "alderaan_impact_med",
        ]
    ].copy()
    rows["kepid"] = np.nan
    rows["n_planets"] = np.nan
    rows["min_period"] = np.nan
    rows["max_period"] = np.nan
    rows["reason"] = "refit_suspicious_duration_ratio_tail"
    rows["priority"] = np.arange(1, len(rows) + 1)
    return rows


def build_control_manifest(shape: pd.DataFrame, n_per_bin: int) -> pd.DataFrame:
    clean = shape[
        (shape["zeta_median"].between(0.95, 1.05))
        & (shape["e50"] < 0.25)
        & (np.abs(np.log(shape["alderaan_to_koi_duration_ratio"])) <= np.log(1.10))
        & (shape["alderaan_impact_med"] < 0.8)
        & (shape["rho_frac_err"] < 0.2)
    ].copy()
    rows = []
    for (disk, system), sub in clean.groupby(["disk", "system"]):
        sub = sub.sort_values(["koi_model_snr", "kepoi_name"], ascending=[False, True]).drop_duplicates("koi_target")
        for _, row in sub.head(n_per_bin).iterrows():
            rows.append(
                {
                    "koi_target": row["koi_target"],
                    "kepid": np.nan,
                    "kepoi_name": row["kepoi_name"],
                    "disk": disk,
                    "system": system,
                    "P_thick": row.get("P_thick", np.nan),
                    "e50": row["e50"],
                    "zeta_median": row["zeta_median"],
                    "n_flags": 0,
                    "n_planets": np.nan,
                    "min_period": np.nan,
                    "max_period": np.nan,
                    "reason": "clean_control_existing_alderaan",
                    "priority": 100,
                }
            )
    return pd.DataFrame(rows)


def build_missing_validation_manifest(coverage: pd.DataFrame, sample: pd.DataFrame, n_per_bin: int) -> pd.DataFrame:
    missing = coverage[coverage["coverage_state"] == "no_results_file"].copy()
    rows = []
    for (disk, system), sub in missing.groupby(["disk", "system"]):
        sub = sub.merge(
            sample[["koi_target", "kepoi_name", "koi_model_snr", "koi_prad"]].drop_duplicates("kepoi_name"),
            on=["koi_target", "kepoi_name"],
            how="left",
            suffixes=("", "_sample"),
        )
        systems = (
            sub.groupby("koi_target")
            .agg(
                kepid=("kepid", "first"),
                P_thick=("P_thick", "median"),
                n_planets=("kepoi_name", "count"),
                min_period=("koi_period", "min"),
                max_period=("koi_period", "max"),
                max_snr=("koi_model_snr", "max"),
            )
            .reset_index()
            .sort_values(["max_snr", "koi_target"], ascending=[False, True])
        )
        for _, row in systems.head(n_per_bin).iterrows():
            rows.append(
                {
                    "koi_target": row["koi_target"],
                    "kepid": int(row["kepid"]),
                    "disk": disk,
                    "system": system,
                    "P_thick": row["P_thick"],
                    "n_planets": int(row["n_planets"]),
                    "min_period": row["min_period"],
                    "max_period": row["max_period"],
                    "reason": "missing_alderaan_validation_high_snr",
                    "priority": 25,
                }
            )
    return pd.DataFrame(rows)


def write_markdown(full: pd.DataFrame, validation: pd.DataFrame, missing: pd.DataFrame) -> None:
    path = output_dir() / "alderaan_refit_plan.md"
    lines = [
        "# ALDERAAN Refit / Missing-System Plan",
        "",
        "This plan separates three cases:",
        "",
        "- systems with no local ALDERAAN results at all;",
        "- systems with existing ALDERAAN results but suspicious duration-ratio tails;",
        "- clean existing systems to rerun as controls.",
        "",
        "## Counts",
        "",
        full.groupby(["batch", "reason"]).size().rename("systems").reset_index().to_markdown(index=False),
        "",
        "## Validation Batch",
        "",
        validation[
            [
                "koi_target",
                "disk",
                "system",
                "reason",
                "priority",
                "e50",
                "zeta_median",
                "n_flags",
                "n_planets",
            ]
        ]
        .head(80)
        .to_markdown(index=False, floatfmt=".4f"),
        "",
        "## Full Missing Counts",
        "",
        missing.groupby(["disk", "system"]).agg(systems=("koi_target", "count"), planets=("n_planets", "sum")).reset_index().to_markdown(index=False),
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    main()
