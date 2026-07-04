from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

from common import load_config, output_dir, root_path


def main() -> None:
    parser = argparse.ArgumentParser(description="Prepare ALDERAAN validation batches.")
    sub = parser.add_subparsers(dest="command", required=True)

    sub.add_parser("check-env")

    prep = sub.add_parser("prepare")
    prep.add_argument("--sample", default=None, help="canonical_sample.csv path")
    prep.add_argument("--n-per-bin", type=int, default=3)
    prep.add_argument("--n-high-e", type=int, default=5, help="Add this many high-e thin-single stress-test targets.")
    prep.add_argument("--config", default=None)

    dl = sub.add_parser("write-download-script")
    dl.add_argument("--targets", default=None)
    dl.add_argument("--config", default=None)

    run = sub.add_parser("write-run-script")
    run.add_argument("--targets", default=None)
    run.add_argument("--config", default=None)

    args = parser.parse_args()
    cfg = load_config(getattr(args, "config", None))

    if args.command == "check-env":
        check_env()
    elif args.command == "prepare":
        prepare(cfg, args.sample, args.n_per_bin, args.n_high_e)
    elif args.command == "write-download-script":
        write_download_script(cfg, args.targets)
    elif args.command == "write-run-script":
        write_run_script(cfg, args.targets)


def check_env() -> None:
    modules = [
        "astropy",
        "lightkurve",
        "pymc3",
        "exoplanet",
        "dynesty",
        "batman",
        "celerite2",
        "ldtk",
        "arviz",
    ]
    import importlib.util

    for mod in modules:
        print(f"{mod:12s} {'OK' if importlib.util.find_spec(mod) else 'MISSING'}")


def prepare(cfg: dict, sample_path: str | None, n_per_bin: int, n_high_e: int) -> None:
    out = output_dir()
    sample = Path(sample_path) if sample_path else out / "canonical_sample_diagnostic.csv"
    if not sample.exists():
        raise FileNotFoundError(
            "Missing classified sample for validation target selection. "
            f"Run diagnose_sample.py --allow-fallback-gmm first, or pass --sample explicitly: {sample}"
        )
    df = pd.read_csv(sample)
    targets = select_validation_targets(df, n_per_bin, n_high_e)

    project = root_path(cfg, "alderaan_project")
    assert project is not None
    for rel in ["Catalogs", "Data", "Results", "Figures", "Scripts"]:
        (project / rel).mkdir(parents=True, exist_ok=True)

    catalog = build_alderaan_catalog(df[df["koi_target"].isin(targets["koi_target"])], cfg)
    catalog_path = project / "Catalogs" / "sagear_validation_catalog.csv"
    catalog.to_csv(catalog_path)

    targets_path = project / "sagear_validation_targets.csv"
    targets.to_csv(targets_path, index=False)

    write_download_script(cfg, str(targets_path))
    write_run_script(cfg, str(targets_path))

    print(f"Wrote ALDERAAN validation project: {project}")
    print(f"Wrote target list: {targets_path}")
    print(f"Wrote catalog: {catalog_path}")


def select_validation_targets(df: pd.DataFrame, n_per_bin: int, n_high_e: int) -> pd.DataFrame:
    required = ["koi_period", "koi_time0bk", "koi_depth", "koi_duration"]
    valid = df.copy()
    for col in required:
        valid = valid[np.isfinite(pd.to_numeric(valid[col], errors="coerce"))]
    valid = valid[(valid["koi_depth"] > 0) & (valid["koi_duration"] > 0)]
    rows = []
    for disk in ["thin", "thick"]:
        for system in ["single", "multi"]:
            sub = valid[(valid["disk"] == disk) & (valid["system"] == system)].copy()
            sub = sub.sort_values(["P_thick", "kepid"], ascending=[disk == "thin", True])
            targets = sub.drop_duplicates("koi_target").head(n_per_bin)
            for _, row in targets.iterrows():
                rows.append(
                    {
                        "koi_target": row["koi_target"],
                        "kepid": int(row["kepid"]),
                        "disk": disk,
                        "system": system,
                        "P_thick": row.get("P_thick", np.nan),
                        "reason": f"{disk}_{system}_validation",
                    }
                )
    rows.extend(select_high_e_thin_single_targets(valid, n_high_e))
    return pd.DataFrame(rows).drop_duplicates("koi_target")


def select_high_e_thin_single_targets(valid: pd.DataFrame, n_high_e: int) -> list[dict]:
    if n_high_e <= 0:
        return []
    top_path = output_dir() / "eccentricity_diagnostics_thin_single_top_outliers.csv"
    if not top_path.exists():
        return []
    top = pd.read_csv(top_path)
    if "kepoi_name" not in top.columns or "e_value" not in top.columns:
        return []
    candidates = valid[(valid["disk"] == "thin") & (valid["system"] == "single")].merge(
        top[["kepoi_name", "e_value"]],
        on="kepoi_name",
        how="inner",
    )
    candidates = candidates.sort_values(["e_value", "kepid"], ascending=[False, True]).drop_duplicates("koi_target")
    rows = []
    for _, row in candidates.head(n_high_e).iterrows():
        rows.append(
            {
                "koi_target": row["koi_target"],
                "kepid": int(row["kepid"]),
                "disk": row["disk"],
                "system": row["system"],
                "P_thick": row.get("P_thick", np.nan),
                "reason": "thin_single_high_e_stress_test",
            }
        )
    return rows


def build_alderaan_catalog(df: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    rows = []
    defaults = cfg["alderaan"]
    required = ["koi_period", "koi_time0bk", "koi_depth", "koi_duration"]
    for col in required:
        df = df[np.isfinite(pd.to_numeric(df[col], errors="coerce"))]
    df = df[(df["koi_depth"] > 0) & (df["koi_duration"] > 0)]
    for target, grp in df.groupby("koi_target"):
        grp = grp.sort_values("koi_period")
        npl = len(grp)
        # Limb darkening is a stellar (system-wide) property, not per-planet.
        # ALDERAAN's io.parse_catalog hard-requires identical LD_U1/LD_U2 across
        # every planet in a system and raises ValueError otherwise. Computing it
        # per-row let differing/missing koi_ldm_coeff values across sibling
        # planets produce inconsistent rows within the same system - found live
        # via "There are inconsistencies with LD_U1" crashes on the GCP run,
        # 2026-07-04. Compute once per system: prefer the first row with a
        # finite koi_ldm_coeff pair, falling back to defaults only if none of
        # the system's planets have one.
        u1 = u2 = None
        for _, cand in grp.iterrows():
            cu1, cu2 = limb_darkening_from_koi(cand, defaults)
            if u1 is None or (
                np.isfinite(cand.get("koi_ldm_coeff1", np.nan))
                and np.isfinite(cand.get("koi_ldm_coeff2", np.nan))
            ):
                u1, u2 = cu1, cu2
                if np.isfinite(cand.get("koi_ldm_coeff1", np.nan)) and np.isfinite(
                    cand.get("koi_ldm_coeff2", np.nan)
                ):
                    break
        for _, p in grp.iterrows():
            rows.append(
                {
                    "koi_id": target,
                    "kic_id": int(p["kepid"]),
                    "npl": int(npl),
                    "period": float(p["koi_period"]),
                    "epoch": float(p["koi_time0bk"]),
                    "depth": float(p["koi_depth"]),
                    "duration": float(p["koi_duration"]),
                    # 0.5 is a neutral seed when koi_impact is missing: ALDERAAN fits the
                    # impact parameter b under a prior, so this only initializes the search
                    # rather than fixing b, and keeps a missing b from blocking cataloging.
                    "impact": float(p["koi_impact"]) if np.isfinite(p["koi_impact"]) else 0.5,
                    "limbdark_1": float(u1),
                    "limbdark_2": float(u2),
                }
            )
    return pd.DataFrame(rows)


def limb_darkening_from_koi(row: pd.Series, defaults: dict) -> tuple[float, float]:
    c1 = row.get("koi_ldm_coeff1", np.nan)
    c2 = row.get("koi_ldm_coeff2", np.nan)
    if np.isfinite(c1) and np.isfinite(c2):
        return float(c1), float(c2)
    return float(defaults["default_limbdark_1"]), float(defaults["default_limbdark_2"])


def write_download_script(cfg: dict, targets_path: str | None) -> None:
    project = root_path(cfg, "alderaan_project")
    repo = root_path(cfg, "alderaan_repo")
    assert project is not None
    targets = pd.read_csv(targets_path) if targets_path else pd.read_csv(project / "sagear_validation_targets.csv")
    script = project / "Scripts" / "download_validation_lightcurves.ps1"
    data_dir = project / "Data"
    data_dir.mkdir(parents=True, exist_ok=True)
    lines = [
        "$ErrorActionPreference = 'Continue'",
        f"$DataDir = '{data_dir}'",
        "New-Item -ItemType Directory -Force -Path $DataDir | Out-Null",
    ]
    if repo and (repo / "bin" / "get_kepler_data.py").exists():
        lines += [
            f"$AlderaanGet = '{repo / 'bin' / 'get_kepler_data.py'}'",
            f"Push-Location '{data_dir}'",
        ]
        for kepid in targets["kepid"].drop_duplicates():
            lines.append(f"python $AlderaanGet {int(kepid)} -c long -t lightcurve -o get_{int(kepid)}_lc.sh --cmdtype curl")
            lines.append(f"bash get_{int(kepid)}_lc.sh")
        lines.append("Pop-Location")
    else:
        lines.append("# Clone ALDERAAN first: git clone https://github.com/gjgilbert/alderaan external\\alderaan")
        lines.append("# Then rerun this command generator.")
    script.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote download script: {script}")


def write_run_script(cfg: dict, targets_path: str | None) -> None:
    project = root_path(cfg, "alderaan_project")
    repo = root_path(cfg, "alderaan_repo")
    assert project is not None
    targets = pd.read_csv(targets_path) if targets_path else pd.read_csv(project / "sagear_validation_targets.csv")
    run_id = cfg["alderaan"]["run_id"]
    mission = cfg["alderaan"]["mission"]
    script = project / "Scripts" / "run_validation_alderaan.ps1"
    lines = [
        "$ErrorActionPreference = 'Stop'",
        "# Run this from an activated ALDERAAN conda environment.",
        f"$ProjectDir = '{project}'",
    ]
    if repo and repo.exists():
        lines.append(f"Push-Location '{repo}'")
        for target in targets["koi_target"].drop_duplicates():
            lines += [
                f"python bin/detrend_and_estimate_ttvs.py --mission {mission} --target {target} --run_id {run_id} --project_dir $ProjectDir --data_dir \"$ProjectDir\\Data\\\" --catalog sagear_validation_catalog.csv",
                f"python bin/analyze_autocorrelated_noise.py --mission {mission} --target {target} --run_id {run_id} --project_dir $ProjectDir",
                f"python bin/fit_transit_shape_simultaneous_nested.py --mission {mission} --target {target} --run_id {run_id} --project_dir $ProjectDir",
            ]
        lines.append("Pop-Location")
    else:
        lines.append("# Clone ALDERAAN first: git clone https://github.com/gjgilbert/alderaan external\\alderaan")
    script.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote ALDERAAN run script: {script}")


if __name__ == "__main__":
    main()
