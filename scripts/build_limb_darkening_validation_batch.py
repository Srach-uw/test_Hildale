"""Build a small GCP ALDERAAN rerun bundle with reference limb-darkening priors.

The goal is an A/B validation, not a production rerun:

1. Choose a compact set of high-leverage systems whose current eccentricity
   posteriors drive the population mismatch, plus a few controls.
2. Rebuild the ALDERAAN input catalog for those full KOI systems, replacing the
   current KOI-derived limb-darkening centers with ALDERAAN's bundled
   Kepler-Gaia reference values.
3. Emit a self-contained cloud bundle that runs under a new run_id, preserving
   the existing results for direct comparison.
"""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path

import numpy as np
import pandas as pd

from common import koi_target, load_config, read_koi, root_path


RUN_ID = "sagear_ld_reference_validation"

# Availability audit performed against the public MAST Kepler light-curve
# archive on 2026-07-10. Counts are numbers of short-cadence files, not quarters.
SHORT_CADENCE_FILE_COUNTS = {
    7051180: 41,
    5695396: 36,
    8684730: 14,
    10418224: 0,
    7529266: 9,
    9846348: 1,
    12644822: 0,
    4544670: 0,
    6061119: 0,
    6526710: 0,
    7585481: 0,
    1871056: 24,
    10864656: 23,
    7951018: 0,
    5864975: 0,
    11499228: 0,
    11074835: 11,
    11098013: 4,
    12206313: 0,
    6937529: 0,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-targets", type=int, default=24)
    parser.add_argument("--random-seed", type=int, default=20260710)
    parser.add_argument("--out-dir", default="cloud_ld_validation_batch")
    parser.add_argument("--min-ld-delta", type=float, default=0.10)
    parser.add_argument(
        "--allow-reference-sentinel",
        action="store_true",
        help="Allow ALDERAAN bundled reference rows with u1=u2=0.1; default excludes them as likely fallback values.",
    )
    return parser.parse_args()


def repo_root() -> Path:
    return Path(__file__).resolve().parent


def load_reference_ld(cfg: dict) -> pd.DataFrame:
    alderaan_repo = root_path(cfg, "alderaan_repo")
    if alderaan_repo is None:
        alderaan_repo = repo_root().parent / "external" / "alderaan"
    ref_path = alderaan_repo / "Catalogs" / "kepler_dr25_gaia_dr2_crossmatch.csv"
    ref = pd.read_csv(ref_path)
    ref["kic_id"] = pd.to_numeric(ref["kic_id"], errors="coerce").astype("Int64")
    ref["limbdark_1"] = pd.to_numeric(ref["limbdark_1"], errors="coerce")
    ref["limbdark_2"] = pd.to_numeric(ref["limbdark_2"], errors="coerce")
    # Limb darkening is stellar/system-level. Use the median across planets for
    # a KIC to avoid planet-row quirks while preserving the stellar prior center.
    return (
        ref.dropna(subset=["kic_id", "limbdark_1", "limbdark_2"])
        .groupby("kic_id", as_index=False)
        .agg(
            ref_limbdark_1=("limbdark_1", "median"),
            ref_limbdark_2=("limbdark_2", "median"),
            ref_ld_rows=("limbdark_1", "size"),
        )
    )


def rank_validation_targets(
    summary: pd.DataFrame,
    context: pd.DataFrame,
    leverage: pd.DataFrame,
    min_ld_delta: float,
    n_targets: int,
    allow_reference_sentinel: bool,
    random_seed: int,
) -> pd.DataFrame:
    merged = summary.merge(
        context[
            [
                "kepoi_name",
                "ld_abs_max_delta",
                "ld_du1",
                "ld_du2",
                "limbdark_1",
                "limbdark_2",
                "ref_limbdark_1",
                "ref_limbdark_2",
            ]
        ],
        on="kepoi_name",
        how="left",
    )
    merged = merged.merge(
        leverage[["kepoi_name", "delta_loglike_map_minus_min"]],
        on="kepoi_name",
        how="left",
    )
    merged["population"] = merged["disk"].astype(str) + "_" + merged["system"].astype(str) + "s"
    merged["abs_zeta_offset"] = (pd.to_numeric(merged["zeta_median"], errors="coerce") - 1.0).abs()
    merged["e50"] = pd.to_numeric(merged["e50"], errors="coerce")
    merged["ld_abs_max_delta"] = pd.to_numeric(merged["ld_abs_max_delta"], errors="coerce")
    merged["delta_loglike_map_minus_min"] = pd.to_numeric(
        merged["delta_loglike_map_minus_min"], errors="coerce"
    )

    # Target only rows with existing ALDERAAN FITS from our cloud run for the
    # primary A/B test. Archive rows are useful later but do not guarantee raw
    # paired ALDERAAN shapes in the same format.
    new = merged[merged["posterior_source"].eq("new_alderaan")].copy()
    new = new[np.isfinite(new["ld_abs_max_delta"])]
    ref_sentinel = new["ref_limbdark_1"].round(4).eq(0.1) & new["ref_limbdark_2"].round(4).eq(0.1)
    if not allow_reference_sentinel:
        new = new[~ref_sentinel].copy()
    high_ld = new[new["ld_abs_max_delta"] >= min_ld_delta].copy()

    selected: list[pd.DataFrame] = []

    def take(label: str, query: pd.Series, n: int, sort_cols: list[str]) -> None:
        pool = high_ld[query].copy()
        if pool.empty:
            return
        pool = pool.sort_values(sort_cols, ascending=[False] * len(sort_cols))
        already = set(pd.concat(selected)["koi_target"]) if selected else set()
        pool = pool[~pool["koi_target"].isin(already)]
        picked = pool.drop_duplicates("koi_target").head(n).copy()
        if len(picked):
            picked["selection_reason"] = label
            selected.append(picked)

    take(
        "high_leverage_thin_single_short_zeta",
        high_ld["population"].eq("thin_singles") & (high_ld["zeta_median"] < 0.85),
        8,
        ["delta_loglike_map_minus_min", "e50"],
    )
    take(
        "high_leverage_thick_single_short_zeta",
        high_ld["population"].eq("thick_singles") & (high_ld["zeta_median"] < 0.85),
        3,
        ["delta_loglike_map_minus_min", "e50"],
    )
    take(
        "high_leverage_multis",
        high_ld["population"].isin(["thin_multis", "thick_multis"]),
        4,
        ["delta_loglike_map_minus_min", "abs_zeta_offset"],
    )
    take(
        "thin_single_long_zeta_control",
        high_ld["population"].eq("thin_singles") & (high_ld["zeta_median"] > 1.15),
        2,
        ["delta_loglike_map_minus_min", "e50"],
    )

    # Stable controls: low-ish e, zeta close to 1, but still with a meaningful
    # LD offset. If these move as much as the high-leverage systems, the LD
    # priors are globally important; if not, the tail is more specific.
    stable = high_ld[
        (high_ld["e50"] < 0.20)
        & (high_ld["zeta_median"].between(0.9, 1.1))
        & high_ld["population"].isin(["thin_singles", "thin_multis"])
    ].copy()
    if not stable.empty:
        already = set(pd.concat(selected)["koi_target"]) if selected else set()
        stable = stable[~stable["koi_target"].isin(already)]
        stable = stable.sort_values(["ld_abs_max_delta", "delta_loglike_map_minus_min"], ascending=False)
        stable = stable.drop_duplicates("koi_target").head(3).copy()
        stable["selection_reason"] = "stable_ld_offset_control"
        selected.append(stable)

    # The targeted tail is intentionally enriched for suspected failures. Add
    # one deterministic, otherwise-unselected system per population so the
    # factorial effects are not interpreted from pathological cases alone.
    already = set(pd.concat(selected)["koi_target"]) if selected else set()
    random_pool = new[~new["koi_target"].isin(already)].drop_duplicates("koi_target").copy()
    random_controls = []
    for offset, population in enumerate(["thin_singles", "thick_singles", "thin_multis", "thick_multis"]):
        pool = random_pool[random_pool["population"].eq(population)]
        if pool.empty:
            continue
        pick = pool.sample(n=1, random_state=random_seed + offset).copy()
        pick["selection_reason"] = f"random_{population}_control"
        random_controls.append(pick)
    if random_controls:
        selected.append(pd.concat(random_controls, ignore_index=True))

    if not selected:
        raise RuntimeError("No validation targets selected")
    out = pd.concat(selected, ignore_index=True)
    out = out.drop_duplicates("koi_target").head(n_targets).copy()
    out.insert(0, "validation_priority", np.arange(1, len(out) + 1))
    return out


def expand_to_full_koi_systems(base_catalog: pd.DataFrame, selected: pd.DataFrame, cfg: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build complete ALDERAAN systems before applying the population period cut.

    Sagear states that every planet in a system was sampled simultaneously.
    The production missing-target catalog was built after the 1-100 day cut,
    which omitted valid longer-period siblings in some systems. This expansion
    restores all non-false-positive KOIs with finite transit seeds while keeping
    each arm's limb-darkening center system-level.
    """
    raw = read_koi(cfg).copy()
    raw["koi_target"] = raw["kepoi_name"].map(koi_target)
    targets = set(selected["koi_target"].astype(str))
    raw = raw[
        raw["koi_target"].isin(targets)
        & raw["koi_disposition"].isin(["CONFIRMED", "CANDIDATE"])
    ].copy()

    seed_cols = ["koi_period", "koi_time0bk", "koi_depth", "koi_duration"]
    for col in seed_cols + ["koi_impact"]:
        raw[col] = pd.to_numeric(raw[col], errors="coerce")
    raw["seedable"] = raw[seed_cols].notna().all(axis=1) & (raw["koi_depth"] > 0) & (raw["koi_duration"] > 0)

    original_ld = (
        base_catalog.groupby("koi_id", as_index=False)
        .agg(limbdark_1=("limbdark_1", "first"), limbdark_2=("limbdark_2", "first"))
        .set_index("koi_id")
    )
    rows: list[dict[str, object]] = []
    inventory_rows: list[dict[str, object]] = []
    for target in sorted(targets):
        system = raw[raw["koi_target"].eq(target)].sort_values("koi_period")
        if target not in original_ld.index:
            raise RuntimeError(f"Original ALDERAAN limb-darkening center missing for {target}")
        u1 = float(original_ld.at[target, "limbdark_1"])
        u2 = float(original_ld.at[target, "limbdark_2"])
        seedable = system[system["seedable"]].copy()
        if seedable.empty:
            raise RuntimeError(f"No seedable KOIs for selected target {target}")
        npl = len(seedable)
        for _, planet in system.iterrows():
            inventory_rows.append(
                {
                    "koi_target": target,
                    "kepoi_name": planet["kepoi_name"],
                    "kepid": int(planet["kepid"]),
                    "period_days": planet["koi_period"],
                    "disposition": planet["koi_disposition"],
                    "seedable": bool(planet["seedable"]),
                    "included_in_alderaan_system": bool(planet["seedable"]),
                    "inside_population_period_cut": bool(
                        np.isfinite(planet["koi_period"]) and 1.0 <= planet["koi_period"] <= 100.0
                    ),
                }
            )
        for _, planet in seedable.iterrows():
            rows.append(
                {
                    "koi_id": target,
                    "kic_id": int(planet["kepid"]),
                    "npl": int(npl),
                    "period": float(planet["koi_period"]),
                    "epoch": float(planet["koi_time0bk"]),
                    "depth": float(planet["koi_depth"]),
                    "duration": float(planet["koi_duration"]),
                    "impact": float(planet["koi_impact"]) if np.isfinite(planet["koi_impact"]) else 0.5,
                    "limbdark_1": u1,
                    "limbdark_2": u2,
                }
            )
    catalog = pd.DataFrame(rows).sort_values(["koi_id", "period"]).reset_index(drop=True)
    inventory = pd.DataFrame(inventory_rows).sort_values(["koi_target", "period_days"]).reset_index(drop=True)
    return catalog, inventory


def build_reference_catalog(base_catalog: pd.DataFrame, selected: pd.DataFrame, ref_ld: pd.DataFrame) -> pd.DataFrame:
    targets = set(selected["koi_target"])
    cat = base_catalog[base_catalog["koi_id"].isin(targets)].copy()
    if cat.empty:
        raise RuntimeError("No selected targets found in base ALDERAAN catalog")
    cat["kic_id"] = pd.to_numeric(cat["kic_id"], errors="coerce").astype("Int64")
    ref = ref_ld.copy()
    merged = cat.merge(ref, on="kic_id", how="left")
    missing = merged["ref_limbdark_1"].isna() | merged["ref_limbdark_2"].isna()
    if missing.any():
        missing_targets = sorted(merged.loc[missing, "koi_id"].unique())
        raise RuntimeError(f"Reference limb darkening missing for selected targets: {missing_targets}")
    sentinel = merged["ref_limbdark_1"].round(4).eq(0.1) & merged["ref_limbdark_2"].round(4).eq(0.1)
    if sentinel.any():
        sentinel_targets = sorted(merged.loc[sentinel, "koi_id"].unique())
        raise RuntimeError(
            "Selected targets include ALDERAAN reference u1=u2=0.1 sentinel-like values; "
            f"rerun with --allow-reference-sentinel only for an explicit fallback-prior stress test: {sentinel_targets}"
        )
    merged["original_limbdark_1"] = merged["limbdark_1"]
    merged["original_limbdark_2"] = merged["limbdark_2"]
    merged["limbdark_1"] = merged["ref_limbdark_1"]
    merged["limbdark_2"] = merged["ref_limbdark_2"]
    keep = [
        "koi_id",
        "kic_id",
        "npl",
        "period",
        "epoch",
        "depth",
        "duration",
        "impact",
        "limbdark_1",
        "limbdark_2",
    ]
    return merged[keep].sort_values(["koi_id", "period"]).reset_index(drop=True)


def write_bundle(
    out_dir: Path,
    selected: pd.DataFrame,
    reference_catalog: pd.DataFrame,
    original_catalog: pd.DataFrame,
    inventory: pd.DataFrame,
) -> None:
    template_dir = repo_root() / "cloud_missing_batch"
    out_dir.mkdir(parents=True, exist_ok=True)
    for name in [
        "run_batch.sh",
        "run_one_target.sh",
        "setup_vm.sh",
        "summarize_progress.sh",
        "validate_bundle.py",
        "pack_results.sh",
        "create_gcp_spot_vm.sh",
        "patch_alderaan_repro.py",
        "patch_alderaan_paper_priors.py",
    ]:
        src = template_dir / name
        if src.exists():
            shutil.copy2(src, out_dir / name)

    target_rows = []
    base_by_target = original_catalog.groupby("koi_id")
    for i, row in selected.iterrows():
        grp = base_by_target.get_group(row["koi_target"])
        target_rows.append(
            {
                "target_index": i,
                "koi_target": row["koi_target"],
                "kepid": int(row["kepid"]),
                "disk": row["disk"],
                "system": row["system"],
                "P_thick": row["P_thick"],
                "validation_priority": int(row["validation_priority"]),
                "selection_reason": row["selection_reason"],
                "driving_kepoi_name": row["kepoi_name"],
                "driving_e50": row["e50"],
                "driving_zeta_median": row["zeta_median"],
                "driving_ld_abs_max_delta": row["ld_abs_max_delta"],
                "target_sample_planets": len(grp),
                "all_kepoi_names": ",".join(grp["koi_id"].astype(str)),
                "min_period": grp["period"].min(),
                "max_period": grp["period"].max(),
                "short_cadence_file_count": int(SHORT_CADENCE_FILE_COUNTS.get(int(row["kepid"]), -1)),
            }
        )

    targets = pd.DataFrame(target_rows)
    targets.to_csv(out_dir / "targets_ld_reference_validation.csv", index=False)
    original_catalog.to_csv(out_dir / "sagear_original_full_system_catalog.csv")
    reference_catalog.to_csv(out_dir / "sagear_ld_reference_catalog.csv")
    inventory.to_csv(out_dir / "full_system_inventory.csv", index=False)
    selected.to_csv(out_dir / "ld_reference_validation_selection.csv", index=False)

    sc_targets = targets[targets["short_cadence_file_count"] > 0].copy()
    sc_targets.to_csv(out_dir / "targets_short_cadence_validation.csv", index=False)

    repeat_parts = []
    for reason, n in [
        ("high_leverage_thin_single_short_zeta", 4),
        ("high_leverage_thick_single_short_zeta", 1),
        ("high_leverage_multis", 1),
        ("thin_single_long_zeta_control", 1),
        ("stable_ld_offset_control", 1),
    ]:
        repeat_parts.append(targets[targets["selection_reason"].eq(reason)].head(n))
    repeat_targets = pd.concat(repeat_parts, ignore_index=True).drop_duplicates("koi_target")
    repeat_targets.to_csv(out_dir / "targets_repeatability_validation.csv", index=False)

    readme = f"""# Factorial ALDERAAN Validation

This bundle supersedes the original one-arm 20-target screen. It isolates four
known non-equivalences instead of changing them together:

1. Full-system fitting: all non-false-positive, seedable KOIs are fit, including
   companions outside the 1-100 day population cut.
2. Limb darkening: original KOI centers versus the bundled Kepler-Gaia centers,
   with identical targets, cadence, full-system catalog, ALDERAAN commit and RNG seed.
3. Cadence: long-only versus long+short for the {len(sc_targets)} targets with
   short-cadence files in the MAST audit. A value of -1 in the manifest means
   cadence availability was not audited for the added random controls.
4. Transit-fit priors: pinned public ALDERAAN versus the priors printed in
   Sagear Table 1, on the {len(repeat_targets)} repeatability systems only.

The pinned ALDERAAN commit is `7443dff16b7f9092e14a6f0cc1f8948d457c9e0b`.
Run `bash setup_vm.sh` once in this bundle before launching any arm.

Recommended exploratory sequence on the GCP VM:

```bash
cd ~/sagear_ld_validation_batch
bash setup_vm.sh
JOBS=6 nohup bash run_validation_matrix.sh > validation_matrix.log 2>&1 &
```

The matrix contains:

- `original_lc`: {len(targets)} full systems, original limb darkening, long cadence;
- `reference_lc`: the same {len(targets)} systems with reference limb darkening;
- `original_lcsc` and `reference_lcsc`: the {len(sc_targets)} short-cadence systems;
- `original_lc_repeat`: {len(repeat_targets)} systems repeated with a different seed.
- `paper_priors_original_lc`: those same {len(repeat_targets)} systems using
  Sagear Table 1's C0/C1 and Rp/R* priors in a separate patched clone.

This is {2 * len(targets) + 2 * len(sc_targets) + 2 * len(repeat_targets)} target fits for the
current selection, not 20. `run_reference_lc.sh` launches only the reference-LD
arm and is exploratory by itself; the full matrix is required for attribution.

Progress:

```bash
cd ~/sagear_ld_validation_batch
pgrep -af "run_batch|run_one_target|parallel|detrend_and_estimate_ttvs|analyze_autocorrelated_noise|fit_transit_shape" || true
tail -n 100 validation_matrix.log
find projects -name '*-results.fits' | wc -l
```

Pack all arms only after the matrix completes:

```bash
bash pack_ld_validation_results.sh
```

Interpretation is paired. Compare `reference_lc-original_lc` for limb darkening,
`original_lc_repeat-original_lc` for sampler variability, and LC+SC minus LC
within each limb-darkening arm for cadence. Do not compare the reference arm
directly to the historical cloud FITS because those fits omitted some system
companions and were not seeded.

Compare `paper_priors_original_lc-original_lc` only on the repeatability target
subset. An effect is persuasive only when it is larger than the paired
repeatability variation; this arm tests a manuscript/public-code ambiguity and
is not the canonical default.
"""
    (out_dir / "README_LD_VALIDATION.md").write_text(readme, encoding="utf-8")

    arm_specs = {
        "original_lc": ("targets_ld_reference_validation.csv", "sagear_original_full_system_catalog.csv", "long", 0),
        "reference_lc": ("targets_ld_reference_validation.csv", "sagear_ld_reference_catalog.csv", "long", 0),
        "original_lcsc": ("targets_short_cadence_validation.csv", "sagear_original_full_system_catalog.csv", "both", 0),
        "reference_lcsc": ("targets_short_cadence_validation.csv", "sagear_ld_reference_catalog.csv", "both", 0),
        "original_lc_repeat": ("targets_repeatability_validation.csv", "sagear_original_full_system_catalog.csv", "long", 1000003),
        "paper_priors_original_lc": ("targets_repeatability_validation.csv", "sagear_original_full_system_catalog.csv", "long", 0),
    }
    for arm, (target_file, catalog_file, cadence_mode, seed_offset) in arm_specs.items():
        (out_dir / f"run_{arm}.sh").write_text(
            f"""#!/usr/bin/env bash
set -euo pipefail
export TARGET_CSV="{target_file}"
export CATALOG_SOURCE="{catalog_file}"
export CATALOG_NAME="sagear_validation_catalog.csv"
export RUN_ID="sagear_validation_{arm}"
export PROJECT_DIR="$PWD/projects/{arm}"
export ALDERAAN_REPO="${{ALDERAAN_REPO:-$HOME/alderaan_sagear_pinned}}"
export CADENCE_MODE="{cadence_mode}"
export SEED_OFFSET="{seed_offset}"
export JOBS="${{JOBS:-6}}"
bash run_batch.sh
""",
            encoding="utf-8",
        )
        if arm == "paper_priors_original_lc":
            arm_path = out_dir / f"run_{arm}.sh"
            arm_text = arm_path.read_text(encoding="utf-8").replace(
                'export ALDERAAN_REPO="${ALDERAAN_REPO:-$HOME/alderaan_sagear_pinned}"',
                'export ALDERAAN_REPO="${ALDERAAN_PAPER_PRIOR_REPO:-$HOME/alderaan_sagear_paper_priors}"',
            )
            arm_path.write_text(arm_text, encoding="utf-8")

    (out_dir / "run_ld_validation.sh").write_text(
        "#!/usr/bin/env bash\nset -euo pipefail\nbash run_reference_lc.sh\n",
        encoding="utf-8",
    )
    (out_dir / "run_validation_matrix.sh").write_text(
        """#!/usr/bin/env bash
set -euo pipefail
for arm in original_lc reference_lc original_lcsc reference_lcsc original_lc_repeat paper_priors_original_lc; do
  echo "[$(date -Is)] starting $arm"
  bash "run_${arm}.sh"
done
echo "[$(date -Is)] validation matrix complete"
""",
        encoding="utf-8",
    )
    (out_dir / "pack_ld_validation_results.sh").write_text(
        """#!/usr/bin/env bash
set -euo pipefail
OUT="alderaan_factorial_validation_$(date +%Y%m%d_%H%M%S).tar.gz"
tar -czf "$OUT" projects logs provenance *.csv README_LD_VALIDATION.md run_*.sh patch_alderaan_*.py
echo "$OUT"
""",
        encoding="utf-8",
    )


def main() -> None:
    args = parse_args()
    root = repo_root()
    cfg = load_config(root / "config.json")
    outputs = root / "outputs"

    summary = pd.read_csv(outputs / "eccentricity_posterior_summary_merged_paired_exact_qcprimary.csv")
    context = pd.read_csv(outputs / "limb_darkening_population_context_QC_PRIMARY.csv")
    leverage = pd.read_csv(outputs / "rayleigh_per_planet_leverage_PAIRED_EXACT_qcprimary.csv")
    filtered_base_catalog = pd.read_csv(root / "cloud_missing_batch" / "sagear_missing_catalog_FIXED.csv")
    ref_ld = load_reference_ld(cfg)

    selected = rank_validation_targets(
        summary,
        context,
        leverage,
        args.min_ld_delta,
        args.n_targets,
        args.allow_reference_sentinel,
        args.random_seed,
    )
    original_catalog, inventory = expand_to_full_koi_systems(filtered_base_catalog, selected, cfg)
    catalog = build_reference_catalog(original_catalog, selected, ref_ld)

    out_dir = root / args.out_dir
    write_bundle(out_dir, selected, catalog, original_catalog, inventory)

    # Additional audit artifacts in outputs for easy comparison.
    selected.to_csv(outputs / "ld_reference_validation_targets.csv", index=False)
    catalog.to_csv(outputs / "ld_reference_validation_catalog.csv")
    original_catalog.to_csv(outputs / "ld_original_full_system_validation_catalog.csv")
    inventory.to_csv(outputs / "ld_validation_full_system_inventory.csv", index=False)

    print(f"Selected {len(selected)} validation targets")
    print(selected[["validation_priority", "koi_target", "kepoi_name", "population", "selection_reason", "e50", "zeta_median", "ld_abs_max_delta"]].to_string(index=False))
    print()
    print(f"Wrote bundle: {out_dir}")
    print(f"Run id: {RUN_ID}")


if __name__ == "__main__":
    main()
