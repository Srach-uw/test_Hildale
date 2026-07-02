from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from classifier_threshold_diagnostics import build_probability_fields
from common import load_config, output_dir
from toomre_diagnostics import ensure_velocity_aliases


def main() -> None:
    parser = argparse.ArgumentParser(description="List systems whose disk labels change across near-match classifiers.")
    parser.add_argument("--config", default=None)
    parser.add_argument("--sample", default=None)
    parser.add_argument("--out-prefix", default="classifier_disagreement_diagnostics")
    args = parser.parse_args()

    cfg = load_config(args.config)
    out = output_dir()
    sample_path = Path(args.sample) if args.sample else out / "canonical_sample_diagnostic.csv"
    sample = ensure_velocity_aliases(pd.read_csv(sample_path))
    fields = {name: probs for name, probs, _note in build_probability_fields(sample, cfg)}

    classifiers = {
        "direct_gmm_p0p50": ("planet_host_direct_gmm", 0.50),
        "direct_gmm_p0p61": ("planet_host_direct_gmm", 0.61),
        "geom_gmm_p0p525": ("planet_host_geom_gmm", 0.525),
        "chem_1thin_1thick_p0p49": ("chem_1thin_1thick_prior_0p425", 0.49),
        "chem_2thin_2thick_p0p46": ("chem_2thin_2thick_prior_0p350", 0.46),
    }
    rows = sample[
        [
            c
            for c in [
                "kepid",
                "kepoi_name",
                "koi_target",
                "system",
                "koi_period",
                "koi_prad",
                "V_phi",
                "V_perp",
                "V_phi_geom",
                "V_perp_geom",
            ]
            if c in sample.columns
        ]
    ].copy()

    for label, (field_name, threshold) in classifiers.items():
        probs = fields.get(field_name)
        if probs is None:
            continue
        rows[f"P_thick_{label}"] = probs
        rows[f"disk_{label}"] = label_disks(probs, threshold)

    disk_cols = [c for c in rows.columns if c.startswith("disk_")]
    p_cols = [c for c in rows.columns if c.startswith("P_thick_")]
    rows["n_distinct_labels"] = rows[disk_cols].nunique(axis=1)
    rows["p_thick_min"] = rows[p_cols].min(axis=1)
    rows["p_thick_max"] = rows[p_cols].max(axis=1)
    rows["p_thick_range"] = rows["p_thick_max"] - rows["p_thick_min"]

    planet_path = out / f"{args.out_prefix}_planet_labels.csv"
    rows.to_csv(planet_path, index=False)

    changed = rows[rows["n_distinct_labels"] > 1].copy()
    changed_path = out / f"{args.out_prefix}_changed_planets.csv"
    changed.to_csv(changed_path, index=False)

    system = summarize_systems(rows, disk_cols)
    system_path = out / f"{args.out_prefix}_system_summary.csv"
    system.to_csv(system_path, index=False)

    pairwise = pairwise_summary(rows, disk_cols)
    pairwise_path = out / f"{args.out_prefix}_pairwise_summary.csv"
    pairwise.to_csv(pairwise_path, index=False)

    print("=== Disagreement Summary ===")
    print(
        pd.DataFrame(
            [
                {
                    "planets": len(rows),
                    "changed_planets": len(changed),
                    "changed_fraction": len(changed) / len(rows) if len(rows) else np.nan,
                    "systems": rows["kepid"].nunique(),
                    "changed_systems": changed["kepid"].nunique(),
                }
            ]
        ).to_string(index=False)
    )
    print("\n=== Changed planets by system ===")
    print(changed.groupby("system").size().rename("changed_planets").reset_index().to_string(index=False))
    print("\n=== Pairwise classifier changes ===")
    print(pairwise.to_string(index=False))
    print(f"\nWrote: {planet_path}")
    print(f"Wrote: {changed_path}")
    print(f"Wrote: {system_path}")
    print(f"Wrote: {pairwise_path}")


def label_disks(probs: np.ndarray, threshold: float) -> np.ndarray:
    labels = np.where(probs > threshold, "thick", "thin").astype(object)
    labels[~np.isfinite(probs)] = "unclassified"
    return labels


def summarize_systems(rows: pd.DataFrame, disk_cols: list[str]) -> pd.DataFrame:
    grouped = rows.groupby(["kepid", "koi_target", "system"], dropna=False)
    out_rows = []
    for (kepid, target, system), grp in grouped:
        row = {"kepid": kepid, "koi_target": target, "system": system, "n_planets": len(grp)}
        for col in disk_cols:
            vals = sorted(set(grp[col].dropna()))
            row[col] = ",".join(vals)
        row["any_classifier_disagreement"] = any("," in str(row[col]) for col in disk_cols) or grp["n_distinct_labels"].max() > 1
        row["max_p_thick_range"] = float(grp["p_thick_range"].max())
        out_rows.append(row)
    return pd.DataFrame(out_rows).sort_values(["any_classifier_disagreement", "max_p_thick_range"], ascending=[False, False])


def pairwise_summary(rows: pd.DataFrame, disk_cols: list[str]) -> pd.DataFrame:
    out = []
    for i, left in enumerate(disk_cols):
        for right in disk_cols[i + 1 :]:
            changed = rows[rows[left] != rows[right]]
            out.append(
                {
                    "left_classifier": left.removeprefix("disk_"),
                    "right_classifier": right.removeprefix("disk_"),
                    "changed_planets": len(changed),
                    "changed_systems": changed["kepid"].nunique(),
                    "changed_singles": int((changed["system"] == "single").sum()),
                    "changed_multis": int((changed["system"] == "multi").sum()),
                }
            )
    return pd.DataFrame(out).sort_values("changed_planets")


if __name__ == "__main__":
    main()
