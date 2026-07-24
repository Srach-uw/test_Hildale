from __future__ import annotations

"""Audit and plot the Toomre coordinates against Sagear's published host table.

This is deliberately separate from the exploratory classifier plots.  The
published machine-readable table is the authority for host labels and the
velocities used for the host markers.  The audit also quantifies whether the
local direct Angus and Astropy coordinate columns reproduce those values.
"""

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def finite_mask(frame: pd.DataFrame, columns: list[str]) -> pd.Series:
    return np.isfinite(frame[columns].apply(pd.to_numeric, errors="coerce")).all(axis=1)


def summarize_coordinates(merged: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    published_vperp = np.hypot(merged["vr_used_kms"], merged["vz_used_kms"])
    comparisons = {
        "direct_angus": (-merged["V_phi"], merged["V_perp"]),
        "old_astropy": (merged["V_phi_astropy"], merged["V_perp_astropy"]),
        "position_derived": (merged["V_phi_geom"], merged["V_perp_geom"]),
    }
    for source, (vphi, vperp) in comparisons.items():
        dx = pd.to_numeric(vphi, errors="coerce") - merged["vphi_used_kms"]
        dy = pd.to_numeric(vperp, errors="coerce") - published_vperp
        for axis, delta in [("vphi", dx), ("vperp", dy)]:
            finite = np.isfinite(delta.to_numpy(float))
            values = np.abs(delta.to_numpy(float)[finite])
            rows.append(
                {
                    "source": source,
                    "axis": axis,
                    "n": int(finite.sum()),
                    "median_delta_kms": float(np.median(delta.to_numpy(float)[finite])) if finite.any() else np.nan,
                    "median_abs_delta_kms": float(np.median(values)) if values.size else np.nan,
                    "p95_abs_delta_kms": float(np.quantile(values, 0.95)) if values.size else np.nan,
                    "max_abs_delta_kms": float(values.max()) if values.size else np.nan,
                }
            )
    return pd.DataFrame(rows)


def make_plot(hosts: pd.DataFrame, background: pd.DataFrame | None, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(10, 7.4), constrained_layout=True)
    if background is not None and not background.empty:
        ax.scatter(
            background["vphi_plot"],
            background["vperp_plot"],
            s=2,
            c="0.72",
            alpha=0.16,
            linewidths=0,
            rasterized=True,
            label="Angus+2022 stars",
        )
    for disk, color, label in [
        ("thin", "#008b8b", "Published thin-disk hosts"),
        ("thick", "#9b1b1b", "Published thick-disk hosts"),
    ]:
        sub = hosts[hosts["disk_published"] == disk]
        ax.scatter(
            sub["vphi_used_kms"],
            sub["vperp_published"],
            marker="+",
            s=42,
            c=color,
            linewidths=1.35,
            alpha=0.86,
            label=label,
        )
    center = -229.0
    x = np.linspace(-350, -80, 900)
    for radius in [50, 100, 150, 200]:
        y2 = radius**2 - (x - center) ** 2
        ok = y2 >= 0
        ax.plot(x[ok], np.sqrt(y2[ok]), color="0.55", ls="--", lw=1.1, alpha=0.7)
    ax.set_xlim(-330, -100)
    ax.set_ylim(0, 200)
    ax.set_xlabel(r"$V_\phi$ [km s$^{-1}$]")
    ax.set_ylabel(r"$\sqrt{V_R^2 + V_Z^2}$ [km s$^{-1}$]")
    ax.set_title("Toomre diagram using Sagear published host velocities and labels")
    ax.legend(loc="upper left", frameon=True)
    fig.savefig(path, dpi=240)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", default=None)
    parser.add_argument("--inventory", default=None)
    parser.add_argument("--current-sample", default=None)
    parser.add_argument("--published-hosts", default=None)
    parser.add_argument("--background", default=None)
    args = parser.parse_args()

    root = Path(args.root) if args.root else Path(__file__).resolve().parent
    out = root / "outputs"
    inventory_path = Path(args.inventory) if args.inventory else out / "sagear2026_planet_inventory_pre_visual_qc.csv"
    current_path = Path(args.current_sample) if args.current_sample else out / "canonical_sample_diagnostic.csv"
    published_path = Path(args.published_hosts) if args.published_hosts else out / "sagear2026_published_kinematic_hosts.csv"

    inventory = pd.read_csv(inventory_path)
    current = pd.read_csv(current_path).drop_duplicates("kepid")
    inventory_host = inventory.drop_duplicates("kepid")
    published = pd.read_csv(published_path)
    if published["kepid"].duplicated().any():
        raise ValueError("Published host table must contain one row per KIC host")

    required_pub = ["kepid", "vphi_used_kms", "vr_used_kms", "vz_used_kms", "p_thick_published", "disk_published"]
    missing = sorted(set(required_pub) - set(published.columns))
    if missing:
        raise ValueError(f"Published host table missing columns: {missing}")
    required_current = ["kepid", "V_phi", "V_perp", "V_phi_astropy", "V_perp_astropy", "V_phi_geom", "V_perp_geom", "disk"]
    missing = sorted(set(required_current) - set(current.columns))
    if missing:
        raise ValueError(f"Current sample missing coordinate columns: {missing}")

    merged = published.merge(current[required_current], on="kepid", how="inner", validate="one_to_one")
    inventory_cols = [c for c in ["kepid", "disk", "disk_published", "p_thick_published"] if c in inventory_host]
    merged = merged.merge(inventory_host[inventory_cols], on="kepid", how="left", suffixes=("", "_inventory"))
    if merged.empty:
        raise ValueError("No published hosts have current coordinates")
    merged["vperp_published"] = np.hypot(merged["vr_used_kms"], merged["vz_used_kms"])
    merged["label_agrees_current"] = merged["disk_published"].eq(merged["disk"])
    inventory_label = merged.get("disk_published_inventory", merged.get("disk_inventory"))
    merged["label_agrees_inventory"] = merged["disk_published"].eq(inventory_label)

    summary = summarize_coordinates(merged)
    summary.to_csv(out / "published_toomre_coordinate_audit.csv", index=False)
    merged.to_csv(out / "published_toomre_host_join_audit.csv", index=False)

    background = None
    if args.background:
        bg = pd.read_csv(args.background)
        if {"V_phi_astropy", "V_perp_astropy"}.issubset(bg.columns):
            background = bg[finite_mask(bg, ["V_phi_astropy", "V_perp_astropy"])].copy()
            background["vphi_plot"] = background["V_phi_astropy"]
            background["vperp_plot"] = background["V_perp_astropy"]
    make_plot(merged, background, out / "toomre_sagear_published_truth.png")

    metadata = {
        "inventory": str(inventory_path),
        "current_sample": str(current_path),
        "published_hosts": str(published_path),
        "n_published_hosts": int(len(published)),
        "n_joined_hosts": int(len(merged)),
        "n_missing_current_coordinates": int(len(published) - len(merged)),
        "current_label_agreement": int(merged["label_agrees_current"].sum()),
        "inventory_label_agreement": int(merged["label_agrees_inventory"].sum()),
        "coordinate_authority": "published vphi_used_kms, vr_used_kms, vz_used_kms",
        "warning": "Do not use the direct Angus proxy Toomre panel as a Sagear replication figure.",
    }
    (out / "published_toomre_truth_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print(summary.to_string(index=False))
    print(
        f"published_hosts={len(published)} joined={len(merged)} "
        f"current_labels_agree={int(merged['label_agrees_current'].sum())} "
        f"inventory_labels_agree={int(merged['label_agrees_inventory'].sum())}"
    )


if __name__ == "__main__":
    main()
