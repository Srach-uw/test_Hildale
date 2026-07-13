from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import load_config, output_dir, read_sagear2026_kinematic_hosts


EXPECTED = {"hosts": 1888, "thin": 1515, "thick": 373, "measured_velocity": 585}


def validate_published_hosts(hosts: pd.DataFrame) -> dict[str, int]:
    counts = {
        "hosts": len(hosts),
        "thin": int((hosts["disk_published"] == "thin").sum()),
        "thick": int((hosts["disk_published"] == "thick").sum()),
        "measured_velocity": int(hosts["has_measured_velocity"].sum()),
    }
    if counts != EXPECTED:
        raise ValueError(f"Published table count contract failed: {counts} != {EXPECTED}")
    threshold = np.where(hosts["p_thick_published"] > 0.5, "thick", "thin")
    if not np.array_equal(threshold, hosts["disk_published"].to_numpy()):
        raise ValueError("Published disk labels do not match P_thick > 0.5")
    return counts


def reconcile(hosts: pd.DataFrame, sample_path: Path) -> tuple[pd.DataFrame, dict[str, float]]:
    sample = pd.read_csv(sample_path)
    if "kepid" not in sample:
        raise ValueError(f"Sample has no kepid column: {sample_path}")
    sample["kepid"] = pd.to_numeric(sample["kepid"], errors="coerce").astype("Int64")
    own_cols = ["kepid"]
    for candidate in ("disk", "p_thick", "vr", "vphi", "vz", "V_R", "V_phi", "V_z"):
        if candidate in sample and candidate not in own_cols:
            own_cols.append(candidate)
    ours = sample[own_cols].dropna(subset=["kepid"]).drop_duplicates("kepid").copy()
    ours["kepid"] = ours["kepid"].astype(int)
    merged = hosts.merge(ours, on="kepid", how="outer", indicator=True)
    both = merged["_merge"] == "both"
    stats: dict[str, float] = {
        "sample_planets": float(len(sample)),
        "sample_hosts": float(len(ours)),
        "overlap_hosts": float(both.sum()),
        "sample_only_hosts": float((merged["_merge"] == "right_only").sum()),
        "published_only_hosts": float((merged["_merge"] == "left_only").sum()),
    }
    if "disk" in merged:
        ours_disk = merged["disk"].astype(str).str.lower()
        valid = both & ours_disk.isin({"thin", "thick"})
        stats["disk_agreement_fraction"] = float((ours_disk[valid] == merged.loc[valid, "disk_published"]).mean())
        stats["disk_disagreements"] = float((ours_disk[valid] != merged.loc[valid, "disk_published"]).sum())
        merged["disk_label_agrees"] = np.where(valid, ours_disk == merged["disk_published"], pd.NA)
    if "p_thick" in merged:
        valid = both & np.isfinite(pd.to_numeric(merged["p_thick"], errors="coerce"))
        x = pd.to_numeric(merged.loc[valid, "p_thick"], errors="coerce")
        y = merged.loc[valid, "p_thick_published"]
        stats["p_thick_correlation"] = float(x.corr(y))
        stats["p_thick_mae"] = float(np.mean(np.abs(x - y)))
    return merged, stats


def relabel_planets(hosts: pd.DataFrame, sample_path: Path) -> tuple[pd.DataFrame, dict[str, int]]:
    planets = pd.read_csv(sample_path)
    if "system" not in planets:
        raise ValueError("Planet sample lacks the canonical pre-cut system label")
    planets["kepid"] = pd.to_numeric(planets["kepid"], errors="coerce").astype("Int64")
    published = hosts[["kepid", "disk_published", "p_thick_published"]].copy()
    relabeled = planets.merge(published, on="kepid", how="inner", validate="many_to_one").copy()
    if "disk" in relabeled:
        relabeled = relabeled.rename(columns={"disk": "disk_reconstructed"})
    if "p_thick" in relabeled:
        relabeled = relabeled.rename(columns={"p_thick": "p_thick_reconstructed"})
    relabeled = relabeled.copy()
    relabeled["disk"] = relabeled["disk_published"]
    # Disk labels come from the published host table, but multiplicity remains
    # the pre-cut architecture carried by the canonical planet sample. Keep an
    # overlap-only recount solely to audit the historical notebook failure mode.
    counts_per_host = relabeled.groupby("kepid")["kepid"].transform("size")
    relabeled["system_overlap_recount"] = np.where(counts_per_host == 1, "single", "multi")
    relabeled["system_published_overlap"] = relabeled["system"].astype(str)
    relabeled["multiplicity_recount_disagrees"] = (
        relabeled["system_overlap_recount"] != relabeled["system_published_overlap"]
    )
    counts: dict[str, int] = {
        "published_relabel_planets": len(relabeled),
        "published_relabel_hosts": int(relabeled["kepid"].nunique()),
        "published_relabel_multiplicity_recount_disagreements": int(
            relabeled["multiplicity_recount_disagrees"].sum()
        ),
    }
    for disk in ("thin", "thick"):
        for system in ("single", "multi"):
            key = f"published_relabel_{disk}_{system}_planets"
            counts[key] = int(((relabeled["disk"] == disk) & (relabeled["system_published_overlap"] == system)).sum())
    return relabeled, counts


def plot_toomre(hosts: pd.DataFrame, path: Path) -> None:
    transverse = np.hypot(hosts["vr_used_kms"], hosts["vz_used_kms"])
    fig, ax = plt.subplots(figsize=(8.2, 7.2), constrained_layout=True)
    thin = hosts["disk_published"] == "thin"
    ax.scatter(hosts.loc[thin, "vphi_used_kms"], transverse[thin], s=9, alpha=0.42, color="#168b91", label=f"Thin ({thin.sum()})")
    ax.scatter(hosts.loc[~thin, "vphi_used_kms"], transverse[~thin], s=13, alpha=0.62, color="#a5231f", label=f"Thick ({(~thin).sum()})")
    ax.set(xlabel=r"$V_\phi$ [km s$^{-1}$]", ylabel=r"$\sqrt{V_R^2 + V_Z^2}$ [km s$^{-1}$]", title="Sagear et al. (2026) Published Host Classification")
    ax.legend(frameon=False)
    ax.grid(alpha=0.16)
    fig.savefig(path, dpi=220)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Audit the replication against Sagear et al. 2026 published Table 1.")
    parser.add_argument("--config", default=None)
    parser.add_argument("--sample", action="append", default=[], help="Planet-level CSV to reconcile; may be repeated.")
    args = parser.parse_args()
    cfg = load_config(args.config)
    out = output_dir()
    hosts = read_sagear2026_kinematic_hosts(cfg)
    counts = validate_published_hosts(hosts)
    hosts.to_csv(out / "sagear2026_published_kinematic_hosts.csv", index=False)
    plot_toomre(hosts, out / "toomre_sagear2026_published_truth.png")

    rows = [{"metric": key, "value": value, "sample": "published_table1"} for key, value in counts.items()]
    for raw in args.sample:
        sample_path = Path(raw).expanduser().resolve()
        merged, stats = reconcile(hosts, sample_path)
        stem = sample_path.stem.replace(" ", "_")
        merged.to_csv(out / f"sagear2026_host_reconciliation_{stem}.csv", index=False)
        rows.extend({"metric": key, "value": value, "sample": stem} for key, value in stats.items())
        relabeled, planet_counts = relabel_planets(hosts, sample_path)
        relabeled.to_csv(out / f"sagear2026_published_relabel_{stem}.csv", index=False)
        rows.extend({"metric": key, "value": value, "sample": stem} for key, value in planet_counts.items())
    summary = pd.DataFrame(rows)
    summary.to_csv(out / "sagear2026_publication_audit_summary.csv", index=False)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
