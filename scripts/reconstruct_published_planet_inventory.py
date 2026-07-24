from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from common import (
    load_config,
    output_dir,
    read_berger2018_stellar_table,
    read_berger_table2,
    read_koi,
    read_sagear2026_kinematic_hosts,
    root_path,
)


PAPER_COUNTS = {
    ("thin", "single"): 1121,
    ("thick", "single"): 275,
    ("thin", "multi"): 862,
    ("thick", "multi"): 207,
}
PRE_VISUAL_COUNTS = {
    ("thin", "single"): 1121,
    ("thick", "single"): 275,
    ("thin", "multi"): 883,
    ("thick", "multi"): 212,
}


def normalize_planet_catalog(koi: pd.DataFrame) -> pd.DataFrame:
    required = {
        "kepid",
        "kepoi_name",
        "koi_disposition",
        "koi_period",
        "koi_count",
    }
    missing = sorted(required - set(koi.columns))
    if missing:
        raise ValueError(f"KOI catalog is missing required columns: {missing}")
    out = koi.copy()
    out["kepid"] = pd.to_numeric(out["kepid"], errors="coerce").astype("Int64")
    out["koi_period"] = pd.to_numeric(out["koi_period"], errors="coerce")
    out["koi_count"] = pd.to_numeric(out["koi_count"], errors="coerce")
    out["kepoi_name"] = out["kepoi_name"].astype(str).str.strip()
    out["koi_target"] = out["kepoi_name"].str.extract(
        r"^(K\d{5})", expand=False
    )
    return out


def build_pre_visual_inventory(
    koi: pd.DataFrame, hosts: pd.DataFrame
) -> pd.DataFrame:
    koi = normalize_planet_catalog(koi)
    host_ids = set(pd.to_numeric(hosts["kepid"], errors="coerce").dropna().astype(int))
    selected = koi.loc[
        koi["kepid"].isin(host_ids)
        & koi["koi_disposition"].astype(str).str.upper().ne("FALSE POSITIVE")
        & koi["koi_period"].between(1.0, 100.0, inclusive="both")
    ].copy()
    if selected["kepoi_name"].duplicated().any():
        duplicates = sorted(
            selected.loc[selected["kepoi_name"].duplicated(), "kepoi_name"].unique()
        )
        raise ValueError(f"Duplicate KOI planet identifiers: {duplicates[:10]}")
    if selected["koi_target"].isna().any():
        bad = selected.loc[selected["koi_target"].isna(), "kepoi_name"].tolist()
        raise ValueError(f"Cannot derive KOI target for: {bad[:10]}")
    if selected["koi_count"].isna().any():
        bad = selected.loc[selected["koi_count"].isna(), "kepoi_name"].tolist()
        raise ValueError(f"Missing pre-cut KOI multiplicity for: {bad[:10]}")

    selected["system"] = np.where(selected["koi_count"] == 1, "single", "multi")
    selected = selected.merge(
        hosts[
            [
                "kepid",
                "disk_published",
                "p_thick_published",
                "has_measured_velocity",
            ]
        ],
        on="kepid",
        how="inner",
        validate="many_to_one",
    )
    selected["disk"] = selected["disk_published"]
    selected["P_thick"] = selected["p_thick_published"]
    selected["paper_category_target"] = [
        PAPER_COUNTS[(disk, system)]
        for disk, system in zip(selected["disk"], selected["system"])
    ]
    selected["pre_visual_category_count"] = selected.groupby(
        ["disk", "system"]
    )["kepoi_name"].transform("size")
    selected["category_excess_before_visual_qc"] = (
        selected["pre_visual_category_count"] - selected["paper_category_target"]
    )
    selected["paper_membership_status"] = np.where(
        selected["category_excess_before_visual_qc"] == 0,
        "category_count_exact_before_visual_qc",
        "candidate_in_category_with_unpublished_final_rejections",
    )
    return selected.sort_values(
        ["disk", "system", "koi_target", "koi_period", "kepoi_name"]
    ).reset_index(drop=True)


def add_stellar_metadata(
    inventory: pd.DataFrame,
    berger2020: pd.DataFrame,
    berger2018: pd.DataFrame | None,
) -> pd.DataFrame:
    out = inventory.merge(
        berger2020.drop_duplicates("kepid"),
        on="kepid",
        how="left",
        validate="many_to_one",
    )
    if berger2018 is not None:
        keep = [
            column
            for column in [
                "kepid",
                "berger2018_teff",
                "berger2018_rad",
                "berger2018_evol",
                "berger2018_bin",
            ]
            if column in berger2018.columns
        ]
        out = out.merge(
            berger2018[keep].drop_duplicates("kepid"),
            on="kepid",
            how="left",
            validate="many_to_one",
        )
    return out


def result_target_sets(archive_dir: Path, cloud_dir: Path) -> tuple[set[str], set[str]]:
    archive = {
        path.name.removesuffix("-results.fits")
        for path in archive_dir.glob("*-results.fits")
        if path.stat().st_size > 0
    }
    cloud = {
        path.parent.name
        for path in cloud_dir.glob("*/*-results.fits")
        if path.stat().st_size > 0
    }
    return archive, cloud


def add_result_coverage(
    inventory: pd.DataFrame, archive_targets: set[str], cloud_targets: set[str]
) -> pd.DataFrame:
    out = inventory.copy()
    out["archive_result_target_available"] = out["koi_target"].isin(archive_targets)
    out["cloud_result_target_available"] = out["koi_target"].isin(cloud_targets)
    out["any_result_target_available"] = (
        out["archive_result_target_available"]
        | out["cloud_result_target_available"]
    )
    out["preferred_result_source"] = np.select(
        [
            out["archive_result_target_available"],
            out["cloud_result_target_available"],
        ],
        ["original_archive", "cloud_missing_run"],
        default="missing",
    )
    return out


def count_audit(inventory: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for (disk, system), group in inventory.groupby(["disk", "system"]):
        paper = PAPER_COUNTS[(disk, system)]
        rows.append(
            {
                "disk": disk,
                "system": system,
                "pre_visual_planets": len(group),
                "paper_planets": paper,
                "implied_final_rejections": len(group) - paper,
                "hosts": group["kepid"].nunique(),
                "result_target_planets": int(group["any_result_target_available"].sum()),
                "result_target_hosts": int(
                    group.loc[group["any_result_target_available"], "kepid"].nunique()
                ),
                "missing_result_target_planets": int(
                    (~group["any_result_target_available"]).sum()
                ),
                "missing_result_target_hosts": int(
                    group.loc[~group["any_result_target_available"], "kepid"].nunique()
                ),
            }
        )
    total = {
        "disk": "all",
        "system": "all",
        "pre_visual_planets": len(inventory),
        "paper_planets": sum(PAPER_COUNTS.values()),
        "implied_final_rejections": len(inventory) - sum(PAPER_COUNTS.values()),
        "hosts": inventory["kepid"].nunique(),
        "result_target_planets": int(inventory["any_result_target_available"].sum()),
        "result_target_hosts": int(
            inventory.loc[inventory["any_result_target_available"], "kepid"].nunique()
        ),
        "missing_result_target_planets": int(
            (~inventory["any_result_target_available"]).sum()
        ),
        "missing_result_target_hosts": int(
            inventory.loc[~inventory["any_result_target_available"], "kepid"].nunique()
        ),
    }
    rows.append(total)
    return pd.DataFrame(rows)


def validate_count_contract(inventory: pd.DataFrame, audit: pd.DataFrame) -> None:
    observed = {
        (disk, system): len(group)
        for (disk, system), group in inventory.groupby(["disk", "system"])
    }
    if observed != PRE_VISUAL_COUNTS:
        raise ValueError(
            f"Published-host pre-visual count contract failed: "
            f"{observed} != {PRE_VISUAL_COUNTS}"
        )
    total = audit.loc[(audit["disk"] == "all") & (audit["system"] == "all")].iloc[0]
    if int(total["pre_visual_planets"]) != 2491:
        raise ValueError("Expected 2491 pre-visual-QC planets")
    if int(total["hosts"]) != 1888:
        raise ValueError("Expected all 1888 published hosts")
    if int(total["implied_final_rejections"]) != 26:
        raise ValueError("Expected exactly 26 implied final rejections")


def write_markdown(audit: pd.DataFrame, path: Path) -> None:
    lines = [
        "# Sagear published-host planet inventory reconstruction",
        "",
        "The exact published 1,888-host table was joined to the KOI catalog.",
        "Planets were retained when their disposition was not FALSE POSITIVE and",
        "their period was between 1 and 100 days. Multiplicity was assigned from",
        "the pre-cut KOI `koi_count`, not recounted after filtering.",
        "",
        "| disk | system | pre-visual planets | paper planets | implied removals | "
        "result-covered planets | missing-result planets |",
        "|---|---|---:|---:|---:|---:|---:|",
    ]
    for _, row in audit.iterrows():
        lines.append(
            f"| {row['disk']} | {row['system']} | "
            f"{int(row['pre_visual_planets'])} | {int(row['paper_planets'])} | "
            f"{int(row['implied_final_rejections'])} | "
            f"{int(row['result_target_planets'])} | "
            f"{int(row['missing_result_target_planets'])} |"
        )
    lines.extend(
        [
            "",
            "The two single-planet categories match the paper exactly before",
            "visual posterior QC. The remaining excess is exactly 21 thin-multi",
            "and 5 thick-multi planets, totaling 26. The public article does not",
            "identify those 26 planets, so individual membership remains unknown.",
            "The result-coverage columns are target-level availability checks;",
            "planet-level period matching must still be performed by the extractor.",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Reconstruct Sagear's planet inventory from the published hosts."
    )
    parser.add_argument("--config", default=None)
    parser.add_argument(
        "--archive-results-dir",
        default="inputs/alderaan_archive_raw_20260612/ALDERAAN_posteriors",
    )
    parser.add_argument(
        "--cloud-results-dir",
        default="alderaan_project/Results/sagear_missing",
    )
    args = parser.parse_args()

    cfg = load_config(args.config)
    hosts = read_sagear2026_kinematic_hosts(cfg)
    inventory = build_pre_visual_inventory(read_koi(cfg), hosts)
    berger2018_path = root_path(cfg, "berger2018_stellar")
    berger2018 = (
        read_berger2018_stellar_table(berger2018_path)
        if berger2018_path is not None and berger2018_path.exists()
        else None
    )
    inventory = add_stellar_metadata(
        inventory, read_berger_table2(cfg), berger2018
    )

    base = Path(__file__).resolve().parent
    archive_dir = Path(args.archive_results_dir)
    cloud_dir = Path(args.cloud_results_dir)
    if not archive_dir.is_absolute():
        archive_dir = base / archive_dir
    if not cloud_dir.is_absolute():
        cloud_dir = base / cloud_dir
    archive_targets, cloud_targets = result_target_sets(archive_dir, cloud_dir)
    inventory = add_result_coverage(inventory, archive_targets, cloud_targets)
    audit = count_audit(inventory)
    validate_count_contract(inventory, audit)

    out = output_dir()
    inventory.to_csv(
        out / "sagear2026_planet_inventory_pre_visual_qc.csv", index=False
    )
    audit.to_csv(
        out / "sagear2026_planet_inventory_count_audit.csv", index=False
    )
    inventory.loc[~inventory["any_result_target_available"]].to_csv(
        out / "sagear2026_planets_missing_alderaan_target_results.csv",
        index=False,
    )
    write_markdown(
        audit, out / "sagear2026_planet_inventory_reconstruction.md"
    )
    print(audit.to_string(index=False))


if __name__ == "__main__":
    main()
