from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from common import load_config, output_dir, read_koi


TARGETS = {
    "thin_singles_planets": 1121,
    "thick_singles_planets": 275,
    "thin_multi_planets": 862,
    "thick_multi_planets": 207,
    "thin_multi_hosts": 394,
    "thick_multi_hosts": 98,
}


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Audit whether Sagear single/multi labels are defined before or after quality cuts."
    )
    parser.add_argument("--config", default=None)
    parser.add_argument("--sample", default=None)
    parser.add_argument("--old-astropy-sample", default=None)
    args = parser.parse_args()

    cfg = load_config(args.config)
    out = output_dir()
    sample_path = Path(args.sample) if args.sample else out / "canonical_sample_diagnostic.csv"
    old_path = Path(args.old_astropy_sample) if args.old_astropy_sample else out / "canonical_sample_diagnostic_old_astropy.csv"

    raw = read_koi(cfg)
    direct = pd.read_csv(sample_path)
    samples = [("direct", direct)]
    if old_path.exists():
        samples.append(("old_astropy", pd.read_csv(old_path)))

    definitions = build_system_definitions(raw, direct)
    rows = []
    moved_rows = []
    for label, sample in samples:
        for def_name, system_by_kepid in definitions.items():
            assigned = assign_system(sample, def_name, system_by_kepid)
            rows.append(count_row(label, def_name, sample, assigned))
            moved_rows.extend(moved_planets(label, def_name, sample, assigned))

    counts = pd.DataFrame(rows).sort_values(["l1_error", "classifier", "system_definition"])
    moved = pd.DataFrame(moved_rows)
    counts_path = out / "system_definition_counts.csv"
    moved_path = out / "system_definition_moved_planets.csv"
    counts.to_csv(counts_path, index=False)
    moved.to_csv(moved_path, index=False)

    print(counts.to_string(index=False))
    print(f"\nWrote: {counts_path}")
    print(f"Wrote: {moved_path}")


def build_system_definitions(raw: pd.DataFrame, sample: pd.DataFrame) -> dict[str, pd.Series | None]:
    definitions: dict[str, pd.Series | None] = {"after_all_cuts_current": None}

    koi_count = pd.to_numeric(sample.drop_duplicates("kepid").set_index("kepid")["koi_count"], errors="coerce")
    definitions["koi_count_gt1_catalog"] = pd.Series(np.where(koi_count > 1, "multi", "single"), index=koi_count.index)

    disposition = raw["koi_disposition"].isin(["CONFIRMED", "CANDIDATE"])
    period = (raw["koi_period"] >= 1.0) & (raw["koi_period"] <= 100.0)
    confirmed = raw["koi_disposition"].eq("CONFIRMED")

    definitions["raw_confirmed_candidate_count_gt1"] = count_to_system(raw[disposition])
    definitions["period_confirmed_candidate_count_gt1"] = count_to_system(raw[disposition & period])
    definitions["raw_confirmed_only_count_gt1"] = count_to_system(raw[confirmed])
    definitions["period_confirmed_only_count_gt1"] = count_to_system(raw[confirmed & period])

    return definitions


def count_to_system(df: pd.DataFrame) -> pd.Series:
    counts = df.groupby("kepid")["kepoi_name"].count()
    return pd.Series(np.where(counts > 1, "multi", "single"), index=counts.index)


def assign_system(sample: pd.DataFrame, def_name: str, system_by_kepid: pd.Series | None) -> pd.Series:
    if system_by_kepid is None:
        return sample["system"].astype(str)
    out = sample["kepid"].map(system_by_kepid)
    if out.isna().any():
        missing = sorted(sample.loc[out.isna(), "kepid"].astype(int).unique())
        raise ValueError(
            f"System definition {def_name} lacks multiplicity for KIC identifiers: {missing[:10]}"
        )
    return out.astype(str)


def count_row(classifier: str, definition: str, sample: pd.DataFrame, system: pd.Series) -> dict[str, object]:
    row: dict[str, object] = {"classifier": classifier, "system_definition": definition}
    for disk in ["thin", "thick"]:
        for sys_label, macro in [("single", "singles"), ("multi", "multi")]:
            sub = sample[(sample["disk"] == disk) & (system == sys_label)]
            if sys_label == "single":
                key = f"{disk}_singles_planets"
            else:
                key = f"{disk}_multi_planets"
            row[key] = int(len(sub))
            row[f"delta_{key}"] = int(len(sub)) - TARGETS[key]
            if sys_label == "multi":
                hkey = f"{disk}_multi_hosts"
                row[hkey] = int(sub["kepid"].nunique())
                row[f"delta_{hkey}"] = int(sub["kepid"].nunique()) - TARGETS[hkey]
    keys = [
        "thin_singles_planets",
        "thick_singles_planets",
        "thin_multi_planets",
        "thick_multi_planets",
        "thin_multi_hosts",
        "thick_multi_hosts",
    ]
    row["l1_error"] = int(sum(abs(row[f"delta_{k}"]) for k in keys))
    row["moved_from_current_planets"] = int((system != sample["system"].astype(str)).sum())
    return row


def moved_planets(classifier: str, definition: str, sample: pd.DataFrame, system: pd.Series) -> list[dict[str, object]]:
    changed = sample[system != sample["system"].astype(str)].copy()
    if changed.empty:
        return []
    changed["new_system"] = system.loc[changed.index].to_numpy()
    keep = [
        "kepid",
        "kepoi_name",
        "disk",
        "system",
        "new_system",
        "koi_count",
        "koi_period",
        "koi_model_snr",
        "koi_prad",
        "P_thick",
    ]
    rows = []
    for row in changed[[c for c in keep if c in changed.columns]].to_dict("records"):
        row["classifier"] = classifier
        row["system_definition"] = definition
        rows.append(row)
    return rows


if __name__ == "__main__":
    main()
