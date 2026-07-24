"""Validate that a posterior summary belongs to the intended Sagear inventory.

Older ALDERAAN summaries in this project were assembled from overlapping archive
runs.  This audit makes that failure visible by comparing planet identifiers,
host identifiers, duplicates, and population counts against one explicit
inventory.  It never silently trims or repairs a summary.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


REQUIRED = {"kepoi_name", "koi_target", "disk", "system"}


def _require(frame: pd.DataFrame, label: str) -> None:
    missing = sorted(REQUIRED - set(frame.columns))
    if missing:
        raise ValueError(f"{label} is missing required columns: {missing}")


def audit_inventory(inventory: pd.DataFrame, summary: pd.DataFrame) -> dict[str, object]:
    _require(inventory, "inventory")
    _require(summary, "summary")
    inv_planets = set(inventory["kepoi_name"].astype(str))
    sum_planets = set(summary["kepoi_name"].astype(str))
    inv_hosts = set(inventory["koi_target"].astype(str))
    sum_hosts = set(summary["koi_target"].astype(str))

    inv_groups = (
        inventory.assign(kepoi_name=inventory["kepoi_name"].astype(str))
        .groupby(["disk", "system"])
        .size()
        .astype(int)
        .to_dict()
    )
    sum_groups = (
        summary.assign(kepoi_name=summary["kepoi_name"].astype(str))
        .groupby(["disk", "system"])
        .size()
        .astype(int)
        .to_dict()
    )
    return {
        "inventory_planets": len(inv_planets),
        "summary_planets": len(sum_planets),
        "inventory_hosts": len(inv_hosts),
        "summary_hosts": len(sum_hosts),
        "inventory_duplicate_planet_rows": int(inventory["kepoi_name"].duplicated().sum()),
        "summary_duplicate_planet_rows": int(summary["kepoi_name"].duplicated().sum()),
        "summary_planets_inside_inventory": len(sum_planets & inv_planets),
        "summary_planets_outside_inventory": len(sum_planets - inv_planets),
        "inventory_planets_without_summary": len(inv_planets - sum_planets),
        "summary_hosts_inside_inventory": len(sum_hosts & inv_hosts),
        "summary_hosts_outside_inventory": len(sum_hosts - inv_hosts),
        "inventory_groups": {f"{k[0]}_{k[1]}": v for k, v in inv_groups.items()},
        "summary_groups": {f"{k[0]}_{k[1]}": v for k, v in sum_groups.items()},
        "summary_posterior_sources": (
            summary["posterior_source"].astype(str).value_counts().to_dict()
            if "posterior_source" in summary
            else {}
        ),
        "summary_transit_fit_sources": (
            summary["transit_fit_source"].astype(str).value_counts().to_dict()
            if "transit_fit_source" in summary
            else {}
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inventory", default="outputs/sagear2026_planet_inventory_pre_visual_qc.csv")
    parser.add_argument("--summary", required=True)
    parser.add_argument("--output", default=None)
    parser.add_argument("--strict", action="store_true", help="fail if any summary planet is outside inventory or duplicated")
    args = parser.parse_args()

    inventory = pd.read_csv(args.inventory)
    summary = pd.read_csv(args.summary)
    result = audit_inventory(inventory, summary)
    result["inventory"] = str(Path(args.inventory).resolve())
    result["summary"] = str(Path(args.summary).resolve())
    result["strict_pass"] = (
        result["summary_planets_outside_inventory"] == 0
        and result["summary_duplicate_planet_rows"] == 0
    )
    output = Path(args.output) if args.output else Path(args.summary).with_suffix(".inventory_audit.json")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, default=int) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, default=int))
    if args.strict and not result["strict_pass"]:
        raise SystemExit("summary inventory contract failed")


if __name__ == "__main__":
    main()
