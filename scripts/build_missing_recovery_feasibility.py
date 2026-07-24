"""Summarize the recoverability of the exact published-inventory gap."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


POPULATION_ORDER = ["thick_multis", "thick_singles", "thin_multis", "thin_singles"]


def recovery_class(row: pd.Series) -> str:
    attempted = row["prior_launch_status"] == "previously_launched_no_usable_result"
    low_snr = float(row["max_snr"]) < 20.0
    if not attempted and not low_snr:
        return "A_never_attempted_snr_ge_20"
    if not attempted and low_snr:
        return "B_never_attempted_snr_lt_20"
    if attempted and not low_snr:
        return "C_prior_failed_snr_ge_20"
    return "D_prior_failed_snr_lt_20"


def build(targets: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    required = {
        "koi_target",
        "disk",
        "system",
        "missing_planets",
        "max_snr",
        "prior_launch_status",
    }
    missing = required - set(targets.columns)
    if missing:
        raise ValueError(f"Missing required target columns: {sorted(missing)}")

    detail = targets.copy()
    detail["population"] = detail["disk"] + "_" + detail["system"].replace(
        {"multi": "multis", "single": "singles"}
    )
    detail["recovery_class"] = detail.apply(recovery_class, axis=1)
    detail["canonical_action"] = detail["recovery_class"].map(
        {
            "A_never_attempted_snr_ge_20": "canonical_run",
            "B_never_attempted_snr_lt_20": "canonical_run_then_visual_qc",
            "C_prior_failed_snr_ge_20": "canonical_run_then_seeded_retry",
            "D_prior_failed_snr_lt_20": "canonical_run_then_seeded_retry_and_visual_qc",
        }
    )

    summary = (
        detail.groupby(["population", "recovery_class", "canonical_action"], observed=True)
        .agg(systems=("koi_target", "nunique"), planets=("missing_planets", "sum"))
        .reset_index()
    )
    summary["population"] = pd.Categorical(
        summary["population"], categories=POPULATION_ORDER, ordered=True
    )
    summary = summary.sort_values(["population", "recovery_class"]).reset_index(drop=True)
    return detail, summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--targets",
        type=Path,
        default=Path("cloud_published_inventory_missing_batch/targets_missing_launchable.csv"),
    )
    parser.add_argument(
        "--output-detail",
        type=Path,
        default=Path("outputs/missing_recovery_feasibility_targets.csv"),
    )
    parser.add_argument(
        "--output-summary",
        type=Path,
        default=Path("outputs/missing_recovery_feasibility_summary.csv"),
    )
    parser.add_argument(
        "--output-report",
        type=Path,
        default=Path("outputs/missing_recovery_feasibility.md"),
    )
    args = parser.parse_args()

    detail, summary = build(pd.read_csv(args.targets))
    for path in (args.output_detail, args.output_summary, args.output_report):
        path.parent.mkdir(parents=True, exist_ok=True)
    detail.to_csv(args.output_detail, index=False)
    summary.to_csv(args.output_summary, index=False)

    attempted = detail["prior_launch_status"].eq("previously_launched_no_usable_result")
    lines = [
        "# Missing-posterior recovery feasibility",
        "",
        f"- Exact missing systems: {detail.koi_target.nunique()}",
        f"- Exact missing planets: {int(detail.missing_planets.sum())}",
        f"- Never previously attempted: {int((~attempted).sum())} systems",
        f"- Previously attempted without usable FITS: {int(attempted.sum())} systems",
        f"- Targets with maximum catalog S/N below 20: {int((detail.max_snr < 20).sum())} systems",
        "",
        "S/N is a triage flag, not an exclusion rule. Every target receives the canonical run first.",
        "Only targets still lacking a nonempty result FITS receive the deterministic changed-seed retry.",
        "",
        summary.to_markdown(index=False),
        "",
    ]
    args.output_report.write_text("\n".join(lines), encoding="utf-8")
    print(summary.to_string(index=False))
    print(f"systems,{detail.koi_target.nunique()}")
    print(f"planets,{int(detail.missing_planets.sum())}")


if __name__ == "__main__":
    main()
