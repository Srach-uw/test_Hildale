from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


PAPER_COUNTS = {
    ("thin", "single"): 1121,
    ("thick", "single"): 275,
    ("thin", "multi"): 862,
    ("thick", "multi"): 207,
}


def _numeric(frame: pd.DataFrame, column: str) -> pd.Series:
    if column not in frame:
        return pd.Series(np.nan, index=frame.index, dtype=float)
    return pd.to_numeric(frame[column], errors="coerce")


def _boolean(frame: pd.DataFrame, column: str) -> pd.Series:
    if column not in frame:
        return pd.Series(False, index=frame.index, dtype=bool)
    values = frame[column]
    if values.dtype == bool:
        return values.fillna(False)
    return (
        values.astype("string")
        .str.strip()
        .str.lower()
        .isin({"true", "1", "yes"})
    )


def validate_contract(frame: pd.DataFrame) -> pd.DataFrame:
    required = {"kepid", "kepoi_name", "disk", "system", "posterior_status"}
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"Posterior manifest is missing required columns: {missing}")
    if frame["kepoi_name"].duplicated().any():
        raise ValueError("Posterior manifest contains duplicate planet identifiers")

    observed = frame.groupby(["disk", "system"]).size()
    rows = []
    for key, paper_count in PAPER_COUNTS.items():
        pre_visual = int(observed.get(key, 0))
        rows.append(
            {
                "disk": key[0],
                "system": key[1],
                "pre_visual_planets": pre_visual,
                "paper_planets": paper_count,
                "implied_visual_or_convergence_rejections": pre_visual - paper_count,
            }
        )
    audit = pd.DataFrame(rows)
    if int(audit["pre_visual_planets"].sum()) != 2491:
        raise ValueError("Expected 2,491 reconstructed pre-visual-QC planets")
    if int(audit["paper_planets"].sum()) != 2465:
        raise ValueError("Expected 2,465 published final planets")
    if int(audit["implied_visual_or_convergence_rejections"].sum()) != 26:
        raise ValueError("Expected exactly 26 count-implied final rejections")
    if int(frame["kepid"].nunique()) != 1888:
        raise ValueError("Expected the 1,888 published host stars")
    return audit


def rank_unpublished_rejection_candidates(
    frame: pd.DataFrame, audit: pd.DataFrame
) -> pd.DataFrame:
    candidates = frame.loc[frame["system"].eq("multi")].copy()
    missing_posterior = ~candidates["posterior_status"].eq("posterior_available")
    primary_exclude = _boolean(candidates, "qc_primary_exclude")
    impact = _numeric(candidates, "koi_impact")
    snr = _numeric(candidates, "koi_model_snr")
    transits = _numeric(candidates, "koi_num_transits")

    # The score is a reproducible diagnostic ranking, not recovered authorship data.
    candidates["candidate_score"] = (
        1000.0 * missing_posterior.astype(float)
        + 500.0 * primary_exclude.astype(float)
        + 10.0 * impact.fillna(0).clip(lower=0, upper=2)
        + 20.0 / np.sqrt(snr.clip(lower=1).fillna(1))
        + 10.0 / np.sqrt(transits.clip(lower=1).fillna(1))
    )
    candidates["missing_posterior"] = missing_posterior
    candidates["primary_qc_exclude"] = primary_exclude
    candidates["candidate_rank_within_disk"] = (
        candidates.groupby("disk")["candidate_score"]
        .rank(method="first", ascending=False)
        .astype(int)
    )
    removal_targets = {
        row["disk"]: int(row["implied_visual_or_convergence_rejections"])
        for _, row in audit.loc[audit["system"].eq("multi")].iterrows()
    }
    candidates["count_constrained_candidate_set"] = [
        rank <= removal_targets[disk]
        for disk, rank in zip(
            candidates["disk"], candidates["candidate_rank_within_disk"]
        )
    ]
    candidates["identity_status"] = np.where(
        candidates["count_constrained_candidate_set"],
        "sensitivity_candidate_not_author_confirmed",
        "not_selected_by_diagnostic_ranking",
    )
    keep = [
        "kepid",
        "kepoi_name",
        "koi_target",
        "disk",
        "system",
        "posterior_status",
        "qc_primary_exclude",
        "qc_reasons",
        "koi_period",
        "koi_impact",
        "koi_model_snr",
        "koi_num_transits",
        "missing_posterior",
        "primary_qc_exclude",
        "candidate_score",
        "candidate_rank_within_disk",
        "count_constrained_candidate_set",
        "identity_status",
    ]
    return candidates[[column for column in keep if column in candidates]].sort_values(
        ["disk", "candidate_rank_within_disk", "kepoi_name"]
    )


def write_report(audit: pd.DataFrame, candidates: pd.DataFrame, path: Path) -> None:
    selected = candidates.loc[candidates["count_constrained_candidate_set"]]
    lines = [
        "# Final visual-QC contract audit",
        "",
        "## Confirmed from the published article",
        "",
        "- The public Table 1 is a host-level kinematic table with 1,888 stars.",
        "- Transit fits are visually inspected after ALDERAAN fitting.",
        "- Less than 2% are removed for quality or non-convergent posteriors.",
        "- The final sample contains 2,465 planets.",
        "",
        "## Reconstructed arithmetic",
        "",
        "| disk | system | before visual QC | published | implied removals |",
        "|---|---|---:|---:|---:|",
    ]
    for _, row in audit.iterrows():
        lines.append(
            f"| {row['disk']} | {row['system']} | "
            f"{int(row['pre_visual_planets'])} | {int(row['paper_planets'])} | "
            f"{int(row['implied_visual_or_convergence_rejections'])} |"
        )
    lines.extend(
        [
            "",
            "The difference is 26/2,491 = 1.04%, consistent with the stated",
            "\"less than 2%\" rejection. The count arithmetic implies no removals",
            "from either singles category, 5 from thick multis, and 21 from thin",
            "multis, subject to the KOI catalog epoch used in this reconstruction.",
            "",
            "## What remains unknown",
            "",
            "The article does not publish the 26 planet identifiers. Table 1 cannot",
            "encode this planet-level decision because it contains one row per host.",
            "The candidate ledger ranks plausible removals from missing posterior",
            "products, explicit QC flags, impact parameter, signal-to-noise, and",
            "transit count. Its selected rows are sensitivity candidates only and",
            "must not be described as Sagear's actual rejected planets.",
            "",
            f"Diagnostic candidate set: {len(selected)} planets "
            f"({int((selected['disk'] == 'thick').sum())} thick, "
            f"{int((selected['disk'] == 'thin').sum())} thin).",
            "",
            "## Operational decision",
            "",
            "Fit every recoverable missing system. Do not pre-delete 26 planets",
            "before fitting. After extraction, report the complete pre-visual branch,",
            "the deterministic count-constrained sensitivity branch, and any",
            "author-confirmed branch separately.",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    manifest = pd.read_csv(args.manifest, low_memory=False)
    output = Path(args.output_dir).resolve()
    output.mkdir(parents=True, exist_ok=True)
    audit = validate_contract(manifest)
    candidates = rank_unpublished_rejection_candidates(manifest, audit)

    audit.to_csv(output / "visual_qc_count_contract.csv", index=False)
    candidates.to_csv(output / "visual_qc_candidate_ledger.csv", index=False)
    write_report(audit, candidates, output / "visual_qc_contract_audit.md")
    summary = {
        "pre_visual_planets": int(audit["pre_visual_planets"].sum()),
        "published_planets": int(audit["paper_planets"].sum()),
        "implied_rejections": int(
            audit["implied_visual_or_convergence_rejections"].sum()
        ),
        "implied_rejection_fraction": float(
            audit["implied_visual_or_convergence_rejections"].sum()
            / audit["pre_visual_planets"].sum()
        ),
        "candidate_identities_author_confirmed": False,
    }
    (output / "visual_qc_contract_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
