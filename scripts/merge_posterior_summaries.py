from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from common import output_dir


def main() -> None:
    parser = argparse.ArgumentParser(description="Merge existing and newly extracted ALDERAAN e,omega posterior summaries.")
    parser.add_argument(
        "--base",
        default=None,
        help="Existing posterior summary CSV. Defaults to eccentricity_posterior_summary_old_astropy_rawcc.csv.",
    )
    parser.add_argument("--new", required=True, help="New posterior summary CSV, e.g. cloud extraction output.")
    parser.add_argument("--out", default=None, help="Merged output CSV path.")
    parser.add_argument("--coverage-out", default=None, help="Merged coverage output CSV path.")
    parser.add_argument("--sample", default=None, help="Canonical sample path for coverage.")
    args = parser.parse_args()

    out_dir = output_dir()
    base_path = Path(args.base) if args.base else out_dir / "eccentricity_posterior_summary_old_astropy_rawcc.csv"
    new_path = Path(args.new)
    out_path = Path(args.out) if args.out else out_dir / "eccentricity_posterior_summary_merged.csv"
    coverage_path = Path(args.coverage_out) if args.coverage_out else out_dir / "eccentricity_posterior_coverage_merged.csv"
    sample_path = Path(args.sample) if args.sample else out_dir / "canonical_sample_old_astropy_rawcc.csv"

    base = pd.read_csv(base_path)
    new = pd.read_csv(new_path)
    sample = pd.read_csv(sample_path)

    base["posterior_source"] = "existing_archive"
    new["posterior_source"] = "new_alderaan"
    merged = pd.concat([base, new], ignore_index=True, sort=False)
    merged["_source_rank"] = merged["posterior_source"].map({"existing_archive": 0, "new_alderaan": 1}).fillna(0)
    merged = (
        merged.sort_values(["kepoi_name", "_source_rank"])
        .drop_duplicates("kepoi_name", keep="last")
        .drop(columns=["_source_rank"])
        .sort_values(["disk", "system", "koi_target", "koi_period"])
        .reset_index(drop=True)
    )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    coverage_path.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(out_path, index=False)
    coverage = coverage_summary(sample, merged)
    coverage.to_csv(coverage_path, index=False)

    print(f"Base rows: {len(base)}")
    print(f"New rows: {len(new)}")
    print(f"Merged unique posterior rows: {len(merged)}")
    print(f"Wrote: {out_path}")
    print(f"Wrote: {coverage_path}")
    print("\nCoverage:")
    print(coverage.to_string(index=False))


def coverage_summary(sample: pd.DataFrame, summary: pd.DataFrame) -> pd.DataFrame:
    rows = []
    summary_ids = set(summary["kepoi_name"])
    for (disk, system), sub in sample.groupby(["disk", "system"]):
        posterior_planets = int(sub["kepoi_name"].isin(summary_ids).sum())
        rows.append(
            {
                "disk": disk,
                "system": system,
                "sample_planets": int(len(sub)),
                "posterior_planets": posterior_planets,
                "missing_planets": int(len(sub) - posterior_planets),
                "coverage_fraction": posterior_planets / len(sub) if len(sub) else 0.0,
            }
        )
    return pd.DataFrame(rows).sort_values(["disk", "system"])


if __name__ == "__main__":
    main()
