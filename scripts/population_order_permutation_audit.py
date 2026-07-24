from __future__ import annotations

import argparse
import itertools
from pathlib import Path

import numpy as np
import pandas as pd


PUBLISHED = {
    "thick_singles": (0.066, 0.045, 0.096),
    "thin_singles": (0.022, 0.017, 0.029),
    "thick_multis": (0.033, 0.015, 0.065),
    "thin_multis": (0.030, 0.023, 0.031),
}
PUBLISHED_ORDER = list(PUBLISHED)
REQUIRED = ["population", "expected_e", "expected_e_lo", "expected_e_hi"]


def load_fit(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path)
    missing = [column for column in REQUIRED if column not in frame.columns]
    if missing:
        raise ValueError(f"{path} is missing columns: {missing}")
    if set(frame["population"]) != set(PUBLISHED):
        raise ValueError(f"Unexpected population labels in {path}")
    return frame.set_index("population").loc[PUBLISHED_ORDER].reset_index()


def score_permutation(frame: pd.DataFrame, source_order: tuple[str, ...]) -> dict[str, object]:
    source = frame.set_index("population").loc[list(source_order)].reset_index()
    values = source["expected_e"].to_numpy(float)
    lows = source["expected_e_lo"].to_numpy(float)
    highs = source["expected_e_hi"].to_numpy(float)
    paper = np.array([PUBLISHED[name][0] for name in PUBLISHED_ORDER])
    paper_lo = np.array([PUBLISHED[name][1] for name in PUBLISHED_ORDER])
    paper_hi = np.array([PUBLISHED[name][2] for name in PUBLISHED_ORDER])
    overlap = (highs >= paper_lo) & (lows <= paper_hi)
    return {
        "source_order": "|".join(source_order),
        "target_order": "|".join(PUBLISHED_ORDER),
        "interval_overlaps": int(overlap.sum()),
        "rmse": float(np.sqrt(np.mean((values - paper) ** 2))),
        "mean_absolute_error": float(np.mean(np.abs(values - paper))),
        "max_absolute_error": float(np.max(np.abs(values - paper))),
        "thick_singles_value": values[0],
        "thin_singles_value": values[1],
        "thick_multis_value": values[2],
        "thin_multis_value": values[3],
        "overlap_flags": "|".join("yes" if value else "no" for value in overlap),
    }


def audit(frame: pd.DataFrame) -> pd.DataFrame:
    rows = [score_permutation(frame, permutation) for permutation in itertools.permutations(PUBLISHED_ORDER)]
    return (
        pd.DataFrame(rows)
        .sort_values(["rmse", "mean_absolute_error", "source_order"])
        .reset_index(drop=True)
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Enumerate hidden population-array order permutations.")
    parser.add_argument("--outputs", default="outputs")
    parser.add_argument("--equal", default=None)
    parser.add_argument("--dynesty", default=None)
    parser.add_argument("--tag", default="", help="Optional suffix for output artifacts.")
    args = parser.parse_args()
    root = Path(args.outputs).expanduser().resolve()
    equal_path = Path(args.equal) if args.equal else root / (
        "rayleigh_population_fit_transit_selection_manuscript_reciprocal_"
        "exact_equal_raw_diagnostic_20260724.csv"
    )
    dynesty_path = Path(args.dynesty) if args.dynesty else root / (
        "rayleigh_population_fit_transit_selection_manuscript_reciprocal_"
        "exact_dynesty_weighted_20260724.csv"
    )
    outputs = []
    for method, path in [("equal_raw", equal_path), ("dynesty_weighted", dynesty_path)]:
        table = audit(load_fit(path))
        table.insert(0, "method", method)
        outputs.append(table)
    result = pd.concat(outputs, ignore_index=True)
    suffix = f"_{args.tag}" if args.tag else ""
    result.to_csv(root / f"population_order_permutation_audit{suffix}.csv", index=False)

    best = result.groupby("method", sort=False).head(1)
    lines = [
        "# Population-array order audit",
        "",
        "This audit enumerates all 24 assignments of the four locally fitted",
        "population results to the four published Table 3 categories. It is a",
        "diagnostic for an undocumented array-order error, not a license to relabel",
        "the physical disk classifications.",
        "",
        "| method | best source order assigned to published order | overlaps / 4 | RMSE |",
        "|---|---|---:|---:|",
    ]
    for row in best.itertuples(index=False):
        lines.append(
            f"| {row.method} | `{row.source_order}` | {row.interval_overlaps} / 4 | {row.rmse:.5f} |"
        )
    lines.extend(
        [
            "",
            "The published target order is `thick_singles | thin_singles |",
            "thick_multis | thin_multis`. A best permutation that differs from the",
            "identity indicates a possible downstream result-array assignment issue,",
            "but cannot distinguish an author-side bug from an ordering mistake in",
            "this reconstruction without the original population code.",
        ]
    )
    (root / f"population_order_permutation_audit{suffix}.md").write_text(
        "\n".join(lines) + "\n", encoding="ascii"
    )
    print(best[["method", "source_order", "interval_overlaps", "rmse"]].to_string(index=False))


if __name__ == "__main__":
    main()
