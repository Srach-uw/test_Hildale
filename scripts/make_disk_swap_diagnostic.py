from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def swap_disk_labels(frame: pd.DataFrame) -> pd.DataFrame:
    required = {"disk", "system"}
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"summary is missing required columns: {missing}")

    labels = set(frame["disk"].dropna().astype(str))
    if labels != {"thin", "thick"}:
        raise ValueError(f"expected only thin/thick disk labels, found {sorted(labels)}")

    out = frame.copy()
    out["disk_published_unswapped"] = out["disk"]
    out["disk"] = out["disk"].map({"thin": "thick", "thick": "thin"})
    out["disk_assignment_mode"] = "inverted_for_diagnostic_only"
    return out


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Create an explicitly fenced thin/thick label-swap diagnostic. "
            "This output is not a canonical scientific sample."
        )
    )
    parser.add_argument("--summary", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    source = Path(args.summary).expanduser().resolve()
    destination = Path(args.out).expanduser().resolve()
    if source == destination:
        raise ValueError("diagnostic output must not overwrite the source summary")

    swapped = swap_disk_labels(pd.read_csv(source))
    destination.parent.mkdir(parents=True, exist_ok=True)
    swapped.to_csv(destination, index=False)
    print(f"Wrote {len(swapped)} explicitly diagnostic rows: {destination}")


if __name__ == "__main__":
    main()
