#!/usr/bin/env python
"""Emit validated KOI target/KIC pairs from a recovery target CSV."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def load_pairs(path: Path) -> list[tuple[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"koi_target", "kepid"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"{path} must contain koi_target and kepid columns")
        pairs = [
            (row["koi_target"].strip(), row["kepid"].strip())
            for row in reader
        ]
    if not pairs or any(not target or not kepid for target, kepid in pairs):
        raise ValueError(f"{path} contains an empty target or kepid")
    if len({target for target, _ in pairs}) != len(pairs):
        raise ValueError(f"{path} contains duplicate koi_target rows")
    return pairs


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("targets", type=Path)
    args = parser.parse_args()
    for target, kepid in load_pairs(args.targets):
        print(f"{target},{kepid}")


if __name__ == "__main__":
    main()
