from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd

from common import load_config, output_dir


def main() -> None:
    parser = argparse.ArgumentParser(description="Audit Sagear manuscript macro count consistency.")
    parser.add_argument("--config", default=None)
    parser.add_argument("--main-tex", default=None)
    parser.add_argument("--out-prefix", default="target_consistency_diagnostics")
    args = parser.parse_args()

    cfg = load_config(args.config)
    root = Path(cfg["_root"])
    tex_path = Path(args.main_tex) if args.main_tex else root / "main.tex"

    macros = parse_macros(tex_path)
    macro_table = pd.DataFrame([{"macro": k, "value": v} for k, v in sorted(macros.items())])
    checks = consistency_checks(macros)

    out = output_dir()
    macro_path = out / f"{args.out_prefix}_macros.csv"
    checks_path = out / f"{args.out_prefix}_checks.csv"
    macro_table.to_csv(macro_path, index=False)
    checks.to_csv(checks_path, index=False)

    print("=== Sagear Macro Counts ===")
    print(macro_table.to_string(index=False))
    print("\n=== Consistency Checks ===")
    print(checks.to_string(index=False))
    print(f"\nWrote: {macro_path}")
    print(f"Wrote: {checks_path}")


def parse_macros(tex_path: Path) -> dict[str, int]:
    text = tex_path.read_text(errors="replace")
    return {
        match.group(1): int(match.group(2))
        for match in re.finditer(r"\\newcommand\\(all\w+)\{(\d+)\s*\}", text)
    }


def consistency_checks(macros: dict[str, int]) -> pd.DataFrame:
    rows = []

    def add(name: str, left: int, right: int, note: str) -> None:
        rows.append(
            {
                "check": name,
                "left_value": left,
                "right_value": right,
                "delta_left_minus_right": left - right,
                "consistent": left == right,
                "note": note,
            }
        )

    thin_planets = macros.get("allthinplanets", 0)
    thick_planets = macros.get("allthickplanets", 0)
    subgroup_planets = (
        macros.get("allthinsingles", 0)
        + macros.get("allthicksingles", 0)
        + macros.get("allthinmultiplanets", 0)
        + macros.get("allthickmultiplanets", 0)
    )
    thin_subgroups = macros.get("allthinsingles", 0) + macros.get("allthinmultiplanets", 0)
    thick_subgroups = macros.get("allthicksingles", 0) + macros.get("allthickmultiplanets", 0)
    thin_stars = macros.get("allthinstars", 0)
    thick_stars = macros.get("allthickstars", 0)
    subgroup_hosts_like = (
        macros.get("allthinsingles", 0)
        + macros.get("allthicksingles", 0)
        + macros.get("allthinmultistars", 0)
        + macros.get("allthickmultistars", 0)
    )

    add("disk_planets_vs_allplanets", thin_planets + thick_planets, macros.get("allplanets", 0), "allthinplanets + allthickplanets vs allplanets")
    add("subgroup_planets_vs_allplanets", subgroup_planets, macros.get("allplanets", 0), "singles + multi planet subgroup macros vs allplanets")
    add("thin_subgroups_vs_thin_planets", thin_subgroups, thin_planets, "thin singles + thin multi planets vs allthinplanets")
    add("thick_subgroups_vs_thick_planets", thick_subgroups, thick_planets, "thick singles + thick multi planets vs allthickplanets")
    add("disk_stars_vs_allstars", thin_stars + thick_stars, macros.get("allstars", 0), "allthinstars + allthickstars vs allstars")
    add("subgroup_hosts_like_vs_allstars", subgroup_hosts_like, macros.get("allstars", 0), "single hosts + multi hosts from subgroup macros vs allstars")
    add("thin_host_subgroups_vs_thin_stars", macros.get("allthinsingles", 0) + macros.get("allthinmultistars", 0), thin_stars, "thin single hosts + thin multi hosts vs allthinstars")
    add("thick_host_subgroups_vs_thick_stars", macros.get("allthicksingles", 0) + macros.get("allthickmultistars", 0), thick_stars, "thick single hosts + thick multi hosts vs allthickstars")

    return pd.DataFrame(rows)


if __name__ == "__main__":
    main()
