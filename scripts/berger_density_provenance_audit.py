"""Audit the Berger stellar-density provenance used by the Sagear replication.

This is deliberately an audit, not a silent replacement of the canonical
density prior. Sagear's published methods prose names Berger et al. (2018),
while the public source comments and conclusion point to Berger et al. (2020).
The published materials do not resolve that provenance conflict, so the two
sources are reported separately.
"""

from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd


def parse_b20(path: Path) -> pd.DataFrame:
    rows = []
    with gzip.open(path, "rt") as handle:
        for line in handle:
            try:
                rows.append(
                    {
                        "kepid": int(line[0:8]),
                        "mass": float(line[9:14]),
                        "logg": float(line[51:57]),
                        "radius": float(line[91:98]),
                        "rho_log": float(line[116:122]),
                        "rho_log_upper": float(line[123:129]),
                        "rho_log_lower": float(line[130:136]),
                    }
                )
            except (TypeError, ValueError):
                continue
    frame = pd.DataFrame(rows).drop_duplicates("kepid")
    if not frame.empty:
        frame["kepid"] = frame["kepid"].astype("string").str.strip()
    return frame


def parse_b18(path: Path) -> pd.DataFrame:
    """Read either the original fixed-width Berger 2018 table or its TSV export."""
    if path.suffix == ".gz":
        rows = []
        with gzip.open(path, "rt", errors="replace") as handle:
            for line in handle:
                try:
                    rows.append(
                        {
                            "kepid": line[0:8].strip(),
                            "radius_b18": float(line[67:74]),
                            "radius_b18_err_upper": float(line[75:81]),
                            "radius_b18_err_lower": float(line[82:88]),
                            "evol_b18": int(line[95:96]),
                            "bin_b18": int(line[97:98]),
                        }
                    )
                except (TypeError, ValueError):
                    continue
        return pd.DataFrame(rows).drop_duplicates("kepid")
    lines = path.read_text(errors="replace").splitlines()
    header = next(i for i, line in enumerate(lines) if line.startswith("KIC\t"))
    frame = pd.read_csv(Path(path), sep="\t", comment="#", skiprows=header, low_memory=False)
    frame = frame.rename(
        columns={
            "KIC": "kepid",
            "R*": "radius_b18",
            "E_R*": "radius_b18_err_upper",
            "e_R*": "radius_b18_err_lower",
            "Evol": "evol_b18",
            "Bin": "bin_b18",
        }
    )[[
        "kepid", "radius_b18", "radius_b18_err_upper",
        "radius_b18_err_lower", "evol_b18", "bin_b18"
    ]].drop_duplicates("kepid")
    frame["kepid"] = frame["kepid"].astype("string").str.strip()
    frame = frame[pd.to_numeric(frame["kepid"], errors="coerce").notna()].copy()
    return frame


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", default=".")
    parser.add_argument("--output", default="outputs/berger_density_provenance_audit")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    output = root / args.output
    output.parent.mkdir(parents=True, exist_ok=True)

    b20_path = root.parent / "table2.dat.gz"
    b18_path = root / "inputs" / "sagear_2026_published" / "berger2018_table1.dat.gz"
    if not b18_path.exists():
        b18_path = root.parent / "data" / "berger2018_table1_min.tsv"
    paper_path = root / "inputs" / "sagear_2026_published" / "journal_assets_20260723" / "Sagear_2026_AJ_172_42_published.txt"
    config = json.loads((root / "config.json").read_text())
    b20 = parse_b20(b20_path)
    b18 = parse_b18(b18_path)

    inventory_path = root / "outputs" / "sagear2026_planet_inventory_pre_visual_qc.csv"
    if inventory_path.exists():
        inventory = pd.read_csv(inventory_path, low_memory=False)
        key = "kepid" if "kepid" in inventory.columns else "KIC"
        inventory = inventory.rename(columns={key: "kepid"})
        hosts = inventory[["kepid"]].drop_duplicates()
        hosts["kepid"] = hosts["kepid"].astype("string")
    else:
        hosts = pd.DataFrame(columns=["kepid"])

    joined = hosts.merge(b20, on="kepid", how="left").merge(b18, on="kepid", how="left")
    mass = pd.to_numeric(joined["mass"], errors="coerce")
    radius = pd.to_numeric(joined["radius_b18"], errors="coerce")
    rho_log = pd.to_numeric(joined["rho_log"], errors="coerce")
    joined["rho_b20_from_mass_b18_radius"] = mass / radius**3
    joined["delta_logrho_b20_vs_hybrid"] = rho_log - np.log10(
        pd.to_numeric(joined["rho_b20_from_mass_b18_radius"], errors="coerce")
    )
    finite_delta = joined["delta_logrho_b20_vs_hybrid"].replace([np.inf, -np.inf], np.nan).dropna()

    paper = paper_path.read_text(errors="replace")
    result = {
        "berger2020_table_rows": int(len(b20)),
        "berger2020_unique_kic": int(b20["kepid"].nunique()),
        "berger2020_has_density_columns": bool(b20[["rho_log", "rho_log_upper", "rho_log_lower"]].notna().all().all()),
        "berger2018_rows": int(len(b18)),
        "berger2018_has_density_column": False,
        "inventory_unique_hosts": int(len(hosts)),
        "inventory_hosts_with_b20": int(joined["rho_log"].notna().sum()),
        "inventory_hosts_with_b18_radius": int(joined["radius_b18"].notna().sum()),
        "hybrid_density_comparison_n": int(len(finite_delta)),
        "hybrid_delta_logrho_median": float(finite_delta.median()) if len(finite_delta) else None,
        "hybrid_delta_logrho_p16": float(finite_delta.quantile(0.16)) if len(finite_delta) else None,
        "hybrid_delta_logrho_p84": float(finite_delta.quantile(0.84)) if len(finite_delta) else None,
        "paper_mentions_berger2018_density": "stellar sample and stellar densities presented by T. A. Berger\net al. (2018)" in paper,
        "paper_mentions_berger2020_metallicity_age": "T. A. Berger et al. (2020)" in paper,
        "config_density_source": config["paths"].get("berger_table2"),
    }
    (output.with_suffix(".json")).write_text(json.dumps(result, indent=2) + "\n")
    joined.to_csv(output.with_suffix(".csv"), index=False)
    lines = [
        "# Berger Density Provenance Audit",
        "",
        "## Finding",
        "Sagear's published methods prose names Berger et al. (2018) for the stellar-density prior, while the public source comments and conclusion point to Berger et al. (2020). Berger 2018 Table 1 supplies radius, evolution, and binary flags, but no stellar-density column, so the public materials do not uniquely identify an exact B18 density input.",
        "",
        "## Evidence",
        f"- Berger 2020 `table2.dat.gz`: {len(b20):,} rows, with mass, logg, radius, and asymmetric log-density columns.",
        f"- Berger 2018 table: {len(b18):,} rows, with radius/evolution/binary fields and no density field.",
        f"- Exact published-inventory hosts audited: {len(hosts):,}; B20 density coverage: {int(joined['rho_log'].notna().sum()):,}; B18 radius coverage: {int(joined['radius_b18'].notna().sum()):,}.",
        "- The operational extractor uses the B20 published density and its asymmetric uncertainties because that is the only public source with the required density posterior inputs and it is supported by the source comments/conclusion; this remains an explicit operational choice, not proof of Sagear's hidden run configuration.",
        "",
        "## Diagnostic comparison",
        f"For hosts with both sources, a hybrid density using B20 mass and B18 radius differs from B20's published density by median delta_log10(rho) = {result['hybrid_delta_logrho_median']:.4f} dex (16th-84th percentile {result['hybrid_delta_logrho_p16']:.4f} to {result['hybrid_delta_logrho_p84']:.4f}). This is a sensitivity diagnostic, not evidence that Sagear used the hybrid.",
        "",
        "## Decision",
        "Keep B20 as the operational branch, but do not describe it as definitively Sagear-equivalent. Treat B18-radius/B20-mass and the KG-RADII-derived B18 diagnostic as explicit sensitivity branches, and request author clarification about the apparent 2018/2020 provenance mismatch before calling the replication exact.",
    ]
    output.with_suffix(".md").write_text("\n".join(lines) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
