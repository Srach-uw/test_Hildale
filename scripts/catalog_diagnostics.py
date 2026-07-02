from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import output_dir


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compare the current Gaia-Kepler/Berger-2020-style sample against Berger+2018."
    )
    parser.add_argument("--sample", default=None, help="Canonical sample CSV to inspect.")
    parser.add_argument("--berger2018", default=None, help="VizieR TSV for J/ApJ/866/99/table1.")
    parser.add_argument("--out-prefix", default="catalog_diagnostics")
    args = parser.parse_args()

    out_dir = output_dir()
    root = out_dir.parent.parent
    sample_path = Path(args.sample) if args.sample else out_dir / "canonical_sample_diagnostic.csv"
    berger_path = Path(args.berger2018) if args.berger2018 else root / "data" / "berger2018_table1_min.tsv"

    sample = pd.read_csv(sample_path)
    berger = read_vizier_tsv(berger_path)
    berger = normalize_berger2018(berger)

    merged = sample.merge(berger, on="kepid", how="left", suffixes=("", "_berger2018"))
    merged["has_berger2018"] = merged["berger2018_radius"].notna()

    summary = sample_summary(merged, sample_path, berger_path)
    missing = missing_by_population(merged)
    radius = radius_comparison(merged)
    outliers = annotate_existing_outliers(out_dir, berger)

    summary_path = out_dir / f"{args.out_prefix}_summary.csv"
    missing_path = out_dir / f"{args.out_prefix}_missing_by_population.csv"
    radius_path = out_dir / f"{args.out_prefix}_radius_comparison.csv"
    outlier_path = out_dir / f"{args.out_prefix}_thin_single_outliers_berger2018.csv"
    plot_path = out_dir / f"{args.out_prefix}_radius_comparison.png"

    summary.to_csv(summary_path, index=False)
    missing.to_csv(missing_path, index=False)
    radius.to_csv(radius_path, index=False)
    outliers.to_csv(outlier_path, index=False)
    make_radius_plot(radius, plot_path)

    print("=== Catalog Diagnostics ===")
    print(summary.to_string(index=False))
    print("\n=== Missing Berger 2018 by population ===")
    print(missing.to_string(index=False))
    print(f"\nWrote: {summary_path}")
    print(f"Wrote: {missing_path}")
    print(f"Wrote: {radius_path}")
    print(f"Wrote: {outlier_path}")
    print(f"Wrote: {plot_path}")


def read_vizier_tsv(path: Path) -> pd.DataFrame:
    lines = path.read_text(errors="replace").splitlines()
    header_idx = next((i for i, line in enumerate(lines) if line.startswith("KIC\t")), None)
    if header_idx is None:
        header_idx = next((i for i, line in enumerate(lines) if line.startswith("Np\t")), None)
    if header_idx is None:
        raise ValueError(f"Could not find VizieR data header in {path}")

    # VizieR TSV has a header row, a units row, a dashed separator row, then data.
    data = "\n".join([lines[header_idx]] + lines[header_idx + 3 :])
    from io import StringIO

    return pd.read_csv(StringIO(data), sep="\t")


def normalize_berger2018(df: pd.DataFrame) -> pd.DataFrame:
    rename = {
        "KIC": "kepid",
        "R*": "berger2018_radius",
        "E_R*": "berger2018_radius_err_plus",
        "e_R*": "berger2018_radius_err_minus",
        "Teff": "berger2018_teff",
        "Evol": "berger2018_evol",
        "Bin": "berger2018_bin",
        "Gaia": "berger2018_gaia_dr2",
    }
    out = df.rename(columns=rename)
    keep = [c for c in rename.values() if c in out.columns]
    out = out[keep].copy()
    out["kepid"] = pd.to_numeric(out["kepid"], errors="coerce").astype("Int64")
    for col in out.columns:
        if col != "kepid":
            out[col] = pd.to_numeric(out[col], errors="coerce")
    out = out.dropna(subset=["kepid"]).drop_duplicates("kepid")
    out["kepid"] = out["kepid"].astype(int)
    return out


def sample_summary(merged: pd.DataFrame, sample_path: Path, berger_path: Path) -> pd.DataFrame:
    rows = [
        {
            "sample": str(sample_path),
            "berger2018_table": str(berger_path),
            "planets": len(merged),
            "hosts": merged["kepid"].nunique(),
            "planets_with_berger2018": int(merged["has_berger2018"].sum()),
            "hosts_with_berger2018": int(merged.loc[merged["has_berger2018"], "kepid"].nunique()),
            "planets_missing_berger2018": int((~merged["has_berger2018"]).sum()),
            "hosts_missing_berger2018": int(merged.loc[~merged["has_berger2018"], "kepid"].nunique()),
        }
    ]
    return pd.DataFrame(rows)


def missing_by_population(merged: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (disk, system), sub in merged.groupby(["disk", "system"], dropna=False):
        missing = sub[~sub["has_berger2018"]]
        rows.append(
            {
                "disk": disk,
                "system": system,
                "planets": len(sub),
                "hosts": sub["kepid"].nunique(),
                "missing_planets": len(missing),
                "missing_hosts": missing["kepid"].nunique(),
                "missing_planet_fraction": len(missing) / len(sub) if len(sub) else np.nan,
            }
        )
    return pd.DataFrame(rows).sort_values(["disk", "system"])


def radius_comparison(merged: pd.DataFrame) -> pd.DataFrame:
    if "berger_rad" not in merged.columns:
        return pd.DataFrame([{"note": "Current sample does not contain berger_rad."}])
    cols = [
        "kepid",
        "kepoi_name",
        "disk",
        "system",
        "berger_rad",
        "berger2018_radius",
        "berger2018_teff",
        "berger2018_evol",
        "berger2018_bin",
    ]
    out = merged.loc[merged["has_berger2018"], [c for c in cols if c in merged.columns]].copy()
    out["radius_ratio_current_over_2018"] = pd.to_numeric(out["berger_rad"], errors="coerce") / pd.to_numeric(
        out["berger2018_radius"], errors="coerce"
    )
    return out


def annotate_existing_outliers(out_dir: Path, berger: pd.DataFrame) -> pd.DataFrame:
    top_path = out_dir / "eccentricity_diagnostics_thin_single_top_outliers.csv"
    if not top_path.exists():
        return pd.DataFrame([{"note": "Run eccentricity_diagnostics.py first."}])
    top = pd.read_csv(top_path)
    merged = top.merge(berger, on="kepid", how="left")
    merged["has_berger2018"] = merged["berger2018_radius"].notna()
    if "berger_rad" in merged.columns:
        merged["radius_ratio_current_over_2018"] = pd.to_numeric(merged["berger_rad"], errors="coerce") / pd.to_numeric(
            merged["berger2018_radius"], errors="coerce"
        )
    return merged


def make_radius_plot(radius: pd.DataFrame, path: Path) -> None:
    if "radius_ratio_current_over_2018" not in radius.columns or radius.empty:
        return
    ratio = pd.to_numeric(radius["radius_ratio_current_over_2018"], errors="coerce")
    ratio = ratio[np.isfinite(ratio)]
    ratio = ratio[(ratio > 0) & (ratio < 5)]
    fig, ax = plt.subplots(figsize=(7, 4.5), constrained_layout=True)
    ax.hist(ratio, bins=np.linspace(0.5, 1.5, 50), color="#4f8f8f", alpha=0.8)
    ax.axvline(1.0, color="black", linestyle="--", linewidth=1)
    ax.set_xlabel("current Berger-style radius / Berger 2018 radius")
    ax.set_ylabel("planet rows")
    ax.set_title("Catalog Radius Comparison")
    fig.savefig(path, dpi=220)
    plt.close(fig)


if __name__ == "__main__":
    main()
