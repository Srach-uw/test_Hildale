from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import output_dir


def main() -> None:
    parser = argparse.ArgumentParser(description="Compare geometric-impact and paired-ALDERAAN-impact eccentricity extractions.")
    parser.add_argument("--geometric", required=True, help="Summary CSV extracted with impact_mode=geometric.")
    parser.add_argument("--alderaan", required=True, help="Summary CSV extracted with impact_mode=alderaan.")
    parser.add_argument("--out-prefix", default="impact_mode_comparison")
    args = parser.parse_args()

    geo = pd.read_csv(args.geometric)
    ald = pd.read_csv(args.alderaan)
    detail = geo.merge(
        ald,
        on="kepoi_name",
        suffixes=("_geometric", "_alderaan"),
    )
    detail["delta_e50_alderaan_minus_geometric"] = detail["e50_alderaan"] - detail["e50_geometric"]
    detail["delta_zeta_alderaan_minus_geometric"] = detail["zeta_median_alderaan"] - detail["zeta_median_geometric"]
    detail["e50_abs_delta"] = detail["delta_e50_alderaan_minus_geometric"].abs()

    group_cols = ["disk_geometric", "system_geometric"]
    summary = (
        detail.groupby(group_cols)
        .agg(
            n=("kepoi_name", "count"),
            geometric_e50_median=("e50_geometric", "median"),
            alderaan_e50_median=("e50_alderaan", "median"),
            median_delta_e50=("delta_e50_alderaan_minus_geometric", "median"),
            mean_delta_e50=("delta_e50_alderaan_minus_geometric", "mean"),
            geometric_e50_gt_0p8=("e50_geometric", lambda x: int((x > 0.8).sum())),
            alderaan_e50_gt_0p8=("e50_alderaan", lambda x: int((x > 0.8).sum())),
        )
        .reset_index()
        .rename(columns={"disk_geometric": "disk", "system_geometric": "system"})
    )

    out_dir = output_dir()
    detail_path = out_dir / f"{args.out_prefix}_detail.csv"
    summary_path = out_dir / f"{args.out_prefix}_summary.csv"
    plot_path = out_dir / f"{args.out_prefix}.png"
    detail.to_csv(detail_path, index=False)
    summary.to_csv(summary_path, index=False)
    make_plot(detail, plot_path)
    print(summary.to_string(index=False))
    print(f"\nWrote: {detail_path}")
    print(f"Wrote: {summary_path}")
    print(f"Wrote: {plot_path}")


def make_plot(detail: pd.DataFrame, path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True)
    colors = {
        ("thin", "single"): "#008b8b",
        ("thick", "single"): "#9b0000",
        ("thin", "multi"): "#76b7b2",
        ("thick", "multi"): "#c76f6f",
    }
    for (disk, system), sub in detail.groupby(["disk_geometric", "system_geometric"]):
        color = colors.get((disk, system), "0.4")
        label = f"{disk} {system} (n={len(sub)})"
        axes[0].scatter(sub["e50_geometric"], sub["e50_alderaan"], s=12, alpha=0.5, color=color, label=label)
        axes[1].hist(sub["delta_e50_alderaan_minus_geometric"], bins=np.linspace(-0.8, 0.8, 65), histtype="step", color=color)
    axes[0].plot([0, 0.95], [0, 0.95], color="0.4", ls="--", lw=1)
    axes[0].set_xlabel("geometric-impact e50")
    axes[0].set_ylabel("paired-ALDERAAN-impact e50")
    axes[0].set_xlim(0, 0.95)
    axes[0].set_ylim(0, 0.95)
    axes[1].axvline(0, color="0.4", ls="--", lw=1)
    axes[1].set_xlabel("paired - geometric e50")
    axes[1].set_ylabel("planets")
    axes[0].legend(fontsize=8)
    fig.savefig(path, dpi=220)
    plt.close(fig)


if __name__ == "__main__":
    main()
