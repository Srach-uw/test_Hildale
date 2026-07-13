from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.io import fits

from common import normalize_dynesty_weights, output_dir


def main() -> None:
    parser = argparse.ArgumentParser(description="Make visual dossiers for suspicious ALDERAAN systems.")
    parser.add_argument("--flagged", default=None)
    parser.add_argument("--shape", default=None)
    parser.add_argument("--summary", default=None)
    parser.add_argument(
        "--results-dir",
        required=True,
    )
    parser.add_argument("--max-planets", type=int, default=40)
    args = parser.parse_args()

    out_dir = output_dir() / "alderaan_target_dossiers"
    out_dir.mkdir(parents=True, exist_ok=True)

    flagged = pd.read_csv(Path(args.flagged) if args.flagged else output_dir() / "alderaan_shape_diagnostics_flagged.csv")
    shape = pd.read_csv(Path(args.shape) if args.shape else output_dir() / "alderaan_shape_diagnostics.csv")
    summary = pd.read_csv(Path(args.summary) if args.summary else output_dir() / "eccentricity_posterior_summary.csv")

    targets = flagged.sort_values(["n_flags", "e50"], ascending=False).head(args.max_planets)
    rows = []
    for _, row in targets.iterrows():
        target = row["koi_target"]
        planet_name = row["kepoi_name"]
        fits_path = Path(args.results_dir) / f"{target}-results.fits"
        posterior_row = summary[summary["kepoi_name"] == planet_name]
        if not fits_path.exists() or posterior_row.empty:
            continue
        try:
            plot_path = make_planet_dossier(fits_path, row, posterior_row.iloc[0], out_dir)
            rows.append(
                {
                    "kepoi_name": planet_name,
                    "koi_target": target,
                    "disk": row["disk"],
                    "system": row["system"],
                    "e50": row["e50"],
                    "zeta_median": row["zeta_median"],
                    "n_flags": row["n_flags"],
                    "plot": str(plot_path),
                }
            )
        except Exception as exc:
            rows.append(
                {
                    "kepoi_name": planet_name,
                    "koi_target": target,
                    "disk": row.get("disk"),
                    "system": row.get("system"),
                    "e50": row.get("e50"),
                    "zeta_median": row.get("zeta_median"),
                    "n_flags": row.get("n_flags"),
                    "plot": "",
                    "error": str(exc)[:200],
                }
            )

    index = pd.DataFrame(rows)
    index_path = out_dir / "index.csv"
    index.to_csv(index_path, index=False)
    write_markdown_index(index, out_dir)
    print(f"Wrote {len(index)} dossier rows to: {index_path}")
    print(f"Wrote markdown index: {out_dir / 'index.md'}")


def make_planet_dossier(fits_path: Path, flagged_row: pd.Series, summary_row: pd.Series, out_dir: Path) -> Path:
    planet_idx = int(summary_row["alderaan_planet_index"])
    suffix = f"_{planet_idx}"
    with fits.open(fits_path, memmap=False) as hdul:
        samples = hdul["SAMPLES"].data
        weights = normalize_dynesty_weights(np.asarray(samples["LN_WT"], dtype=float))
        dur_hr = np.asarray(samples[f"DUR14{suffix}"], dtype=float) * 24.0
        ror = np.asarray(samples[f"ROR{suffix}"], dtype=float)
        impact = np.asarray(samples[f"IMPACT{suffix}"], dtype=float)

    posterior = np.load(summary_row["posterior_file"])
    e_grid = posterior["e_grid"]
    e_pdf = posterior["e_pdf"]

    fig, axes = plt.subplots(2, 2, figsize=(10, 7.2), constrained_layout=True)
    fig.suptitle(
        f"{flagged_row['kepoi_name']} | {flagged_row['disk']} {flagged_row['system']} | "
        f"e50={flagged_row['e50']:.3f}, zeta={flagged_row['zeta_median']:.3f}, flags={int(flagged_row['n_flags'])}",
        fontsize=12,
    )

    weighted_hist(axes[0, 0], dur_hr, weights, bins=60, color="#4477aa")
    axes[0, 0].axvline(float(flagged_row["koi_duration_hr"]), color="k", ls="--", lw=1.2, label="KOI")
    axes[0, 0].axvline(float(flagged_row["alderaan_dur14_hr_med"]), color="#aa3377", lw=1.2, label="ALDERAAN med")
    axes[0, 0].set_xlabel("T14 [hr]")
    axes[0, 0].legend(fontsize=8)

    weighted_hist(axes[0, 1], ror, weights, bins=60, color="#228833")
    axes[0, 1].axvline(float(flagged_row["koi_ror"]), color="k", ls="--", lw=1.2, label="KOI")
    axes[0, 1].axvline(float(flagged_row["alderaan_ror_med"]), color="#aa3377", lw=1.2, label="ALDERAAN med")
    axes[0, 1].set_xlabel("Rp/Rs")
    axes[0, 1].legend(fontsize=8)

    weighted_hist(axes[1, 0], impact, weights, bins=60, color="#ccbb44")
    if np.isfinite(float(flagged_row["koi_impact"])):
        axes[1, 0].axvline(float(flagged_row["koi_impact"]), color="k", ls="--", lw=1.2, label="KOI")
    axes[1, 0].axvline(float(flagged_row["alderaan_impact_med"]), color="#aa3377", lw=1.2, label="ALDERAAN med")
    axes[1, 0].set_xlabel("impact parameter")
    axes[1, 0].legend(fontsize=8)

    axes[1, 1].plot(e_grid, e_pdf, color="#882255", lw=1.6)
    axes[1, 1].axvline(float(flagged_row["e50"]), color="k", ls="--", lw=1.2, label="e50")
    axes[1, 1].set_xlabel("eccentricity")
    axes[1, 1].set_ylabel("posterior density")
    axes[1, 1].set_xlim(0, 0.95)
    axes[1, 1].legend(fontsize=8)

    for ax in axes.ravel():
        ax.grid(alpha=0.18)

    safe = str(flagged_row["kepoi_name"]).replace(".", "_")
    path = out_dir / f"{safe}_dossier.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def weighted_hist(ax, values: np.ndarray, weights: np.ndarray, bins: int, color: str) -> None:
    finite = np.isfinite(values) & np.isfinite(weights)
    values = values[finite]
    weights = weights[finite]
    if len(values) == 0:
        return
    weights = weights / weights.sum()
    ax.hist(values, bins=bins, weights=weights, histtype="stepfilled", alpha=0.35, color=color)
    ax.hist(values, bins=bins, weights=weights, histtype="step", lw=1.4, color=color)


def write_markdown_index(index: pd.DataFrame, out_dir: Path) -> None:
    lines = [
        "# ALDERAAN Target Dossiers",
        "",
        "Each PNG summarizes the existing FITS posterior for one suspicious planet: T14, Rp/Rs, impact parameter, and the derived eccentricity posterior.",
        "",
    ]
    for _, row in index.iterrows():
        rel = Path(row["plot"]).name if isinstance(row.get("plot"), str) and row.get("plot") else ""
        lines.append(
            f"- {row['kepoi_name']} ({row['disk']} {row['system']}): "
            f"e50={float(row['e50']):.3f}, zeta={float(row['zeta_median']):.3f}, flags={int(row['n_flags'])}"
            + (f" -> `{rel}`" if rel else f" -> ERROR: {row.get('error', '')}")
        )
    (out_dir / "index.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
