from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.mixture import GaussianMixture

from common import finite_columns, load_config, output_dir, read_angus
from toomre_diagnostics import ensure_velocity_aliases, read_old_angus_astropy_catalog


TARGETS = {
    "thin_single": 1121,
    "thick_single": 275,
    "thin_multi": 862,
    "thick_multi": 207,
}


def main() -> None:
    parser = argparse.ArgumentParser(description="Compare Toomre classifier conventions against Sagear-style counts.")
    parser.add_argument("--sample", default=None)
    parser.add_argument("--config", default=None)
    parser.add_argument("--max-background", type=int, default=45000)
    parser.add_argument("--max-train", type=int, default=50000)
    args = parser.parse_args()

    cfg = load_config(args.config)
    out_dir = output_dir()
    sample_path = Path(args.sample) if args.sample else out_dir / "canonical_sample_diagnostic.csv"
    sample = ensure_velocity_aliases(pd.read_csv(sample_path))
    angus_direct = read_angus(cfg)
    old_angus = read_old_angus_astropy_catalog(cfg)
    if old_angus is None:
        raise RuntimeError("Could not parse raw Angus table for old-Astropy convention.")

    variants = [
        build_variant(
            name="direct_sample_trained",
            title="Direct Angus, planet-sample GMM",
            sample=sample,
            background=angus_direct,
            sample_x="V_phi",
            sample_y="V_perp",
            bg_x=lambda df: df["vy"],
            bg_y=lambda df: np.sqrt(df["vx"] ** 2 + df["vz"] ** 2),
            train_source="sample",
            display_sign=-1.0,
            max_background=args.max_background,
            max_train=args.max_train,
        ),
        build_variant(
            name="direct_kic_trained",
            title="Direct Angus, KIC-wide GMM",
            sample=sample,
            background=angus_direct,
            sample_x="V_phi",
            sample_y="V_perp",
            bg_x=lambda df: df["vy"],
            bg_y=lambda df: np.sqrt(df["vx"] ** 2 + df["vz"] ** 2),
            train_source="background",
            display_sign=-1.0,
            max_background=args.max_background,
            max_train=args.max_train,
        ),
        build_variant(
            name="old_astropy_sample_trained",
            title="Old Astropy, planet-sample GMM",
            sample=sample,
            background=old_angus,
            sample_x="V_phi_astropy",
            sample_y="V_perp_astropy",
            bg_x=lambda df: df["V_phi_astropy"],
            bg_y=lambda df: df["V_perp_astropy"],
            train_source="sample",
            display_sign=1.0,
            max_background=args.max_background,
            max_train=args.max_train,
        ),
        build_variant(
            name="old_astropy_kic_trained",
            title="Old Astropy, KIC-wide GMM",
            sample=sample,
            background=old_angus,
            sample_x="V_phi_astropy",
            sample_y="V_perp_astropy",
            bg_x=lambda df: df["V_phi_astropy"],
            bg_y=lambda df: df["V_perp_astropy"],
            train_source="background",
            display_sign=1.0,
            max_background=args.max_background,
            max_train=args.max_train,
        ),
    ]

    counts = pd.DataFrame([v["counts"] for v in variants])
    counts_path = out_dir / "toomre_classifier_grid_counts.csv"
    counts.to_csv(counts_path, index=False)

    grid_path = out_dir / "toomre_classifier_grid.png"
    plot_grid(variants, grid_path)

    reference_path = (
        Path(cfg["_root"])
        / "external"
        / "old_hildale_project_zip"
        / "Shreshth_Hildale_Project-main"
        / "reference"
        / "Sagear_Fig2_toomre.png"
    )
    comparison_path = out_dir / "toomre_classifier_grid_with_reference.png"
    if reference_path.exists():
        plot_with_reference(reference_path, grid_path, comparison_path)

    print(counts.to_string(index=False))
    print(f"\nWrote: {counts_path}")
    print(f"Wrote: {grid_path}")
    if reference_path.exists():
        print(f"Wrote: {comparison_path}")


def build_variant(
    name: str,
    title: str,
    sample: pd.DataFrame,
    background: pd.DataFrame,
    sample_x: str,
    sample_y: str,
    bg_x,
    bg_y,
    train_source: str,
    display_sign: float,
    max_background: int,
    max_train: int,
) -> dict[str, object]:
    bg = background.copy()
    bg["raw_x"] = bg_x(bg)
    bg["raw_y"] = bg_y(bg)
    bg = bg[finite_columns(bg, ["raw_x", "raw_y"])].copy()
    display_x_for_filter = display_sign * bg["raw_x"]
    bg = bg[(np.abs(display_x_for_filter + 220) < 220) & (bg["raw_y"] < 300)].copy()

    sample_ok = finite_columns(sample, [sample_x, sample_y])
    if train_source == "sample":
        train = sample.loc[sample_ok, [sample_x, sample_y]].to_numpy(float)
    elif train_source == "background":
        train_df = bg
        if len(train_df) > max_train:
            train_df = train_df.sample(max_train, random_state=42)
        train = train_df[["raw_x", "raw_y"]].to_numpy(float)
    else:
        raise ValueError(train_source)

    gmm = GaussianMixture(n_components=2, covariance_type="full", random_state=42, n_init=8)
    gmm.fit(train)
    thick_component = int(np.argmax(gmm.means_[:, 1]))

    sample_probs = np.full(len(sample), np.nan)
    sample_probs[sample_ok.to_numpy()] = gmm.predict_proba(sample.loc[sample_ok, [sample_x, sample_y]].to_numpy(float))[
        :, thick_component
    ]
    assigned = np.where(sample_probs > 0.5, "thick", "thin").astype(object)
    assigned[~np.isfinite(sample_probs)] = "unclassified"

    if len(bg) > max_background:
        bg_plot = bg.sample(max_background, random_state=42).copy()
    else:
        bg_plot = bg.copy()
    bg_plot["P_thick"] = gmm.predict_proba(bg_plot[["raw_x", "raw_y"]].to_numpy(float))[:, thick_component]
    bg_plot["plot_x"] = display_sign * bg_plot["raw_x"]
    bg_plot["plot_y"] = bg_plot["raw_y"]

    hosts = sample.loc[sample_ok].drop_duplicates("kepid").copy()
    host_probs = gmm.predict_proba(hosts[[sample_x, sample_y]].to_numpy(float))[:, thick_component]
    hosts["plot_x"] = display_sign * hosts[sample_x]
    hosts["plot_y"] = hosts[sample_y]
    hosts["assigned_disk"] = np.where(host_probs > 0.5, "thick", "thin")
    hosts["P_thick"] = host_probs

    counts = count_assigned(sample, assigned)
    counts.update(
        {
            "variant": name,
            "train_source": train_source,
            "means": repr(np.round(gmm.means_, 3).tolist()),
            "weights": repr(np.round(gmm.weights_, 3).tolist()),
            "thick_component": thick_component,
        }
    )
    return {"name": name, "title": title, "background": bg_plot, "hosts": hosts, "counts": counts}


def count_assigned(sample: pd.DataFrame, assigned: np.ndarray) -> dict[str, object]:
    row: dict[str, object] = {}
    for disk in ["thin", "thick"]:
        for system in ["single", "multi"]:
            key = f"{disk}_{system}"
            row[key] = int(((assigned == disk) & (sample["system"] == system)).sum())
            row[f"delta_{key}"] = row[key] - TARGETS[key]
    row["l1_planet_count_error"] = int(
        sum(abs(row[f"delta_{key}"]) for key in ["thin_single", "thick_single", "thin_multi", "thick_multi"])
    )
    return row


def plot_grid(variants: list[dict[str, object]], path: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(15, 11), constrained_layout=True)
    for ax, variant in zip(axes.ravel(), variants):
        bg = variant["background"]
        hosts = variant["hosts"]
        sc = ax.scatter(
            bg["plot_x"],
            bg["plot_y"],
            c=bg["P_thick"],
            cmap="plasma",
            vmin=0,
            vmax=1,
            s=3,
            alpha=0.34,
            linewidths=0,
            rasterized=True,
        )
        draw_arcs(ax)
        for disk, color in [("thin", "#008b8b"), ("thick", "#9b1b1b")]:
            sub = hosts[hosts["assigned_disk"] == disk]
            ax.scatter(sub["plot_x"], sub["plot_y"], marker="+", s=38, c=color, linewidths=1.3, alpha=0.9)
        c = variant["counts"]
        ax.set_title(
            f"{variant['title']}\n"
            f"S thin/thick={c['thin_single']}/{c['thick_single']}; "
            f"M thin/thick={c['thin_multi']}/{c['thick_multi']}; L1={c['l1_planet_count_error']}"
        )
        ax.set_xlim(-330, -100)
        ax.set_ylim(0, 200)
        ax.set_xlabel(r"$V_\phi$ [km s$^{-1}$]")
        ax.set_ylabel(r"$\sqrt{V_R^2 + V_Z^2}$ [km s$^{-1}$]")
    fig.colorbar(sc, ax=axes.ravel().tolist(), pad=0.01, label=r"$P_{\rm thick}$")
    fig.savefig(path, dpi=220)
    plt.close(fig)


def draw_arcs(ax) -> None:
    center = -229.0
    x = np.linspace(-350, -80, 900)
    for radius in [50, 100, 150, 200]:
        y2 = radius**2 - (x - center) ** 2
        ok = y2 >= 0
        ax.plot(x[ok], np.sqrt(y2[ok]), color="0.55", ls="--", lw=1.1, alpha=0.70)


def plot_with_reference(reference_path: Path, grid_path: Path, out_path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(18, 7.4), constrained_layout=True)
    for ax, path, title in [
        (axes[0], reference_path, "Sagear reference"),
        (axes[1], grid_path, "Classifier convention grid"),
    ]:
        ax.imshow(mpimg.imread(path))
        ax.set_title(title)
        ax.axis("off")
    fig.savefig(out_path, dpi=170)
    plt.close(fig)


if __name__ == "__main__":
    main()
