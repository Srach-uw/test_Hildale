from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.mixture import GaussianMixture
from scipy.special import logsumexp
from scipy.stats import multivariate_normal

from common import finite_columns, load_config, output_dir, read_angus, write_json
from toomre_diagnostics import ensure_velocity_aliases, read_old_angus_astropy_catalog


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Create Sagear-style Toomre diagrams with a P_thick background field."
    )
    parser.add_argument("--config", default=None)
    parser.add_argument("--sample", default=None)
    parser.add_argument("--classifier-model", default=None, help="Chemical classifier JSON written by diagnose_sample.py.")
    parser.add_argument("--max-background", type=int, default=80000)
    args = parser.parse_args()

    cfg = load_config(args.config)
    out_dir = output_dir()
    sample_path = Path(args.sample) if args.sample else out_dir / "canonical_sample_diagnostic.csv"
    sample = ensure_velocity_aliases(pd.read_csv(sample_path))
    angus = read_angus(cfg)
    old_angus = read_old_angus_astropy_catalog(cfg)

    model_path = Path(args.classifier_model) if args.classifier_model else out_dir / "chemical_classifier_strict.json"
    chemical = None
    chemical_path = None
    if model_path.exists() and old_angus is not None:
        with model_path.open("r", encoding="utf-8") as handle:
            model = json.load(handle)
        chemical = build_chemical_panel(sample, old_angus, model, args.max_background)
        chemical_path = out_dir / "toomre_sagear_style_chemical.png"
        plot_sagear_style(chemical, chemical_path, "Chemically calibrated cylindrical classifier")

    direct = build_panel_data(
        sample=sample,
        background=angus,
        x_col="V_phi",
        y_col="V_perp",
        bg_x_expr=lambda df: df["vy"],
        bg_y_expr=lambda df: np.sqrt(df["vx"] ** 2 + df["vz"] ** 2),
        display_sign=-1.0,
        max_background=args.max_background,
    )
    direct_path = out_dir / "toomre_sagear_style_direct.png"
    plot_sagear_style(direct, direct_path, "Current direct Angus velocity classifier")

    astropy_path = None
    astropy = None
    if old_angus is not None:
        astropy = build_panel_data(
            sample=sample,
            background=old_angus,
            x_col="V_phi_astropy",
            y_col="V_perp_astropy",
            bg_x_expr=lambda df: df["V_phi_astropy"],
            bg_y_expr=lambda df: df["V_perp_astropy"],
            display_sign=1.0,
            max_background=args.max_background,
        )
        astropy_path = out_dir / "toomre_sagear_style_old_astropy.png"
        plot_sagear_style(astropy, astropy_path, "Old-zip Astropy cylindrical convention")

    reference_path = (
        Path(cfg["_root"])
        / "external"
        / "old_hildale_project_zip"
        / "Shreshth_Hildale_Project-main"
        / "reference"
        / "Sagear_Fig2_toomre.png"
    )
    comparison_path = out_dir / "toomre_sagear_reference_comparison.png"
    if reference_path.exists():
        plot_reference_comparison(reference_path, chemical_path or direct_path, astropy_path, comparison_path)
    else:
        comparison_path = None

    summary = summarize_panels(sample, direct, astropy)
    if chemical is not None:
        summary = pd.concat([summary, summarize_panels(sample, chemical, None).assign(panel="chemical_classifier")])
    summary_path = out_dir / "toomre_sagear_style_summary.csv"
    summary.to_csv(summary_path, index=False)

    metadata = {
        "sample": str(sample_path),
        "direct_plot": str(direct_path),
        "chemical_classifier_model": str(model_path) if model_path.exists() else None,
        "chemical_plot": str(chemical_path) if chemical_path else None,
        "old_astropy_plot": str(astropy_path) if astropy_path else None,
        "reference_comparison": str(comparison_path) if comparison_path else None,
        "summary": str(summary_path),
        "note": (
            "The direct panel plots x=-Angus vy to match Sagear's negative V_phi display. "
            "The background color is the GMM-inferred P_thick field trained on the current planet-host sample."
        ),
    }
    metadata_path = out_dir / "toomre_sagear_style_metadata.json"
    write_json(metadata_path, metadata)

    print(summary.to_string(index=False))
    print(f"\nWrote: {direct_path}")
    if chemical_path:
        print(f"Wrote: {chemical_path}")
    if astropy_path:
        print(f"Wrote: {astropy_path}")
    if comparison_path:
        print(f"Wrote: {comparison_path}")
    print(f"Wrote: {summary_path}")


def build_panel_data(
    sample: pd.DataFrame,
    background: pd.DataFrame,
    x_col: str,
    y_col: str,
    bg_x_expr,
    bg_y_expr,
    display_sign: float,
    max_background: int,
) -> dict[str, object]:
    train_mask = finite_columns(sample, [x_col, y_col])
    train = sample.loc[train_mask, [x_col, y_col]].to_numpy(float)
    gmm = GaussianMixture(n_components=2, covariance_type="full", random_state=42, n_init=30)
    gmm.fit(train)
    score = -gmm.means_[:, 0] + 0.25 * gmm.means_[:, 1]
    thick_component = int(np.argmax(score))

    bg = background.copy()
    bg["panel_x_raw"] = bg_x_expr(bg)
    bg["panel_y"] = bg_y_expr(bg)
    bg = bg[finite_columns(bg, ["panel_x_raw", "panel_y"])].copy()
    bg = bg[(bg["panel_y"] >= 0) & (bg["panel_y"] <= 260) & (np.abs(bg["panel_x_raw"]) < 500)].copy()
    if len(bg) > max_background:
        bg = bg.sample(max_background, random_state=42).copy()

    bg_features = bg[["panel_x_raw", "panel_y"]].to_numpy(float)
    bg["P_thick_panel"] = gmm.predict_proba(bg_features)[:, thick_component]
    bg["panel_x"] = display_sign * bg["panel_x_raw"]

    hosts = sample.loc[train_mask].drop_duplicates("kepid").copy()
    hosts["panel_x"] = display_sign * hosts[x_col]
    hosts["panel_y"] = hosts[y_col]
    host_features = hosts[[x_col, y_col]].to_numpy(float)
    hosts["P_thick_panel"] = gmm.predict_proba(host_features)[:, thick_component]

    return {
        "background": bg,
        "hosts": hosts,
        "gmm_means": gmm.means_,
        "gmm_covariances": gmm.covariances_,
        "thick_component": thick_component,
        "display_sign": display_sign,
        "x_col": x_col,
        "y_col": y_col,
    }


def build_chemical_panel(
    sample: pd.DataFrame,
    background: pd.DataFrame,
    model: dict[str, object],
    max_background: int,
) -> dict[str, object]:
    if model.get("model") != "chemically_conditioned_gaussian_components":
        raise ValueError("Sagear-style chemical panel requires conditioned Gaussian-component metadata")
    x_col = str(model["x_column"])
    y_col = str(model["y_column"])
    required = [x_col, y_col]
    bg = background[finite_columns(background, required)].copy()
    bg = bg[(bg[y_col] >= 0) & (bg[y_col] <= 260) & (bg[x_col].abs() < 500)].copy()
    if len(bg) > max_background:
        bg = bg.sample(max_background, random_state=42).copy()
    features = bg[[x_col, y_col]].to_numpy(float)
    prior_high = float(model["prior_high"])
    low_log = multivariate_normal.logpdf(features, mean=model["low_mean"], cov=model["low_covariance"])
    high_log = multivariate_normal.logpdf(features, mean=model["high_mean"], cov=model["high_covariance"])
    joint = np.column_stack([low_log + np.log1p(-prior_high), high_log + np.log(prior_high)])
    bg["P_thick_panel"] = np.exp(joint[:, 1] - logsumexp(joint, axis=1))
    bg["panel_x"] = bg[x_col]
    bg["panel_y"] = bg[y_col]

    hosts = sample[finite_columns(sample, required)].drop_duplicates("kepid").copy()
    hosts["panel_x"] = hosts[x_col]
    hosts["panel_y"] = hosts[y_col]
    hosts["P_thick_panel"] = pd.to_numeric(hosts["P_thick"], errors="coerce")
    return {
        "background": bg,
        "hosts": hosts,
        "gmm_means": np.asarray([model["low_mean"], model["high_mean"]], dtype=float),
        "gmm_covariances": np.asarray([model["low_covariance"], model["high_covariance"]], dtype=float),
        "thick_component": 1,
        "display_sign": 1.0,
        "x_col": x_col,
        "y_col": y_col,
    }


def plot_sagear_style(panel: dict[str, object], path: Path, title: str) -> None:
    bg = panel["background"]
    hosts = panel["hosts"]
    fig, ax = plt.subplots(figsize=(9.5, 7.2), constrained_layout=True)
    scat = ax.scatter(
        bg["panel_x"],
        bg["panel_y"],
        c=bg["P_thick_panel"],
        s=4,
        cmap="plasma",
        vmin=0,
        vmax=1,
        alpha=0.34,
        linewidths=0,
        rasterized=True,
        label="KIC/Angus stars",
    )
    draw_toomre_arcs(ax)
    for disk, color, label in [
        ("thin", "#008b8b", "Kinematic Thin Disk Planet Hosts"),
        ("thick", "#9b1b1b", "Kinematic Thick Disk Planet Hosts"),
    ]:
        sub = hosts[hosts["disk"] == disk]
        ax.scatter(
            sub["panel_x"],
            sub["panel_y"],
            marker="+",
            s=46,
            c=color,
            linewidths=1.45,
            alpha=0.90,
            label=label,
        )
    cbar = fig.colorbar(scat, ax=ax, pad=0.02)
    cbar.set_label(r"$P_{\rm thick}$")
    ax.set_title(title)
    ax.set_xlabel(r"$V_\phi$ [km s$^{-1}$]")
    ax.set_ylabel(r"$\sqrt{V_R^2 + V_Z^2}$ [km s$^{-1}$]")
    ax.set_xlim(-330, -100)
    ax.set_ylim(0, 200)
    ax.legend(loc="upper left", frameon=True)
    fig.savefig(path, dpi=240)
    plt.close(fig)


def draw_toomre_arcs(ax) -> None:
    center = -229.0
    x = np.linspace(-350, -80, 900)
    for radius in [50, 100, 150, 200]:
        y2 = radius**2 - (x - center) ** 2
        ok = y2 >= 0
        ax.plot(x[ok], np.sqrt(y2[ok]), color="0.55", ls="--", lw=1.2, alpha=0.70)


def plot_reference_comparison(
    reference_path: Path,
    direct_path: Path,
    astropy_path: Path | None,
    out_path: Path,
) -> None:
    images = [(reference_path, "Sagear reference")]
    images.append((direct_path, "Chemically calibrated replication"))
    if astropy_path:
        images.append((astropy_path, "Old-zip Astropy convention"))

    fig, axes = plt.subplots(1, len(images), figsize=(7.2 * len(images), 6.2), constrained_layout=True)
    if len(images) == 1:
        axes = [axes]
    for ax, (path, title) in zip(axes, images):
        ax.imshow(mpimg.imread(path))
        ax.set_title(title)
        ax.axis("off")
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def summarize_panels(sample: pd.DataFrame, direct: dict[str, object], astropy: dict[str, object] | None) -> pd.DataFrame:
    rows = []
    for name, panel in [("direct_angus", direct), ("old_astropy", astropy)]:
        if panel is None:
            continue
        hosts = panel["hosts"]
        for disk in ["thin", "thick"]:
            for system in ["single", "multi"]:
                sub = sample[(sample["disk"] == disk) & (sample["system"] == system)]
                hsub = hosts[(hosts["disk"] == disk) & (hosts["system"] == system)]
                rows.append(
                    {
                        "panel": name,
                        "disk": disk,
                        "system": system,
                        "planets": int(len(sub)),
                        "hosts": int(sub["kepid"].nunique()),
                        "median_display_vphi": float(np.nanmedian(hsub["panel_x"])) if len(hsub) else np.nan,
                        "median_vperp": float(np.nanmedian(hsub["panel_y"])) if len(hsub) else np.nan,
                        "p90_vperp": float(np.nanpercentile(hsub["panel_y"], 90)) if len(hsub) else np.nan,
                        "median_panel_p_thick": float(np.nanmedian(hsub["P_thick_panel"])) if len(hsub) else np.nan,
                    }
                )
    return pd.DataFrame(rows)


if __name__ == "__main__":
    main()
