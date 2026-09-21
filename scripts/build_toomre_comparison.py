"""Place Sagear's Toomre figure beside the Table 1 replot."""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.image as mpimg
import matplotlib.pyplot as plt


def build_comparison(paper: Path, replot: Path, output: Path) -> Path:
    """Write a labeled side-by-side image from two existing provenance figures."""
    if not paper.is_file() or not replot.is_file():
        raise FileNotFoundError("paper and replot figure paths must exist")
    output.parent.mkdir(parents=True, exist_ok=True)
    paper_image = mpimg.imread(paper)
    replot_image = mpimg.imread(replot)
    figure, axes = plt.subplots(1, 2, figsize=(18, 7.5), constrained_layout=True)
    for axis, image, title in zip(
        axes,
        [paper_image, replot_image],
        [
            "Sagear et al. Figure 2",
            "This project: replot using published Table 1 velocities and labels",
        ],
    ):
        axis.imshow(image)
        axis.set_title(title, fontsize=13)
        axis.axis("off")
    figure.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(figure)
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--paper", type=Path, default=Path("reference/Sagear_Fig2_toomre.png"))
    parser.add_argument(
        "--replot",
        type=Path,
        default=Path("metadata/recovery_preflight_20260724/toomre_sagear_published_truth.png"),
    )
    parser.add_argument("--output", type=Path, default=Path("figures/toomre_published_comparison.png"))
    args = parser.parse_args()
    print(build_comparison(args.paper, args.replot, args.output))


if __name__ == "__main__":
    main()
