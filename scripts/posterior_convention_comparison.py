"""Compare existing eccentricity posterior conventions on fixed population labels."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from hierarchical_rayleigh import fit_from_mass_matrix, load_population_masses


DEFAULT_PRODUCTS = {
    "dynesty_weighted": "outputs/eccentricity_posterior_summary_uniform_paired_full_qcprimary.csv",
    "equal_nested_diagnostic": "outputs/eccentricity_posterior_summary_equal_nested_uniform_full_50k_qcprimary.csv",
    "postmodel_unweighted": "outputs/eccentricity_posterior_summary_unweighted_merged_50k_qc.csv",
    "density_draw_validation": "outputs/eccentricity_posterior_summary_densitydraw_validation_50k.csv",
}

POPULATIONS = (
    ("thick", "single", "thick_singles"),
    ("thin", "single", "thin_singles"),
    ("thick", "multi", "thick_multis"),
    ("thin", "multi", "thin_multis"),
)


def compare_product(name: str, path: Path, sigma_count: int) -> list[dict[str, object]]:
    if not path.exists():
        return [{"product": name, "population": "all", "status": "missing", "path": str(path)}]
    summary = pd.read_csv(path)
    rows: list[dict[str, object]] = []
    sigmas = np.linspace(1e-4, 1.0, sigma_count)
    for disk, system, population in POPULATIONS:
        sub = summary[(summary.get("disk") == disk) & (summary.get("system") == system)].copy()
        row: dict[str, object] = {
            "product": name,
            "path": str(path),
            "population": population,
            "n": len(sub),
        }
        if len(sub) < 5:
            row["status"] = "too_few_rows"
            rows.append(row)
            continue
        try:
            masses, e_grid = load_population_masses(
                sub,
                apply_transit_selection=True,
                selection_mode="manuscript_reciprocal",
            )
            fit = fit_from_mass_matrix(
                masses,
                e_grid,
                sigmas,
                apply_transit_selection=True,
                selection_mode="manuscript_reciprocal",
            )
            row.update({"status": "ok", **fit})
        except Exception as exc:
            row.update({"status": f"error:{type(exc).__name__}", "error": str(exc)})
        rows.append(row)
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".")
    parser.add_argument("--output", default="outputs/posterior_convention_comparison.csv")
    parser.add_argument("--sigma-count", type=int, default=1200)
    args = parser.parse_args()
    root = Path(args.root).resolve()
    if args.sigma_count < 100:
        parser.error("--sigma-count must be at least 100")
    rows: list[dict[str, object]] = []
    for name, relative in DEFAULT_PRODUCTS.items():
        rows.extend(compare_product(name, root / relative, args.sigma_count))
    result = pd.DataFrame(rows)
    output = root / args.output
    output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output, index=False)
    print(result.to_string(index=False))
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()
