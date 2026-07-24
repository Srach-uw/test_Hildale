from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


ROWS = [
    {
        "contract": "published planet total",
        "sagear_statement": "2465 planets / 1888 stars",
        "local_implementation": "Published macro audit; reconstructed inventory is kept separate",
        "evidence": "published source main.tex macros and selection paragraph",
        "status": "documented discrepancy",
        "risk": "high",
    },
    {
        "contract": "population subgroup targets",
        "sagear_statement": "thin singles 1121; thick singles 275; thin multis 862; thick multis 207",
        "local_implementation": "Targets recorded independently of the 547-fit coverage",
        "evidence": "published source main.tex macros",
        "status": "target locked",
        "risk": "high",
    },
    {
        "contract": "single/multi definition",
        "sagear_statement": "single or multi host categories are defined for the sample, before posterior availability",
        "local_implementation": "pre-cut host multiplicity map; strict validation and alternate-definition audit",
        "evidence": "common.py, reconstruct_published_planet_inventory.py, multiplicity tests",
        "status": "implemented and tested",
        "risk": "high",
    },
    {
        "contract": "disk labels",
        "sagear_statement": "APOGEE chemical thin/thick calibration plus GMM Galactocentric velocities; P_thick threshold",
        "local_implementation": "published inventory labels are canonical for replication; old Toomre classifier is diagnostic only",
        "evidence": "published methods and published_toomre_truth_audit.csv",
        "status": "mostly locked; source data provenance remains explicit",
        "risk": "high",
    },
    {
        "contract": "stellar density prior",
        "sagear_statement": "Berger et al. 2018 in final text; source draft equation cites Berger et al. 2020",
        "local_implementation": "Berger 2020 density branch retained; 2018 citation inconsistency documented",
        "evidence": "Berger provenance audit and published source",
        "status": "unresolved source ambiguity",
        "risk": "high",
    },
    {
        "contract": "transit-shape posterior",
        "sagear_statement": "ALDERAAN nested sampling; shape parameters sampled simultaneously by system; long/short cadence integration",
        "local_implementation": "ALDERAAN FITS; 82-factor validation completed",
        "evidence": "82 FITS factorial and transit_shape_contract_audit.csv",
        "status": "implementation tested, exact unpublished settings not fully recoverable",
        "risk": "high",
    },
    {
        "contract": "printed ALDERAAN priors",
        "sagear_statement": "C0/C1 N(0,1); Rp/Rstar U(1e-5,0.99); b U(0,1+Rp/Rstar)",
        "local_implementation": "pinned public code uses C0/C1 N(0,0.1) and log-uniform Rp/Rstar; patched arm exists",
        "evidence": "journal Table 2 asset, dynesty_helpers.py, paper-prior factorial arm",
        "status": "validated sensitivity only; full effect unresolved",
        "risk": "high",
    },
    {
        "contract": "post-model proposals",
        "sagear_statement": "e uniform [0, 0.95], omega uniform [-pi/2, 3pi/2], exact circular-density relation",
        "local_implementation": "direct extractor uses paired ALDERAAN T14, Rp/R*, b and exact formula",
        "evidence": "published source commented formalism and extract_eccentricity_posteriors_direct.py",
        "status": "implemented and tested",
        "risk": "medium",
    },
    {
        "contract": "posterior export weighting",
        "sagear_statement": "source text says weighted samples are converted to unweighted draws",
        "local_implementation": "Dynesty-weighted and equal-row branches both retained; canonical interpretation is unresolved",
        "evidence": "MacDougall/Gilbert source and ALDERAAN Results.py",
        "status": "critical unresolved fork",
        "risk": "critical",
    },
    {
        "contract": "hierarchical selection correction",
        "sagear_statement": "source equation uses reciprocal transit-probability factor in the likelihood",
        "local_implementation": "manuscript_reciprocal implements that equation; selection modes are explicit",
        "evidence": "published source comments and hierarchical_rayleigh.py",
        "status": "implemented and tested",
        "risk": "medium",
    },
    {
        "contract": "fit exclusions",
        "sagear_statement": "visual inspection; less than 2 percent removed for non-convergence",
        "local_implementation": "mechanical ALDERAAN status/QC plus explicit exclusion manifest; manual decisions unavailable",
        "evidence": "published methods and local QC audit",
        "status": "unresolved unpublished selection",
        "risk": "high",
    },
    {
        "contract": "posterior coverage",
        "sagear_statement": "photoeccentric posteriors for virtually all sample planets",
        "local_implementation": "547 systems / 2395 planet summaries in current archive; missing set is non-random",
        "evidence": "missing_coverage_selection_bias.csv and cloud archive manifest",
        "status": "not complete",
        "risk": "critical",
    },
    {
        "contract": "journal-hosted data assets",
        "sagear_statement": "published article exposes ASCII Table 2 priors and Table 3 results",
        "local_implementation": "downloaded and archived both assets; no per-planet posterior asset is present",
        "evidence": "inputs/sagear_2026_published/journal_assets_20260724/table2_priors_20260724.txt and table3_results_20260724.txt",
        "status": "checked; no posterior release",
        "risk": "high",
    },
    {
        "contract": "hierarchical uncertainty",
        "sagear_statement": "numpyro, two chains, Rhat < 1.05, leave-10-percent validation",
        "local_implementation": "grid Rayleigh diagnostic plus leave-out outputs; final HBM equivalence still needs exact sample contract",
        "evidence": "published methods and local hierarchy diagnostics",
        "status": "partially reproduced",
        "risk": "medium",
    },
]


def build_matrix() -> pd.DataFrame:
    return pd.DataFrame(ROWS, columns=["contract", "sagear_statement", "local_implementation", "evidence", "status", "risk"])


def main() -> None:
    parser = argparse.ArgumentParser(description="Write the Sagear versus local replication contract matrix.")
    parser.add_argument("--out", default="outputs/source_contract_matrix.csv")
    parser.add_argument("--markdown-out", default="outputs/source_contract_matrix.md")
    args = parser.parse_args()
    matrix = build_matrix()
    out = Path(args.out)
    md = Path(args.markdown_out)
    out.parent.mkdir(parents=True, exist_ok=True)
    md.parent.mkdir(parents=True, exist_ok=True)
    matrix.to_csv(out, index=False)
    md.write_text(
        "# Sagear Replication Source Contract Matrix\n\n"
        "This is an audit artifact, not a claim that every local choice is equivalent to Sagear. "
        "A row is only marked implemented when the behavior is explicit and tested.\n\n"
        + matrix.to_markdown(index=False)
        + "\n",
        encoding="utf-8",
    )
    print(matrix.to_string(index=False))
    print(f"Wrote {out}")
    print(f"Wrote {md}")


if __name__ == "__main__":
    main()
