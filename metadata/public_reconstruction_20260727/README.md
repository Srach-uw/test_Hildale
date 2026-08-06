# Public Reconstruction Evidence

This directory contains the compact evidence for the July 27 public-data
reconstruction. It intentionally excludes posterior NPZ files, raw ALDERAAN
archives, and local filesystem paths.

`photoeccentric_density_audit.csv` reports a population-level circular-density
diagnostic for all 2,465 reconstructed planets. It joins the paired-period
posterior ledger to the one-row-per-planet raw-FITS audit. The accompanying JSON
records both input hashes and the valid-row count.
`photoeccentric_impact_width_audit.csv` separates each population at
the exploratory threshold `b84 - b16 = 0.4` using the same weighted raw FITS.
The result describes how density uncertainty changes with impact-posterior
width; it is not an unbiased population split.
`rayleigh_impact_constrained_sensitivity.csv` gives the forward-normalized and
arXiv-v1 reciprocal Rayleigh fits for those constrained subsets after applying
the primary posterior-QC exclusion. It records the number removed from each
cell explicitly.

The two Rayleigh tables are the complete and deterministic-QC sensitivities.
They use paired ALDERAAN transit-shape and period samples, dynesty weights,
fixed Berger et al. (2020) stellar density, and the manuscript selection
correction.

The inventory and contract audits document the 2,465-planet accounting,
published-label reconciliation, method provenance, and tests of possible
population-order mistakes. These files are diagnostic evidence, not a claim of
an independent astrophysical measurement.
