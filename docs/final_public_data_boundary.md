# Final public-data boundary

Updated: 2026-08-13

The reconstruction now uses the exact finite-duration equation described by
MacDougall et al. and printed in Sagear et al.'s later radius analysis. It keeps
the ALDERAAN duration, radius-ratio, impact, and period samples paired and uses
the Dynesty nested weights stored in each FITS result.

## Closed implementation questions

- Replacing paired ALDERAAN periods with catalog periods is negligible. Among
  2,123 source-faithful posterior rows, the median relative period difference
  is `2.18e-6`, the 99th percentile is `7.52e-5`, and the maximum is
  `4.27e-4`. The leading period contribution to density is below `0.086%` in
  the worst case.
- Raw ALDERAAN nested rows cannot be treated as equal-weight posterior draws.
  The result files store `LN_WT`, `LN_LIKE`, and `LN_Z`. The public ALDERAAN
  reader evaluates `exp(LN_WT - LN_Z[-1])` before computing summaries.
- Swapping thin and thick labels produces a numerical resemblance for several
  point estimates, but the population sizes become 1,028 thick singles and 256
  thin singles instead of the published 275 and 1,121. The interval widths also
  follow the swapped sample sizes. This is not a valid interpretation.
- The public Sagear ALDERAAN fork is code-identical to the pinned Gilbert commit
  used for the cloud fits.
- The released 2023 hierarchy matrices reproduce the older M-dwarf population
  scale under their archived likelihood. Applying the same log base and lower
  Rayleigh boundary to the current disk sample does not recover Table 3.

## Stellar-density ambiguity

The disk paper gives conflicting catalog citations. The Methods section names
Berger et al. (2018), while the commented importance-sampling description and
the Conclusions name Berger et al. (2020). The public Berger 2018 table has
stellar radius and radius errors but no stellar density, stellar mass, or
mass-radius covariance. A richer Kepler-Gaia FITS product adds surface gravity
and its uncertainty, but it still does not contain the density posterior needed
to reproduce the exact input.

Source-supported Berger 2018 reconstructions and independent Berger 2020,
California-Kepler Survey, and asteroseismic controls do not reproduce the
published thin-single result. They remain sensitivity tests because none can
recover an unpublished density ledger exactly.

## Remaining request

An exact numerical replication now requires two planet-specific products:

1. the final accepted-fit table with KOI identifiers and rejection reasons;
2. the stellar-density rows or posterior samples supplied to each eccentricity
   calculation, ideally with the resulting `(e, omega)` posterior samples.

The repository therefore reports a public-data replication discrepancy. It
does not treat that discrepancy as evidence against the paper's astrophysical
interpretation.
