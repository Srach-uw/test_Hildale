# K00367 diagnostic input preflight

Date: 2026-09-24

This record covers one completed preparation check for the standalone
photometry-to-density diagnostic. It is not a transit injection, ALDERAAN
rerun, eccentricity extraction, or population result.

## What was checked

The preparation command ran in the pinned Linux ALDERAAN environment with
Lightkurve 2.1.1. It used 14 local Kepler long-cadence PDCSAP files for K00367,
the archived `K00367-results.fits`, and the local full-system catalog presented
to ALDERAAN. The result FITS and catalog agreed on one fitted planet. The
PDCSAP files agreed on KIC 4815520, `OBSMODE=long cadence`, `TIMEDEL`
0.02043359821692 days, and the Kepler BKJD, TDB time convention.

The declared exposure was 29.4243814323648 minutes. It matched `TIMEDEL *
1440` exactly at the stored precision. The command used the pinned
`KeplerQualityFlags.create_quality_mask(..., bitmask="default")` implementation,
not a substitute quality rule.

## Result

| Quantity | Value |
| --- | ---: |
| Raw cadences | 54,399 |
| Cadences passing the default quality mask | 50,856 |
| Finite cadences | 50,642 |
| Retained cadences | 50,642 |
| Cadences excluded around known model transits | 506 |
| Injection-eligible cadences | 50,136 |
| System-wide exclusion half-width | 0.16 days |

The local prepared NPZ had SHA-256
`06bcd8b2d281a53e57edd465a92530cf9942595f81f2817d8ddddb217cc3738a`.
The local manifest had SHA-256
`508a972a92c7859d146027f677085cd8ab3c6d493f5df88b59811b6c584f20c8`.
Both digests matched the independently produced Linux files.

## Interpretation

This closes the input-provenance gate for K00367. It does not yet test whether
the recovery fit returns the injected circular density. That next stage must
render exposure-integrated synthetic transits into this prepared baseline,
exclude `known_transit_mask` during recovery, and keep the transit geometry
free enough to measure circular density rather than impose it.
