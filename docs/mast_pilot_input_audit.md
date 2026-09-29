# Bounded MAST pilot header audit

This utility is a narrow, read-only availability check for the 24 predeclared
targets in the population-noise audit plan. It follows the plan's input gate:
it determines whether the exact products listed in a previously verified MAST
directory listing expose basic FITS identity and cadence metadata. It is not a
new light-curve reduction and it does not authorize a transit fit.

## Inputs and contract

Run `scripts/audit_mast_pilot_fits_headers.py` with explicit paths for:

- `--availability-csv`: the 48-row verified listing, containing one long- and
  one short-cadence row for each pilot target;
- `--summary`: the fixed source-faithful summary at SHA-256
  `020008a34eec7993affb8bb67ff553e32b9e4c5de413800e79b152d7273b9bae`;
- `--pilot-csv`: the frozen 24-target candidate list;
- `--output-csv` and `--output-json`: new, explicit output paths.

The script refuses a changed summary hash, a cohort other than 2,123 source
rows with the four documented ESS exclusions, a kept cohort other than 2,119
planets in 1,583 targets, ambiguous target-to-KIC mappings, duplicate cadence
rows, or a pilot list that lacks source and candidate-role linkage. The summary
labels the weighting method `dynesty`; its underlying ALDERAAN source weights
are read from the `LN_WT` column.

For example, from the repository root:

```powershell
$researchRoot = $env:SAGEAR_RESEARCH_ROOT
if (-not $researchRoot) {
  throw "Set SAGEAR_RESEARCH_ROOT to the local research directory before running this audit."
}

.\.venv\Scripts\python.exe scripts\audit_mast_pilot_fits_headers.py `
  --availability-csv "$researchRoot\tmp\mast_pilot_availability_verified_20260927\mast_pilot_availability.csv" `
  --summary "$researchRoot\outputs\eccentricity_posterior_summary_SOURCE_FAITHFUL_RMS_PERIASTRON_GILBERT_QC_20260807.csv" `
  --pilot-csv metadata\population_noise_audit_20260920\pilot_candidates.csv `
  --output-csv "$researchRoot\tmp\mast_pilot_header_audit.csv" `
  --output-json "$researchRoot\tmp\mast_pilot_header_audit.json"
```

## Network bound

For every listed filename, the URL is reconstructed only from the mapped KIC
directory and the exact basename. The filename must be a safe Kepler product
name with the matching KIC and cadence suffix. The request uses
`Range: bytes=0-65535`, a streamed response, a bounded timeout, and at most
four concurrent requests by default. A server response is accepted only when
it is `206 Partial Content`, has a valid `Content-Range` beginning at zero,
and returns no more than 65,536 bytes. A `200` response, malformed range,
oversized range, short response, or range mismatch is recorded as a failure
without reading a whole product.

The tool opens only that bounded byte buffer with `astropy.io.fits`. It records
available primary `KEPLERID`, `OBSMODE`, and `QUARTER` values plus HDU 1
metadata and its `TTYPE` columns. Missing or inconsistent identity/cadence
headers are explicit `header_incomplete` records rather than silent successes.

## What this does not establish

Header presence does not establish finite flux, valid uncertainties, quality
flags, transit-window coverage, LC/SC overlap, transit timings, native masks,
or replay readiness. It does not download science payloads, modify a local
archive, request a fit, or make a population-level claim. Those questions
remain behind the provenance and residual-replay gates in
`docs/population_noise_audit_plan_20260920.md`.
