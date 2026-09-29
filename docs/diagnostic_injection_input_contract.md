# Diagnostic injection-input contract

`scripts/prepare_diagnostic_injection_input.py` prepares a local Kepler
PDCSAP light-curve input for `scripts/diagnostic_transit_injection.py`. It is a
separate diagnostic stage. It does not run ALDERAAN, reconstruct any unavailable
ALDERAAN code, alter a stored posterior, detrend flux, or subtract a fitted
transit model.

## Inputs and identity checks

The command requires local PDCSAP FITS files, one ALDERAAN result FITS file,
the exact local ALDERAAN system catalog used for that fit, an expected `TARGET`,
an expected KIC number, declared mask rationale, declared cadence details, and
a new or empty output directory. It makes no network request.

Every PDCSAP primary header must contain the requested positive KIC value in
`KEPLERID` or `KICID`. The ALDERAAN result primary header must contain the
requested `TARGET` and a positive `NPL`. For every planet index from zero to
`NPL - 1`, the result must contain exactly one populated `TTIMES_nn` binary
table with finite `INDEX`, `TTIME`, `MODEL`, `OUT_PROB`, and `OUT_FLAG`
values. The exclusion mask uses `MODEL`. An unexpected timing extension, a
missing timing extension, or an identity mismatch stops the run before output
is written.

The explicit local system catalog must use the ALDERAAN input-catalog fields
`koi_id`, `kic_id`, `npl`, `period`, `epoch`, and `duration`. Every row for the
requested `koi_id` must map to the requested KIC. Its row count and every
stored `npl` value must equal the result FITS `NPL`. This establishes that the
catalog presented to the stage contains the complete fitted system, rather than
only the one planet selected for injection. The manifest preserves its hash and
the target's catalog rows.

`INDEX` must contain nonnegative integers in strict order. `TTIME` and `MODEL`
must be strictly increasing. `OUT_PROB` must lie in `[0, 1]`, and `OUT_FLAG`
must be binary. These checks match the timing-readiness contract used for the
pilot and reject malformed transit-time tables before any cadence is selected.

## Cadence treatment

The tool reads the local `LIGHTCURVE` columns `TIME`, `PDCSAP_FLUX`,
`PDCSAP_FLUX_ERR`, `CADENCENO`, and `SAP_QUALITY`, plus the extension `QUARTER`
header or primary `QUARTER` header. When both headers define `QUARTER`, they
must agree. This accepts legitimate Quarter 0 files whose extension omits that
keyword. It invokes exactly the quality-mask API used by the pinned ALDERAAN
environment:

```python
lightkurve.KeplerQualityFlags.create_quality_mask(quality, bitmask="default")
```

The command requires Lightkurve 2.1.1. If that exact version is unavailable,
it fails instead of substituting a similar quality convention. After that
quality decision, it keeps only finite time, flux, and flux-error rows. It
does not rescale or otherwise modify the retained flux values. It sorts the
retained rows by time and rejects an empty or non-increasing time sequence.

The output NPZ preserves `time`, `flux`, `flux_err`, `cadence`, `quality`,
`quarter`, `source_file_index`, and `source_row_index`. `source_files` maps
each file index to the resolved source path. These fields permit every retained
cadence to be traced to a file and original row.

The PDCSAP headers must agree on the standard Kepler time convention: `TIMEUNIT`
in days, `TIMESYS=TDB`, zero `TIMEZERO`, and `BJDREF` or `BJDREFI/BJDREFF`
equal to 2454833.0. Missing, unsupported, or contradictory values stop the
run. Public ALDERAAN result FITS do not carry equivalent timing metadata, so
the command requires an explicit declaration for the `TTIMES` values: currently
`BKJD`, `d`, and `TDB` only. The declaration must match the PDCSAP convention.
It is recorded in the manifest; the tool never guesses a shared time basis.

`--cadence-mode` and `--exposure-minutes` are also required. They are checked
against the FITS headers rather than treated as labels. `OBSMODE` must normalize
to the declared long or short cadence mode in every source file. `TIMEDEL` must
be finite, positive, and identical across files to floating-point precision.
The declared exposure must match `TIMEDEL * 1440` within 0.05 minutes. The
manifest records the declared and observed values. A long-cadence preparation
therefore cannot be relabeled as short cadence later, and the renderer must use
the exposure actually represented by the input.

## Known-transit exclusion

The command reads every `TTIMES_nn` `MODEL` value and excludes a retained
cadence when its time is within or on the declared half-width of any model
time. The caller supplies `duration_days`, `duration_multiple`, and
`timing_uncertainty_margin_days`. The tool computes:

```
half_width_days = duration_days * duration_multiple + timing_uncertainty_margin_days
```

`duration_days` must be at least the longest companion duration in the supplied
full-system catalog, and `duration_multiple` must be at least one. This makes
the one system-wide exclusion policy conservative by construction. The manifest
records the three inputs and the computed half-width.

The output `known_transit_mask` marks excluded rows. `mask_index` lists the
remaining injection-eligible rows in strictly increasing order. A full
photometry-to-recovery test must exclude `known_transit_mask` from both the
injected-cadence selection and the recovery fit mask. Avoiding injection in
those rows alone is not enough. It therefore
satisfies the input contract of `diagnostic_transit_injection.py`, which
requires finite, strictly increasing `time` and `flux` arrays and a non-empty,
strictly increasing integer `mask_index`. If the declared windows remove every
retained cadence, preparation stops.

## Provenance and immutability

The output directory must be absent or empty. A successful run writes:

| File | Contents |
| --- | --- |
| `diagnostic_injection_input.npz` | Injection-contract fields and cadence provenance. |
| `diagnostic_injection_input_manifest.json` | Source file hashes, source row counts, result FITS hash, parsed timing tables, cadence counts, array hashes, quality convention, and software versions. |

The manifest states the number of raw, quality-passing, finite, retained,
known-transit-excluded, and injection-eligible cadences. The command never
replaces an occupied output directory. This prevents a later run from silently
changing the exact baseline or mask associated with an injection result.

## Example

```powershell
python scripts/prepare_diagnostic_injection_input.py `
  --pdcsap-fits q4.fits q5.fits `
  --alderaan-result-fits K00367-results.fits `
  --alderaan-system-catalog sagear_original_full_system_catalog.csv `
  --target K00367 `
  --kic 1234567 `
  --duration-days 0.15 `
  --duration-multiple 1.0 `
  --timing-uncertainty-margin-days 0.02 `
  --cadence-mode long `
  --exposure-minutes 29.4 `
  --timing-time-reference BKJD `
  --timing-timeunit d `
  --timing-timesys TDB `
  --output-dir diagnostic_input_k00367
```

The command is only a preparation step. A subsequent injection and recovery
run must retain this manifest and declare its own fitting configuration. A full
pipeline test injects into raw, pre-detrending PDCSAP flux and then applies the
declared workflow. A post-detrending injection can be useful as a control, but
it must be recorded separately and cannot stand in for the full pipeline test.
