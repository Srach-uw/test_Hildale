# Analysis scripts

Run commands from the repository root. The compact publication audit works
from files committed under `reference/data/`; analyses that use light curves,
posterior FITS, or large catalogs require explicit local paths.

## Start here

```powershell
python scripts/published_sagear_audit.py
python -m pytest -q
python scripts/check_professor_release.py
```

The published disk labels are canonical. Reconstructed GMM labels are retained
only for diagnosing the original Toomre-classification mismatch.

## Pipeline stages

| Stage | Main script | Purpose |
| --- | --- | --- |
| Published sample | `published_sagear_audit.py` | Audits the article tables and headline counts. |
| Sample reconstruction | `diagnose_sample.py` | Records every catalog join and quality-cut attrition. |
| ALDERAAN preparation | `alderaan_batch.py` | Builds system-level target catalogs and runner commands. |
| Posterior extraction | `extract_eccentricity_posteriors_direct.py` | Preserves paired transit samples and applies dynesty weights. |
| Population inference | `hierarchical_rayleigh.py` | Fits Rayleigh populations with explicit selection modes. |
| Factorial validation | `compare_factorial_validation.py` | Compares cadence, limb-darkening, prior, and seed arms. |
| Weight audit | `nested_weighting_forensic_audit.py` | Tests nested-sampling weight contracts. |
| Post-fit QC | `gilbert_postfit_qc_audit.py` | Applies documented ALDERAAN quality criteria. |
| Toomre diagnostics | `toomre_diagnostics.py` | Compares velocity conventions and classifier variants. |

## Recent diagnostic tools

The September audit scripts are grouped by question rather than by pipeline
stage. `audit_factorial_cadence.py` checks whether a saved factorial arm used
short-cadence data. `build_population_noise_inventory.py` constructs a bounded
source manifest before any residual scan. `compare_alderaan_repeat_posteriors.py`
compares only compatible saved FITS systems and records mismatches instead of
forcing a comparison. `plot_population_comparison.py` generates the public
Table 3 comparison figure.

The first two tools may require large local inputs. Supply file paths on the
command line, or set `SAGEAR_RESEARCH_ROOT` when using the documented local
layout. Do not add raw light curves, large posterior archives, or personal
paths to the repository.

`prepare_diagnostic_injection_input.py` and
`diagnostic_transit_injection.py` support a separate, noncanonical
photometry-to-density control. The preparation step fails on incomplete
system identity, timing, quality-mask, cadence, or exposure provenance. The
renderer preserves those arrays and hashes but does not fit a transit or
produce an ALDERAAN posterior. Read both diagnostic protocols in `docs/`
before using either command.

## Scientific contracts

- Assign system multiplicity before removing individual planets.
- Keep `T14`, `Rp/R*`, impact, and period samples row-paired.
- Use `LN_WT` for dynesty posterior weighting.
- Record the stellar-density source and uncertainty model for every row.
- Keep manuscript-literal and forward-normalized selection modes separate.
- Reject mixed provenance and incomplete QC fields in canonical fits.
- Treat point-estimate `e_photo` values as diagnostics, not posterior data.

These contracts are enforced by the focused tests beside each analysis script.

## Local configuration

Copy `config.json` outside the repository or pass an alternate configuration
path. Do not replace committed defaults with personal paths.

```powershell
python scripts/diagnose_sample.py --config C:\path\to\local_config.json
```

Large inputs are intentionally excluded from Git. Their provenance and expected
locations are described in `docs/data_availability.md`.

## Current boundary

The reconstructed inventory contains 2,465 planets, but the strict
source-faithful hierarchy does not reproduce Table 3. The remaining blockers
are the paper's accepted planet ledger, exact stellar-density priors, final
planet posteriors, and population-model inputs. See
`docs/replication_status.md` and `docs/final_public_data_boundary.md` before
interpreting or extending the numerical results.
`build_diagnostic_alderaan_fixture.py` turns one rendered circular injection
into a separate, single-planet MAST-FITS fixture and a matching ALDERAAN
catalog. It checks immutable input hashes and source-row mappings, masks the
known real transits with NaNs, and writes cloned input files only. It does not
run ALDERAAN or write a result FITS.

`verify_diagnostic_fixture_reader.py` runs only in the pinned ALDERAAN
environment. It audits which fixture cadences ALDERAAN itself retains before
any sampling begins.

`audit_circular_density_overlap.py` is a bounded per-planet diagnostic. It
uses row-paired ALDERAAN `P`, `T14`, `ROR`, `IMPACT`, and `LN_WT` samples to
calculate a Dynesty-weighted circular transit density with the existing
MacDougall implementation. It then compares that density with either declared
stellar-density samples or explicitly declared p16/p50/p84 stellar-density
quantiles. The CSV and JSON manifest report signed log-density offsets, 68%
intervals, a configurable tolerance probability, provenance hashes, and the
configuration. It is not a significance test or a population result.

```powershell
python scripts/audit_circular_density_overlap.py `
  --samples C:\path\to\paired_samples.csv `
  --stellar-samples C:\path\to\stellar_density_samples.csv `
  --output-prefix C:\path\to\diagnostic\planet_name
```

`P` and `T14` must be in days; stellar densities must be in solar-density
units. The default stellar sample column is `rho_star_solar`. For a transparent
quantile-only diagnostic, replace `--stellar-samples` with `--stellar-p16`,
`--stellar-p50`, and `--stellar-p84`.
