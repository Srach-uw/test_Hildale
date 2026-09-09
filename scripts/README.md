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
