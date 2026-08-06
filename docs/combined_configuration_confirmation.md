# Combined-Configuration Confirmation Arm

## Why It Exists

The completed 82-fit experiment is a fractional factorial design. It measures
limb darkening, cadence, the printed-prior ambiguity, and sampler variability
in separate paired arms. It does not contain the exact combined configuration
closest to the written methods: reference limb-darkening centers, available
long and short cadence, and the Table 1 transit-prior sensitivity together.

The full matrix found no panel-wide effect large enough to reconcile the
current population mismatch. It also found a small number of systems with
large individual cadence responses. A compact confirmation arm is therefore
appropriate before treating the three choices as independently additive.

## Scope

The arm runs the nine short-cadence-audited systems already present in the
complete matrix:

`K00064`, `K00283`, `K00319`, `K00680`, `K00716`, `K01001`, `K01299`,
`K02533`, and `K02712`.

This includes `K00283` and `K02533`, the systems with the largest cadence
responses. The arm fits each full seedable KOI system simultaneously, just as
the completed validation arms did.

## Fixed Inputs

| component | value |
|---|---|
| Limb darkening | `sagear_ld_reference_catalog.csv` |
| Cadence | `both`, long cadence plus available short cadence |
| Transit priors | `patch_alderaan_paper_priors.py` sensitivity clone |
| Target set | `targets_short_cadence_validation.csv` |
| Public ALDERAAN base | pinned commit `7443dff16b7f9092e14a6f0cc1f8948d457c9e0b` |
| Run identifier | `sagear_validation_paper_priors_reference_lcsc` |

## Execution

This is a new, separate run and must not be copied into the immutable 82-fit
release. From a VM containing the validated `cloud/ld_validation` bundle:

```bash
cd ~/sagear_ld_validation_batch
source ~/miniforge3/etc/profile.d/conda.sh
conda activate alderaan
JOBS=6 nohup bash run_combined_confirmation.sh > combined_confirmation.log 2>&1 &
```

The runner is resumable. It writes results under
`projects/paper_priors_reference_lcsc/` and skips an already completed target.

## Result

Compare this arm only with the already completed `reference_lcsc` arm for the
same planets and systems. This contrast estimates the paper-prior effect
conditional on reference limb darkening and LC+SC. It is not, by itself, a
factorial difference-in-differences interaction estimate because the available
paper-prior baseline arm uses a different eight-target set.

All nine systems completed with exit code 0. The matched comparison contains
13 planets and no direct-posterior QC exclusions.

| parameter | median signed change | median absolute change |
|---|---:|---:|
| eccentricity | +0.00042 | 0.01135 |
| zeta | -0.00019 | 0.00982 |
| impact parameter | -0.00144 | 0.00740 |
| transit duration (hr) | -0.00454 | 0.01110 |
| radius ratio | +0.000017 | 0.000037 |

Five of 13 eccentricity shifts exceed the 95th percentile of the available
repeat-run shifts, showing that the prior patch can matter for individual
difficult systems. The shifts are not coherent in sign, however, and their
median is near zero. This combined configuration therefore does not explain
the population-wide discrepancy.

The conclusion remains conditional on this targeted system set. This is a
narrow configuration check, not a replacement population study. It also does
not substitute for the exact Berger et al. (2018) density construction or
Sagear's final visual-QC inclusion list.

The immutable FITS release is
`data/alderaan_combined_confirmation_20260806/`; compact regenerated results are
in `metadata/combined_confirmation_20260806/`.

Regenerate the comparison from the two immutable releases:

```powershell
python scripts/compare_factorial_validation.py `
  --validation-root data/alderaan_factorial_validation_20260715/results `
  --combined-validation-root data/alderaan_combined_confirmation_20260806 `
  --include-combined-confirmation `
  --metadata-root data/alderaan_factorial_validation_20260715/provenance/target_sets `
  --inventory data/alderaan_factorial_validation_20260715/provenance/input_catalogs/full_system_inventory.csv `
  --run-contract metadata/factorial_validation_20260715/arm_run_contract.csv `
  --sample <live-sagear-reproduction>/outputs/canonical_sample_old_astropy_rawcc.csv `
  --config <live-sagear-reproduction>/config.json `
  --output-dir <output-directory> `
  --n-proposals 150000 --bootstrap-replicates 10000 --seed 20260715
```
