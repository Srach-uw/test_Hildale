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

## Acceptance Rule

Compare this arm only with the already completed `reference_lcsc` arm for the
same planets and systems. The result should be interpreted as follows:

1. If the combined arm is close to `reference_lcsc` at the affected systems,
   the printed-prior interaction is not a credible explanation for the
   population discrepancy.
2. If it changes the affected systems materially, record the result as a
   target-level posterior-construction sensitivity and repeat the comparison
   on independent controls before changing any full-catalog configuration.
3. Neither outcome substitutes for the missing Berger et al. (2018) density
   construction or Sagear's final visual-QC inclusion list. Those remain the
   gating provenance questions for an exact Table 2 replication.
