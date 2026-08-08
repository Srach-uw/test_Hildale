# Hildale Sagear Replication

Code, compact evidence tables, and validation data for an independent
replication of Sagear et al. (2026), *The Orbital Eccentricities of Planets in
the Kinematic Thin and Thick Galactic Disks*.

Reference article:
[doi:10.3847/1538-3881/ae71bf](https://doi.org/10.3847/1538-3881/ae71bf)

## Status

This repository does not reproduce the population eccentricities reported in
the paper. It recovers a 2,465-planet public-data inventory, preserves the
published host classifications and pre-cut multiplicity labels, and constructs
a uniform posterior set. The total matches the paper, but the four population
counts do not, so the final planet membership is not exact.

The valid public-data branch keeps paired ALDERAAN transit-shape samples and
uses dynesty `LN_WT` weights. Its Rayleigh means are 0.234, 0.246, 0.141, and
0.171 for thick singles, thin singles, thick multis, and thin multis. These are
well above the paper's 0.066, 0.022, 0.033, and 0.030. Equal raw-row weighting
moves some values toward the paper but reverses the reported model comparison:
it favors a half-Gaussian rather than Rayleigh. That branch is retained only as
a diagnostic.

An independent control using Gilbert et al.'s released small-planet catalog
and the same real ALDERAAN archive recovers a Beta mean of 0.0487 for all small
planets and 0.0684 for observed singles, consistent with the corresponding
published values near 0.05 and 0.073. This shows that the extraction and
hierarchical machinery can recover an external published result; it does not
identify Sagear et al.'s unreleased analysis choices.

All values in this repository are replication diagnostics, not a new
astrophysical measurement.

See [docs/replication_status.md](docs/replication_status.md) for the scientific
interpretation and remaining information gaps. Compact evidence for the current
reconstruction is in
[`metadata/public_reconstruction_20260727/`](metadata/public_reconstruction_20260727/),
with the final control and weighting audits in
[`metadata/final_forensic_20260808/`](metadata/final_forensic_20260808/).
The [documentation guide](docs/README.md) separates current conclusions from
historical audits and cloud runbooks.

## Main Findings

1. The machine-readable host table contains 1,888 stars: 1,515 labeled thin and
   373 labeled thick. The text's value of 378 thick hosts is inconsistent with
   the total and subgroup arithmetic.
2. Multiplicity must be assigned from the full eligible planet inventory before
   planet-level quality cuts. A superseded 2,474-row reconstruction moved 45
   planets when this rule was corrected; that count is historical evidence,
   not a claim about the final 2,465-row inventory.
3. ALDERAAN duration, radius-ratio, impact-parameter, and period samples must
   remain paired. The extractor now enforces this contract.
4. Dynesty nested points must be weighted with `LN_WT`. Equal weighting partly
   approaches the paper's means but fails its model ordering, so the numerical
   resemblance is not accepted as a replication.
5. The completed 82-fit validation matrix and a separate nine-system combined
   confirmation show that the tested cadence, limb-darkening, printed-prior,
   and sampler-seed choices do not explain the population-wide discrepancy.
   In the combined arm, the median paired eccentricity shift is `+0.00042`
   across 13 planets.
6. A direct audit of all 2,465 reconstructed planets finds median absolute
   circular-density disagreements of 0.159-0.225 dex across the four
   populations. The corresponding posterior widths are 0.611-0.712 dex.
   This localizes the mismatch upstream of the population fit without yet
   identifying a single cause.
7. Broad impact-parameter posteriors are common and broaden the inferred
   circular density. High inferred eccentricity also persists in the
   exploratory narrower-impact subset: after removing five primary-QC
   failures, its 114 thin singles give a forward-normalized hierarchical mean
   of 0.333 (0.317-0.350), compared with 0.022 in the paper. Because selection
   uses a posterior property, this is a diagnostic rather than an unbiased
   population estimate.
8. Applying Gilbert's full post-fit quality cuts to the Sagear sample does not
   recover Table 3. The same cuts do recover Gilbert's published small-planet
   result from the real ALDERAAN archive, providing a positive control.
9. The remaining high-value unknowns are the exact stellar-density and radius
   products, the final rejected-planet list, the intermediate `(e, omega)`
   posteriors, the nested-point resampling convention, and the population-fit
   implementation used for Table 3.

## Reproduce the Checks

Use Python 3.11:

```bash
python -m venv .venv
python -m pip install -r requirements-lock.txt
python -m pytest -q
python scripts/check_professor_release.py
git diff --check
```

ALDERAAN uses a separate pinned environment and repository commit. Cloud
execution is optional and billable. Read
[`docs/gcp_no_charge_safety_checklist.md`](docs/gcp_no_charge_safety_checklist.md)
before creating a VM.

## Repository Layout

| Path | Contents |
| --- | --- |
| `scripts/` | Sample construction, posterior extraction, hierarchy, diagnostics, and tests |
| `cloud/` | Reproducible ALDERAAN execution and recovery helpers |
| `metadata/` | Compact count audits, fit summaries, and provenance tables |
| `data/alderaan_factorial_validation_20260715/` | 82 validation FITS and manifests stored with Git LFS |
| `reference/` | Published article, tables, and source references |
| `docs/` | Current status, focused audits, and runbooks |
| `legacy/` | Superseded early analysis retained for provenance |

Large and regenerable products are excluded: Kepler light curves, posterior
archives, cloud result bundles, virtual environments, checkpoints, temporary
directories, and private chronological worklogs.

Retrieve the validation FITS after cloning:

```bash
git lfs install
git lfs pull
```

## Citation

Citation metadata for this repository and the source article are provided in
[`CITATION.cff`](CITATION.cff). GitHub can render this file as APA or BibTeX.

## Scientific Boundary

The repository supports reproducibility and diagnosis. A match obtained by
relabeling populations, removing influential planets after looking at the
answer, or ignoring nested-sampling weights is not accepted as a replication.
The reported discrepancy remains open pending author products or a new,
falsifiable methodological explanation.
