# Corrected Cadence Validation

## Why Repeat This Comparison

The archived factorial's 18 LC+SC treatments did not enable short cadence.
All their logs state that no short-cadence data were processed. The local
runners now pass `--use_sc True` when `CADENCE_MODE=both`; old outputs remain
unchanged. This repair needs a real-data validation, not a blanket rerun.

## Target Selection

Use the current 2465-planet inventory for disk and multiplicity labels. The
old factorial target list differs: K01299 is now thick rather than thin;
K00283 and K02533 are absent from this inventory. Do not treat the old list
as a representative sample of the current four populations.

Start with K00319 as the historical low-disagreement control and K02712 as
the high-impact-parameter thin-single diagnostic. Their old inventory lists
14 and four SC files respectively, but this is not a verified count of usable
transits. Add K00972 as a separate multi-system pathology investigation after
the control succeeds. Do not pool these selected systems into a population
estimate. Their role is to establish treatment fidelity and diagnose geometry.

## Before Sampling

1. Fetch or locate original PDCSAP light curves and record file hashes,
   KIC identity, cadence, quarters, finite measurements and transit coverage.
2. Stage the pinned ALDERAAN environment and the same full-system catalog,
   limb-darkening values and seed for each paired comparison.
3. Use fresh project directories and run IDs for LC-only and SC-enabled arms.
   Never reuse completed-result skip logic from the old experiment.
4. Apply the runtime cadence instrumentation to the staged transit driver;
   keep the before/after source hashes and all scientific settings.
5. Run the processed-cadence verifier after detrending and after filtering.
   In the SC arm set `ALDERAAN_REQUIRE_SC=1` so absent selected SC data abort
   before nested sampling. Do not impose this requirement on ordinary LC fits.
6. Inspect the runtime per-quarter and per-planet counts. A single SC quarter
   in one planet does not establish adequate SC coverage for every companion.

## Interpretation

Compare paired posterior duration, impact, radius ratio and circular-density
distributions using nested weights. Then use identical stellar-density and
eccentricity extraction settings. Inspect transit residuals and sampler
diagnostics, not only posterior medians. Preserve the full-system labels.

SC replaces LC within available quarters in the upstream loader; the arms
therefore may differ in time coverage as well as integration time. Record
that explicitly and consider a matched-window LC control before attributing
all changes to temporal resolution. A null effect on two selected systems
does not close cadence for the full sample; a large effect does not establish
the paper's population result either.

## Execution Status

Bounded local validation has progressed for K02712: the task-workspace
manifests verify downloaded raw short- and long-cadence light curves for Q3
and Q14, and isolated conditional fits plus the native replay were completed.
This is not a full corrected ALDERAAN run. Do not assume old cloud credits
remain valid or launch paid compute without a current cost and account check.
No canonical posterior merge is authorized by this plan.
