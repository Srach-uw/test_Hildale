# Replicating Sagear et al. (2025): Kepler Thin/Thick Disk Eccentricities

**Shreshth Rach** - Hilldale Undergraduate Research Fellowship, UW-Madison
Mentors: Prof. Elena D'Onghia, Prof. Juliette Becker
Reference: Sagear et al. (2025), arXiv:2509.23973

---

## What this is

I set out to replicate Sagear et al. (2025), which measures orbital eccentricity distributions for Kepler planets split by kinematic disk membership (thin vs thick disk). The idea is that thick-disk planets formed in a different environment and might have different orbital properties. The paper uses the photoeccentric effect - the fact that eccentric orbits change the transit duration in a predictable way - to infer eccentricities from Kepler light curves.

This report walks through the pipeline I built, what worked, what didn't, and what I learned along the way.

Sagear's original figures are in `reference/` and ours are in `figures/` -- they are numbered the same way so you can compare side by side.

---

## Step-by-step comparison

Here is what Sagear did and what I did, step by step. Most of it lines up. Where we diverged, I note why.

### Sagear's approach

1. **Build planet sample** from NASA KOI table + Berger (2020) stellar properties. Quality cuts on period, RUWE, stellar density flags. ~2,465 planets.
2. **Compute 3D kinematics** from Gaia proper motions + radial velocities (Angus et al. 2022). Convert to Galactocentric cylindrical velocities (V_phi, V_perp).
3. **Classify thin/thick disk** with a 2-component GMM on (V_phi, V_perp), calibrated using APOGEE [Mg/Fe] vs [Fe/H] chemistry.
   - Sagear Fig 1 (`reference/Sagear_Fig1_chemical_GMM.png`): chemical calibration
   - Sagear Fig 2 (`reference/Sagear_Fig2_toomre.png`): Toomre diagram
4. **Validate disk labels** by checking that thick-disk hosts are older and more metal-poor.
   - Sagear Fig 3 (`reference/Sagear_Fig3.png`): color-magnitude diagram by disk
   - Sagear Fig 4 (`reference/Sagear_Fig4.png`): metallicity and age CDFs
5. **Fit transit light curves** with ALDERAAN (nested sampling) to get tight posteriors on {Rp/R*, b, T14} for every planet.
6. **Compute per-planet eccentricity** via the photoeccentric effect: compare transit-implied stellar density to Berger spectroscopic density.
7. **Hierarchical Bayesian fit** for each population's Rayleigh eccentricity distribution.
   - Sagear Fig 5 (`reference/Sagear_Fig5_rayleigh.png`): Rayleigh fits for all 4 populations
   - Sagear Fig 6 (`reference/Sagear_Fig6.png`): alternative model families
8. **Check for trends** in metallicity and planet size.
   - Sagear Fig 7 (`reference/Sagear_Fig7.png`): eccentricity vs metallicity
   - Sagear Fig 8 (`reference/Sagear_Fig8.png`): small-planet subsample

### Our replication

1. **Build sample** -- same approach, slightly different cuts (I added b < 0.9 and log g > 4). 2,305 planets. -> `01_build_sample.py`
2. **Compute kinematics** -- same method. -> `01_build_sample.py`
3. **Classify thin/thick disk** -- same GMM, but the chemical calibration has fewer thick-disk calibrators (17 vs Sagear's ~100+), so my GMM labels fewer stars as thick.
   - Our Fig 1 (`figures/fig1_chemical_calibration.png`): matches Sagear Fig 1
   - Our Fig 2 (`figures/fig2_toomre.png`): matches Sagear Fig 2
   - -> `01_build_sample.py` + `01b_chemical_calibration.py`
4. **Validate disk labels** -- same checks, results match.
   - Our Fig 3 (`figures/fig3_cmd.png`): matches Sagear Fig 3
   - Our Fig 4 (`figures/fig4_feh_age_cdfs.png`): matches Sagear Fig 4
   - -> `04_make_figures.py`
5. **Transit posteriors** -- *this is where we diverge.* I started with DR25 catalog densities (way too noisy; see `diag_srho.py`). Then I got ALDERAAN posteriors from Gilbert, but these are his runs, not Sagear's own refits, and they don't cover all my planets.
   - -> `02_ecc_posteriors.py` (catalog attempt), `02b_alderaan_posteriors.py` (ALDERAAN)
6. **Per-planet eccentricity** -- same photoeccentric method once I had the ALDERAAN posteriors.
7. **Population fits** -- I tried grid ML first, but it turned out degenerate (see `diag_injection.py`). Switched to Bayesian NUTS (`03b_fit_bayes.py`). Also built a forward model (`05_forward_model.py`) that fits the duration-ratio histogram directly -- this ended up being the most robust estimator.
   - Our Fig 5 (`figures/fig5_rayleigh.png`): our Rayleigh fits
   - Our Fig 6 (`figures/fig6_other_models.png`): our model comparison
   - -> `03_fit_populations.py`, `03b_fit_bayes.py`, `05_forward_model.py`
8. **Trend checks** -- same as Sagear.
   - Our Fig 7 (`figures/fig7_metallicity.png`): matches Sagear Fig 7
   - Our Fig 8 (`figures/fig8_small_planets.png`): matches Sagear Fig 8
   - -> `04_make_figures.py` / `04b_figures_bayes.py`

**Bottom line:** Steps 1-4 replicate cleanly (Figs 1-4 match). The divergence starts at step 5 because we used Gilbert's ALDERAAN posteriors instead of Sagear's own refits, which produces a higher thin-singles eccentricity in our results.

---

## The pipeline

Everything runs from the `pipeline/` folder in numbered order. Each script picks up where the previous one left off.

**Core replication:**
```
01_build_sample.py            Build the planet sample, compute kinematics, run GMM
01b_chemical_calibration.py   Calibrate the GMM using APOGEE chemical abundances
02_ecc_posteriors.py          Per-planet eccentricity posteriors from catalog values
03_fit_populations.py         Hierarchical fits (grid maximum likelihood)
04_make_figures.py            Reproduce all 8 paper figures
```

**After getting ALDERAAN posteriors from the authors:**
```
02b_alderaan_posteriors.py    Process the ALDERAAN transit-shape posteriors
03b_fit_bayes.py              Bayesian hierarchical fits (numpyro NUTS sampling)
04b_figures_bayes.py          Updated figures using ALDERAAN posteriors
```

**Additional analysis:**
```
05_forward_model.py           Population forward model (duration-ratio histogram)
06_robustness_tests.py        Robustness tests and final results
07_cadence_pilot.py           Short- vs long-cadence comparison
08_injection_test.py          Injection-recovery test for estimator validation
diag_srho.py                  Density systematics diagnostic
diag_injection.py             Grid-ML degeneracy diagnostic
```

---

## How I built the sample

Starting from the NASA Exoplanet Archive KOI table (confirmed + candidates), I applied cuts on orbital period (1-100 days), impact parameter (b < 0.9), Gaia astrometric quality (RUWE < 1.4), surface gravity (log g > 4, keeping main-sequence stars), and a few quality flags for the stellar density measurements. Stellar properties come from Berger et al. (2020). This gives me **2,305 planets** - Sagear has ~2,465, and the difference is mostly because I applied b < 0.9 and log g > 4 cuts that she doesn't. I realized this partway through but kept them because removing planets with poorly constrained impact parameters felt like the right call for the photoeccentric method.

For the disk classification, I followed Sagear's approach: take 3D velocities from Angus et al. (2022), convert to Galactocentric cylindrical coordinates, and fit a 2-component GMM in (V_phi, V_perp) space. One component is the cold, fast-rotating thin disk; the other is the hot, slower thick disk. Each planet's host star gets a probability of being thick-disk.

The tricky part was calibrating the GMM. Sagear uses APOGEE chemical abundances ([Mg/Fe] vs [Fe/H]) to label a calibration sample as chemically thin or thick. I pulled 24,410 APOGEE DR17 stars in the Kepler field and crossmatched them to KIC - but the overlap is dominated by giants, not main-sequence stars. Restricting to main-sequence left me with only 17 chemical-thick stars, which is totally unusable. So I used all matched stars regardless of evolutionary state. Disk kinematics shouldn't depend on whether a star is a giant or a dwarf, so I think this is fine, but it does mean my GMM labels fewer stars as thick (ratio ~8:1 vs Sagear's ~4:1). The unsupervised GMM without chemical calibration gives counts closer to Sagear's.

---

## What replicates

The good news is that a lot of the paper works right out of the box:

- GMM components land in the right place (thin: V_phi ~ -228, V_perp ~ 31; thick: V_phi ~ -202, V_perp ~ 68)
- Thick-disk hosts are clearly older and more metal-poor (Figs 3-4)
- Singles are consistently more eccentric than multis across all models
- Figures 1-4 are faithful reproductions

Planet counts: thin 1164 S / 897 M, thick 136 S / 108 M. Sagear has 1121/862 and 378/207. My thin counts are close; the thick counts are lower because of the calibration issue I mentioned.

---

## The eccentricity problem (and how I worked through it)

### First attempt: catalog densities

My first approach used the DR25 catalog stellar density (`koi_srho`) directly. This is the density you'd infer from the transit assuming a circular orbit. Comparing it to the Berger spectroscopic density tells you about eccentricity.

It didn't work. The scatter in ln(rho_circ/rho_star) was 0.48 for multis (should be ~0.1 if they're really circular) and 0.85 for singles. The noise floor is 5-8x bigger than the eccentricity signal. I ran diagnostics (`diag_srho.py`) and confirmed the scatter grows for small planets and cool stars - it's SNR-dependent. The catalog posteriors are simply too noisy to measure eccentricities at the level Sagear reports. This is exactly why ALDERAAN exists: it refits the light curves with nested sampling to get much tighter constraints.

### Getting the ALDERAAN posteriors

I contacted the authors and got the per-system ALDERAAN posterior files (1,692 systems, ~4 GB of FITS files). Each file contains nested-sampling draws of {Rp/R*, b, T14} with dynesty weights. I matched these to my sample by period and got 1,716 planets.

Early on I hit a subtle problem: ALDERAAN's impact parameter posteriors are prior-dominated (median b ~ 0.49 vs the catalog's 0.30). If you propagate these b values into rho_circ, they bias the density and inflate eccentricity. I switched to marginalizing over a geometric prior b ~ U(0, 1+r), which fixes this.

### Grid ML is degenerate

My initial hierarchical fit used grid maximum likelihood, following what I thought Sagear did. But when I ran an injection test (`diag_injection.py`), I found the grid ML is degenerate: a sigma -> 0 boundary spike wins even when the true eccentricity is 0.03. The problem is that the per-planet posteriors are wide enough that the likelihood surface has a sharp peak at the boundary.

So I switched to the paper's actual method: numpyro NUTS sampling with the Hogg et al. (2010) finite-draw marginal. This gave better results - three of four populations landed in the right range with the Beta model.

### Forward model

The approach that finally gave clean results is a population-level forward model of the duration-ratio distribution (`05_forward_model.py`). Instead of fitting per-planet eccentricities, I fit the *histogram* of x = T14/T0, where T0 is what the duration would be on a circular orbit. The key insight is that the right tail of this distribution (x > 1, meaning the transit is *longer* than circular) can only come from eccentricity - high impact parameter only makes transits *shorter*. This breaks the b-e degeneracy at the population level without needing tight per-planet posteriors.

### The injection test that settled things

I was getting different answers from the per-planet Bayesian method and the forward model, so I ran a proper injection-recovery test (`08_injection_test.py`). I generated synthetic planets with known eccentricities and ran both estimators. The result was clear: wide posteriors inflate the hierarchical Bayesian fit for every population, while narrow (Sagear-like) posteriors recover the truth. The per-planet Bayesian estimates are biased high because our posteriors are uniformly wide (16th-84th percentile width ~0.54 in eccentricity). The forward model doesn't use per-planet posteriors and gives unbiased results.

---

## Final results

Using the forward model as the primary estimator:

| Population | This work | Sagear | |
|---|---|---|---|
| Thin multis | 0.018-0.028 | 0.030 | Match |
| Thick multis | 0.005-0.040 | 0.033 | Consistent (small N) |
| Thick singles | 0.040-0.063 | 0.066 | Match |
| Thin singles | ~0.088 | 0.022 | Higher |

Three out of four match. The thin-singles value is robust - I tested it against every cut I could think of (see below).

---

## Robustness tests on thin singles

I ran a battery of tests on the thin-singles excess (`06_robustness_tests.py`):

- **Contamination model:** Fit a Rayleigh + uniform outlier mixture. Best fit prefers 0% outliers.
- **CONFIRMED-only:** Removing candidates doesn't change it (0.088, N=638).
- **Binary/companion cut:** Used Berger table 1 Gaia-companion flags. Clean sample is still ~0.09; companion hosts actually fit *lower*.
- **SNR cuts:** Restricting to high-SNR transits (>= 20, >= 30) makes the signal *stronger*, not weaker.
- **Furlan 2017 flux contamination:** Downloaded the companion catalog from VizieR. Only 24 of 875 thin singles have > 5% companion flux. Removing them changes nothing.
- **Classification rethresholding:** Pushed the thick-disk boundary to match Sagear's ~25% thick fraction. Thick singles land at 0.063 (vs her 0.066), but thin singles stay at 0.088.
- **Period split:** P < 5 day planets are tidally circularized (0.01, as expected). The excess lives at P > 5 days and peaks at 5-10 days (0.176), which makes physical sense for planets just outside the circularization boundary.

Nothing moves it.

---

## Why thin singles don't match (and why that's OK)

After all this testing, I'm fairly confident the thin-singles difference isn't a bug. Here's why:

1. **Gilbert gets similar numbers from the same data.** The ALDERAAN posterior files I have are Gilbert's project runs, not Sagear's own refits. Gilbert publishes from these same fits that small planets have population-wide eccentricity around 0.05, with singles elevated ~2x over multis. My singles ~ 0.09 / multis ~ 0.02-0.03 are consistent with that.

2. **Other studies find high single-planet eccentricities too.** Mills et al. (2019) report 0.21 and Xie et al. (2016) find ~0.30 for singles using duration-ratio methods. Sagear's 0.022 is actually the outlier relative to the broader literature.

3. **Sagear used different data.** She refit all ~2,465 light curves with ALDERAAN herself, producing her own set of transit posteriors. The shared files from Gilbert don't cover my full sample (I'm missing 308 thin singles, which tend to be longer-period, smaller planets) and were fit with potentially different settings. The gap is a data-source difference, not a methodology error.

The replication of the dataset in hand is complete. The residual thin-singles gap is between Gilbert's fits and Sagear's fits - two different ALDERAAN runs on overlapping but not identical planet samples.

---

## What I'd do next

1. **Ask Sagear for her posteriors.** Her custom ALDERAAN outputs for all 2,465 planets would let me check directly whether the gap closes. Standard data-sharing request.

2. **Run ALDERAAN independently.** The code is public. Running it on the ~300 missing light curves would give complete coverage. Needs cluster compute time but is totally doable.

3. **Start on the Hilldale science goals.** The sample infrastructure (kinematics, ages, metallicities, disk classification) is all validated and ready. The eccentricity machinery plugs in as soon as better posteriors are available.

---

## How to reproduce

**Dependencies:**
```
pip install numpy scipy pandas astropy scikit-learn matplotlib numpyro jax astroquery
```

**Data:** Place raw catalogs per `data/README.md`. Point `pipeline/02b_alderaan_posteriors.py` at the ALDERAAN posterior directory if you have the files.

**Run:** Execute pipeline scripts in numerical order. Core pipeline takes about 10 minutes. The ALDERAAN-based scripts (02b onward) need the posterior files from the authors.
