# Public-source follow-up

Date: 2026-09-24

## Question

Could a public Sagear-owned code repository or the MAST DOI cited by the paper
provide an unexamined fit configuration, catalog ledger, or posterior product
needed to complete the numerical replication?

## What was checked

The paper's acknowledgments cite MAST DOI
[`10.17909/T9059R`](https://doi.org/10.17909/T9059R). The DOI resolves to the
Kepler Input Catalog collection, rather than a project-specific archive. It is
therefore a source for public stellar-catalog inputs, not a deposit of the
paper's selected target list, fit inputs, or eccentricity posteriors.

Sagear's public GitHub account includes a fork of
[ALDERAAN](https://github.com/ssagear/alderaan). At the time of this check, the
fork has one public branch, `main`, no tags, and `main` resolves to
`7443dff16b7f9092e14a6f0cc1f8948d457c9e0b`. That is the same commit pinned by
the cloud runbooks and validation work in this repository. The public fork does
not provide a distinct paper-era branch or release to test.

The public [photoeccentric repository](https://github.com/ssagear/photoeccentric)
is useful methodological context: its tutorial creates exposure-integrated
circular transits, fits a circular geometry, and compares the inferred circular
transit density with stellar density. It uses a separate fitting workflow, so it
is appropriate as an independent diagnostic reference, not evidence of the
exact ALDERAAN configuration used for the disk-paper analysis.

The public [arXiv source archive](https://arxiv.org/abs/2509.23973) for the
disk paper was also inspected on 2026-09-24. It contains the manuscript source,
bibliography, and rendered
figures, but no catalog table, posterior archive, code, or supplementary data
directory. The source makes one important ambiguity explicit: the Methods name
Berger et al. (2018) for the stellar-density prior, while other manuscript text
uses Berger et al. (2020) for related density or stellar-property statements.
The archive does not resolve which per-target density representation entered the
eccentricity calculation.

## Consequence

These sources narrow the public-data boundary but do not remove it. They do not
contain the final accepted/rejected planet ledger, the stellar-density priors
supplied to every fit, or the per-planet eccentricity posterior
products needed to reproduce the published hierarchical likelihood exactly.

The next discriminating public-data test remains the gated circular
injection-and-recovery control described in
[`diagnostic_injection_protocol.md`](diagnostic_injection_protocol.md). It can
test the behavior of the available fitting path, but it cannot substitute for
unreleased paper-specific inputs.

## Follow-up: published host-table integrity

On 2026-09-26, the journal's machine-readable Table 1 was downloaded again
from the published article and compared byte-for-byte with
`reference/data/sagear2026_table1_kinematic_hosts_mrt.txt`. The SHA-256 is
`92280ede0c828413abb7c8314ac6f35b0ccc3e68aabbcf6a94075e84b30e76ed` for both
copies. The table contains 1,888 hosts, comprising 1,515 thin-disk and 373
thick-disk hosts; 585 have measured velocities. These values are now covered
by the parser contract tests.

This closes the kinematic-label branch of the investigation for the available
posterior archive. Every archive row with a published host match has the same
thin/thick label as Table 1 (1,775 of 1,775 rows across 1,310 hosts). An older
selection diagnostic contains reconstructed labels that disagree for some
hosts, but it is not used as the disk-label authority for the posterior
archive. The remaining replication gap cannot be attributed to a mismatch in
the published kinematic association table.
