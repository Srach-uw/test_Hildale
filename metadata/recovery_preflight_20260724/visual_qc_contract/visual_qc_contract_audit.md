# Final visual-QC contract audit

## Confirmed from the published article

- The public Table 1 is a host-level kinematic table with 1,888 stars.
- Transit fits are visually inspected after ALDERAAN fitting.
- Less than 2% are removed for quality or non-convergent posteriors.
- The final sample contains 2,465 planets.

## Reconstructed arithmetic

| disk | system | before visual QC | published | implied removals |
|---|---|---:|---:|---:|
| thin | single | 1121 | 1121 | 0 |
| thick | single | 275 | 275 | 0 |
| thin | multi | 883 | 862 | 21 |
| thick | multi | 212 | 207 | 5 |

The difference is 26/2,491 = 1.04%, consistent with the stated
"less than 2%" rejection. The count arithmetic implies no removals
from either singles category, 5 from thick multis, and 21 from thin
multis, subject to the KOI catalog epoch used in this reconstruction.

## What remains unknown

The article does not publish the 26 planet identifiers. Table 1 cannot
encode this planet-level decision because it contains one row per host.
The candidate ledger ranks plausible removals from missing posterior
products, explicit QC flags, impact parameter, signal-to-noise, and
transit count. Its selected rows are sensitivity candidates only and
must not be described as Sagear's actual rejected planets.

Diagnostic candidate set: 26 planets (5 thick, 21 thin).

## Operational decision

Fit every recoverable missing system. Do not pre-delete 26 planets
before fitting. After extraction, report the complete pre-visual branch,
the deterministic count-constrained sensitivity branch, and any
author-confirmed branch separately.
