# Published-Host Recovery Preflight

Updated: 2026-07-24

The published Table 1 host inventory, rather than the earlier working sample, is now the
reference universe for recovery. It identifies 142 host systems containing 174 population
planets with no usable local ALDERAAN result. Of these systems, 116 were absent from the
earlier 592-target launch and 26 were launched but did not yield a usable result.

The discrepancy cannot yet be assigned to a published-method error. Missing coverage is
not random: the uncovered systems have lower median catalog surface gravity, more
evolutionary-state flags, larger radii, and lower density metadata than covered systems.
This is especially relevant to the thin-single population.

The recovery candidate is deliberately separated from the 82-fit factorial validation.
The factorial matrix tested plausible implementation differences and found effects too
small to explain the population mismatch. The remaining 142 systems must be recovered or
their exclusion must be documented before a population-level replication can be called
complete.

Supporting tables, the recovery feasibility manifest, the published-table provenance
audit, and the independent Toomre reconstruction are in
`metadata/recovery_preflight_20260724/`. The target list is
`missing_recovery_feasibility_targets.csv`; it is a manifest for controlled cloud
execution, not a set of inferred eccentricities.
