# Remaining Reproducibility Questions for Sagear et al.

The final AJ article and machine-readable Table 1 resolved the kinematic host labels,
velocities, `P_thick` values, ALDERAAN release, and the apparent 378-versus-373 host
count inconsistency. The questions below are the remaining blockers to an exact
planet-level numerical replication.

1. Could you share the final planet-level analysis table containing KIC, KOI, period,
   single/multi label, and transit-fit QC inclusion? Table 1 fixes host membership and
   disk labels, but it does not identify the 2,465 individual planets or the visually
   rejected/nonconvergent fits.
2. The methods state that stellar densities are from Berger et al. (2018), while that
   catalog provides radii rather than the homogeneous density column available in
   Berger et al. (2020). Which exact table/columns or mass-radius calculation produced
   each stellar-density prior and uncertainty?
3. Could you share the final paired ALDERAAN transit-shape samples or `(e, omega)`
   posteriors, plus the list of visually rejected fits? These products determine whether
   apparently high-eccentricity leverage points are astrophysical or fit pathologies.
4. Which exact NumPyro model code and run settings were used for the hierarchical
   inference? The paper gives two chains and 1,000 steps with `Rhat < 1.05`, but not the
   warmup count, initialization, random seeds, or complete likelihood implementation.
5. Were asymmetric Berger density errors represented by a split normal, symmetrized,
   or propagated through density samples in the postmodel importance calculation?
6. Can the machine-readable ALDERAAN input catalog be shared, including limb-darkening
   prior centers and the cadence chosen for each system? This would make the remaining
   transit-prior and short-cadence factorial tests directly identifiable.

The public ALDERAAN v0.1.0 archive is code-equivalent to the commit used in the Hilldale
cloud run, so an ALDERAAN version mismatch is no longer considered a blocker.
