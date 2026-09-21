"""Fit diagnostic radius controls without modifying canonical membership."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from audit_likelihood_floor import stable_terms
from hierarchical_rayleigh import (validate_summary_contract, load_population_masses,
                                  posterior_weights_from_ll, weighted_quantile)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['summary', 'inventory', 'output']:
        p.add_argument('--'+name, type=Path, required=True)
    a = p.parse_args()
    s = pd.read_csv(a.summary)
    inv = pd.read_csv(a.inventory)
    validate_summary_contract(s)
    if s.kepoi_name.duplicated().any() or inv.kepoi_name.duplicated().any():
        raise ValueError('Duplicate planet identifiers')
    flags = s.qc_primary_exclude.astype(str).str.lower()
    if not flags.isin(['true', 'false']).all():
        raise ValueError('Invalid QC flags')
    s = s.loc[flags.eq('false')].merge(inv[['kepoi_name','berger_rad','koi_prad']],
        on='kepoi_name', how='left', validate='one_to_one', indicator=True)
    if not s['_merge'].eq('both').all():
        raise ValueError('Unmatched inventory rows')
    s['radius_proxy'] = s.ror50 * s.berger_rad * 109.1
    sigmas = np.linspace(.002, 1, 600)
    rows = []
    a.output.mkdir(parents=True, exist_ok=True)
    for (disk, system), group in s.groupby(['disk', 'system']):
        group = group.reset_index(drop=True)
        matrix, grid = load_population_masses(group, True, 'manuscript_reciprocal')
        terms = stable_terms(matrix, grid, sigmas, 'manuscript_reciprocal')
        for field in ['radius_proxy', 'koi_prad']:
            radius = pd.to_numeric(group[field], errors='coerce')
            valid = np.isfinite(radius) & radius.gt(0)
            for label, mask in [('small', valid & radius.lt(3.5)),
                                ('intermediate', valid & radius.ge(3.5) & radius.le(8)),
                                ('large', valid & radius.gt(8))]:
                if not mask.any():
                    continue
                ll = terms[mask.to_numpy()].sum(axis=0)
                weights = posterior_weights_from_ll(sigmas, ll)
                lo, med, hi = weighted_quantile(sigmas*np.sqrt(np.pi/2), weights, [.16,.5,.84])
                rows.append(dict(disk=disk, system=system, radius_definition=field,
                    group=label, n=int(mask.sum()), missing_radius=int((~valid).sum()),
                    mean=med, lo=lo, hi=hi, grid_edge_mass=float(weights[0]+weights[-1])))
                group.loc[mask, ['kepoi_name', field]].to_csv(
                    a.output/f'{disk}_{system}_{field}_{label}_members.csv', index=False)
        pd.DataFrame(rows).to_csv(a.output/'fits.csv', index=False)
        print(pd.DataFrame(rows).tail(6).to_string(index=False), flush=True)
    (a.output/'scope.json').write_text(json.dumps(dict(
        inputs={str(path):hashlib.sha256(path.read_bytes()).hexdigest()
                for path in [a.summary,a.inventory]}, selected=len(s),
        sigma_grid=[.002,1,600], mode='manuscript_reciprocal',
        caveats=['Diagnostic stored-posterior fits, not new canonical results.',
                 'Radius proxy does not propagate radius uncertainty.',
                 'KOI catalog radius is a separate sensitivity definition.',
                 'Inference and selection assumptions remain under audit.']), indent=2))


if __name__ == '__main__':
    main()
