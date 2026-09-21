"""Rank full-sample replication tensions without changing canonical membership."""

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
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--summary',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--selection-mode', choices=['manuscript_reciprocal', 'legacy_forward_norm'],
                        default='manuscript_reciprocal')
    args=parser.parse_args()
    summary=pd.read_csv(args.summary)
    validate_summary_contract(summary)
    if summary.kepoi_name.duplicated().any():
        raise ValueError('Duplicate planet IDs')
    flags=summary.qc_primary_exclude.astype(str).str.lower()
    if not flags.isin(['true','false']).all():
        raise ValueError('Invalid QC flags')
    selected=summary[flags=='false'].copy()
    args.output.mkdir(parents=True,exist_ok=True)
    sigmas=np.linspace(.002,1,600)
    means=sigmas*np.sqrt(np.pi/2)
    def estimate(ll):
        return weighted_quantile(means,posterior_weights_from_ll(sigmas,ll),[.16,.5,.84])
    rows=[]
    for disk,system,paper in [('thin','single',.022),('thick','single',.066),('thin','multi',.030),('thick','multi',.033)]:
        sub=selected[(selected.disk==disk)&(selected.system==system)].reset_index(drop=True)
        matrix,grid=load_population_masses(sub,True,args.selection_mode)
        terms=stable_terms(matrix,grid,sigmas,args.selection_mode)
        lo,median,hi=estimate(terms.sum(axis=0))
        reference=stable_terms(matrix,grid,np.array([paper,median])/np.sqrt(np.pi/2),args.selection_mode)
        score=reference[:,1]-reference[:,0]
        order=np.argsort(score)[::-1]
        ledger=sub[['kepoi_name','koi_target','disk','system','posterior_source','rho_true_solar','e50']].copy()
        ledger['log_contrast_fit_over_paper']=score
        ledger.sort_values('log_contrast_fit_over_paper',ascending=False).to_csv(args.output/f'{disk}_{system}_leverage.csv',index=False)
        removed=estimate(terms[np.setdiff1d(np.arange(len(sub)),order[:5])].sum(axis=0))[1]
        rng=np.random.default_rng(20260916)
        shifts=[]
        for trial in range(10):
            keep=np.ones(len(sub),bool)
            keep[rng.choice(len(sub),max(1,round(.1*len(sub))),replace=False)]=False
            shifts.append(float(estimate(terms[keep].sum(axis=0))[1]/median-1))
        rows.append(dict(disk=disk,system=system,n=len(sub),mean=median,lo=lo,hi=hi,paper_mean=paper,
            positive_contrast_planets=int((score>0).sum()),net_log_contrast=float(score.sum()),
            top5_net_contrast_fraction=float(score[order[:5]].sum()/score.sum()),
            diagnostic_mean_without_top5=removed,max_absolute_random_leave10pct_shift=max(abs(np.array(shifts))),
            random_leave10pct_shifts=json.dumps(shifts)))
        pd.DataFrame(rows).to_csv(args.output/'populations.csv',index=False)
        print(rows[-1],flush=True)
    (args.output/'scope.json').write_text(json.dumps(dict(summary_sha256=hashlib.sha256(args.summary.read_bytes()).hexdigest(),
        input_rows=len(summary),qc_excluded=int((flags=='true').sum()),selected_rows=len(selected),
        sigma_grid=[.002,1,600],mode=args.selection_mode,
        limitations='Log-domain diagnostic on stored posterior grids. Top-five removal is sensitivity only, not a new sample. Prior/integration assumptions remain unverified. Ten random trials do not establish exhaustive robustness.'),indent=2)+'\n')


if __name__=='__main__':
    main()
