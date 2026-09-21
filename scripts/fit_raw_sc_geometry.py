"""Conditional raw-PDCSAP geometry check with event-specific linear baselines."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
import pandas as pd
from astropy.io import fits
from scipy.optimize import least_squares
from duration_measure_control import density_and_derivative


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['raw-dir','timings','parameters','dependency-path','output']:
        p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--fixed-windows',action='store_true',help='Select and normalize identical points around regularized centers')
    p.add_argument('--offset-starts',type=float,nargs='+',default=[0.])
    p.add_argument('--cadence',choices=['short','long'],default='short')
    p.add_argument('--coverage-reference',type=Path,help='SC directory defining common time coverage and eligible events')
    p.add_argument('--processed-file',type=Path,help='Use archived TIME/FLUX/ERROR arrays instead of raw files')
    a=p.parse_args(); sys.path.insert(0,str(a.dependency_path))
    import batman
    tts=np.loadtxt(a.timings,ndmin=2)
    period=float(np.polyfit(tts[:,0],tts[:,2],1)[0])
    all_t=[]; all_f=[]; all_e=[]; hashes={}
    suffix='*_slc.fits' if a.cadence=='short' else '*_llc.fits'
    exposure=58.85 if a.cadence=='short' else 1765.5
    for path in ([] if a.processed_file else sorted(a.raw_dir.glob(suffix))):
        hashes[str(path)]=hashlib.sha256(path.read_bytes()).hexdigest()
        with fits.open(path) as h:
            d=h[1].data
            use=np.isfinite(d['TIME'])&np.isfinite(d['PDCSAP_FLUX'])&np.isfinite(d['PDCSAP_FLUX_ERR'])&(d['PDCSAP_FLUX_ERR']>0)&(d['SAP_QUALITY']==0)
            all_t.extend(d['TIME'][use]); all_f.extend(d['PDCSAP_FLUX'][use]); all_e.extend(d['PDCSAP_FLUX_ERR'][use])
    if a.processed_file:
        hashes[str(a.processed_file)]=hashlib.sha256(a.processed_file.read_bytes()).hexdigest()
        with fits.open(a.processed_file) as h:
            t,f,e=[np.asarray(h[key].data).ravel() for key in ['TIME','FLUX','ERROR']]
            use=np.isfinite(t)&np.isfinite(f)&np.isfinite(e)&(e>0)
            all_t,all_f,all_e=t[use],f[use],e[use]
    time,flux,error=map(np.asarray,[all_t,all_f,all_e])
    reference_time=None
    if a.coverage_reference:
        reference=[]
        for path in sorted(a.coverage_reference.glob('*_slc.fits')):
            hashes[str(path)]=hashlib.sha256(path.read_bytes()).hexdigest()
            with fits.open(path) as h:
                d=h[1].data
                use=np.isfinite(d['TIME'])&np.isfinite(d['PDCSAP_FLUX'])&np.isfinite(d['PDCSAP_FLUX_ERR'])&(d['PDCSAP_FLUX_ERR']>0)&(d['SAP_QUALITY']==0)
                reference.extend(d['TIME'][use])
        reference_time=np.sort(np.asarray(reference))
        if len(reference_time)<100:
            raise ValueError('Insufficient reference coverage')
        blocks=np.split(reference_time,np.where(np.diff(reference_time)>.01)[0]+1)
        covered=np.zeros(len(time),dtype=bool)
        half=exposure/172800
        for block in blocks:
            covered|=(time-half>=block[0]-58.85/172800)&(time+half<=block[-1]+58.85/172800)
        time,flux,error=time[covered],flux[covered],error[covered]
    pars=pd.read_csv(a.parameters).iloc[0]
    rows=[]; selection_hashes={}
    for label,column in [('regularized',2),('independent',1)]:
        phases=[]; values=[]; errors=[]; event_ids=[]; selected_times=[]
        for event,tt in enumerate(tts[:,column]):
            if reference_time is not None:
                ref_delta=reference_time-tts[event,2]
                if (abs(ref_delta)<.4).sum()<100 or (abs(ref_delta)<.12).sum()<50:
                    continue
            delta=time-tt
            selection_delta=time-tts[event,2] if a.fixed_windows else delta
            use=abs(selection_delta)<.4
            minimum,central=(100,50) if a.cadence=='short' else (12,4)
            if use.sum()<minimum or (abs(selection_delta[use])<.12).sum()<central:
                continue
            scale=np.median(flux[use])
            phases.extend(delta[use]); values.extend(flux[use]/scale)
            errors.extend(error[use]/scale); event_ids.extend([event]*int(use.sum()))
            selected_times.extend(time[use])
        phase,y,err,event_ids=map(np.asarray,[phases,values,errors,event_ids])
        if len(y)==0: raise ValueError('No covered transits')
        selection_hashes[label]=hashlib.sha256(np.column_stack([selected_times,y,err,event_ids]).astype('<f8').tobytes()).hexdigest()
        groups=[np.where(event_ids==j)[0] for j in np.unique(event_ids)]
        designs=[np.column_stack([np.ones(len(idx)),phase[idx]])/err[idx,None] for idx in groups]
        inverses=[np.linalg.pinv(d) for d in designs]
        tp=batman.TransitParams(); tp.t0=0.; tp.per=period; tp.rp=.014; tp.a=8.; tp.inc=89.; tp.ecc=0.; tp.w=90.
        tp.u=[float(pars.limbdark_1),float(pars.limbdark_2)]; tp.limb_dark='quadratic'
        model=batman.TransitModel(tp,phase,supersample_factor=3 if a.cadence=='short' else 31,exp_time=exposure/86400)
        for impact in [0.,.5,.8,.96]:
            def residual(x):
                r,hours,offset=x
                a2=((1+r)**2-impact**2)/np.sin(np.pi*hours/24/period)**2+impact**2
                tp.rp=r; tp.a=np.sqrt(a2); tp.inc=np.degrees(np.arccos(impact/tp.a)); tp.t0=offset/24
                res=(model.light_curve(tp)-y)/err
                for idx,d,inv in zip(groups,designs,inverses):
                    res[idx]-=d@(inv@res[idx])
                return res
            candidates=[least_squares(residual,[.014,h,offset_start],bounds=([.003,2.,-.5],[.05,10.,.5]),
                diff_step=1e-4,max_nfev=200,x_scale='jac') for h in [5.6,6.3] for offset_start in a.offset_starts]
            fit=min(candidates,key=lambda x:x.cost)
            r,h,offset=fit.x
            rho,_=density_and_derivative(h/24,0,0,period,r,impact)
            rows.append(dict(timing=label,fixed_impact=impact,radius_ratio=r,duration_hours=h,
                offset_hours=offset,chi_square=2*fit.cost,success=bool(fit.success),
                rho_circ_solar=float(rho),points=len(y),events=len(groups),event_indices=','.join(str(int(tts[j,0])) for j in np.unique(event_ids)),
                starts=len(candidates),converged_starts=sum(bool(c.success) for c in candidates),
                starting_solution_chisq_range=float(np.ptp([2*c.cost for c in candidates]))))
            print(rows[-1],flush=True)
    if a.fixed_windows and len(set(selection_hashes.values()))!=1:
        raise ValueError('Fixed-window measurement inputs differ between timing modes')
    a.output.mkdir(parents=True,exist_ok=True)
    pd.DataFrame(rows).to_csv(a.output/'fits.csv',index=False)
    for path in [a.timings,a.parameters]: hashes[str(path)]=hashlib.sha256(path.read_bytes()).hexdigest()
    (a.output/'scope.json').write_text(json.dumps(dict(hashes=hashes,period_days=period,
        batman_version=batman.__version__,cadence=a.cadence,exposure_seconds=exposure,coverage_reference=str(a.coverage_reference),fixed_windows=a.fixed_windows,measurement_hashes=selection_hashes,offset_starts=a.offset_starts,
        processed_file=str(a.processed_file),
        scope='Fixed impact grid, fixed LD, saved timings; raw PDCSAP quality zero or named processed file, per-event linear baselines, Gaussian supplied errors. Not ALDERAAN, a posterior, or a chi-square significance test.'),indent=2))


if __name__=='__main__': main()
