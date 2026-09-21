"""Conditional BATMAN transit fits with supplied versus measured quarterly errors.

This is a geometry diagnostic, not a replacement posterior or independent
detrending. Limb darkening and regularized transit times remain fixed.
"""

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from astropy.io import fits
from scipy.optimize import least_squares

from audit_transit_timing_fold import nearest_phase
from duration_measure_control import density_and_derivative


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['directory','noise','output','dependency-path']:
        parser.add_argument('--'+name,type=Path,required=True)
    parser.add_argument('--supersample',type=int,default=11)
    parser.add_argument('--profile',action='store_true')
    parser.add_argument('--diff-step',type=float,default=None)
    args = parser.parse_args()
    sys.path.insert(0,str(args.dependency_path))
    import batman
    paths = [args.directory/'K02712_lc_detrended.fits',
             args.directory/'K02712_00_quick.ttvs',
             args.directory/'K02712_transit_parameters.csv', args.noise]
    with fits.open(paths[0]) as hdul:
        time,flux,error,quarter = [np.asarray(hdul[k].data) for k in ['TIME','FLUX','ERROR','QUARTER']]
    times = np.loadtxt(paths[1])
    period = float(np.polyfit(times[:,0],times[:,2],1)[0])
    phase = nearest_phase(time,times[:,2])/24
    use = np.isfinite(time+flux+error)&(error>0)&(abs(phase)<.4)
    phase,flux,error,quarter = [v[use] for v in [phase,flux,error,quarter]]
    pars = pd.read_csv(paths[2]).iloc[0]
    noise = pd.read_csv(args.noise)
    factors = noise[noise.stage=='detrended'].set_index('quarter').normalized_adjacent_scale
    scales = np.array([factors.loc[q] for q in quarter])
    p = batman.TransitParams()
    p.t0=0.; p.per=period; p.rp=.014; p.a=8.; p.inc=89.; p.ecc=0.; p.w=90.
    p.u=[float(pars.limbdark_1),float(pars.limbdark_2)]; p.limb_dark='quadratic'
    model = batman.TransitModel(p,phase,supersample_factor=args.supersample,exp_time=.02043365)
    groups = [quarter==q for q in np.unique(quarter)]

    def evaluate(x):
        radius, fraction, hours, offset = x
        impact = fraction*(1+radius)
        a2 = ((1+radius)**2-impact**2)/np.sin(np.pi*(hours/24)/period)**2+impact**2
        p.rp=radius; p.a=np.sqrt(a2); p.inc=np.degrees(np.arccos(impact/p.a)); p.t0=offset/24
        return model.light_curve(p)

    rows = []
    for mode,sigma in [('supplied',error),('quarter_adjacent_scale',error*scales)]:
        def residual(x):
            raw = evaluate(x)-flux
            for group in groups:
                raw[group] -= np.average(raw[group],weights=1/sigma[group]**2)
            return raw/sigma
        for start in [.2,.8,.96]:
            fit=least_squares(residual,[.014,start/1.014,6.,0.],
                bounds=([.003,0.,2.,-.5],[.05,.999,10.,.5]),
                diff_step=args.diff_step,x_scale='jac',max_nfev=400,ftol=1e-10,xtol=1e-10,gtol=1e-8)
            radius,fraction,hours,offset=fit.x
            impact=fraction*(1+radius)
            rho,_=density_and_derivative(hours/24,0,0,period,radius,impact)
            rows.append(dict(mode=mode,start_impact=start,success=bool(fit.success),
                evaluations=fit.nfev,chi_square=float(2*fit.cost),impact=impact,
                radius_ratio=radius,duration_hours=hours,time_offset_hours=offset,
                rho_circ_solar=float(rho)))
        if args.profile:
            for impact in [0.,.2,.4,.6,.8,.9,.94,.95,.96,.97,.98]:
                def fixed_residual(x):
                    return residual([x[0],impact/(1+x[0]),x[1],x[2]])
                candidates=[]
                for duration in [5.6,6.3]:
                    fit=least_squares(fixed_residual,[.014,duration,-.025],
                        bounds=([.003,2.,-.5],[.05,10.,.5]),diff_step=args.diff_step,
                        x_scale='jac',max_nfev=300,ftol=1e-10,xtol=1e-10,gtol=1e-8)
                    candidates.append(fit)
                fit=min(candidates,key=lambda f:f.cost)
                radius,hours,offset=fit.x
                rho,_=density_and_derivative(hours/24,0,0,period,radius,impact)
                rows.append(dict(mode=mode,start_impact=np.nan,profile_fixed_impact=True,
                    success=bool(fit.success),evaluations=fit.nfev,chi_square=float(2*fit.cost),
                    impact=impact,radius_ratio=radius,duration_hours=hours,
                    time_offset_hours=offset,rho_circ_solar=float(rho)))
    args.output.mkdir(parents=True,exist_ok=True)
    result=pd.DataFrame(rows)
    result.to_csv(args.output/'fits.csv',index=False)
    provenance=dict(inputs={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        photometry_count=len(flux),quarter_count=len(groups),period_days=period,
        batman_version=batman.__version__,exposure_days=.02043365,supersample_factor=args.supersample,
        finite_difference_step=args.diff_step,profile=args.profile,
        scope='Conditional optimization: fixed LD and saved TTVs, profiled quarter offsets, no GP, no density prior. Not a posterior or significance test; chi-square not comparable across error models.')
    (args.output/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    print(result.to_string(index=False))


if __name__=='__main__':
    main()
