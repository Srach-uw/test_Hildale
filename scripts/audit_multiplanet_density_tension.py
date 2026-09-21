"""Audit period matching and paired geometry within influential multi systems."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits
from audit_influential_transit_shapes import quantiles
from extract_eccentricity_posteriors_direct import _sample_frame, nested_sample_weights, paired_period_samples
from duration_measure_control import density_and_derivative


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for key in ['summary','catalog','output']:
        parser.add_argument('--'+key,type=Path,required=True)
    args=parser.parse_args()
    summary=pd.read_csv(args.summary)
    catalog=pd.read_csv(args.catalog,comment='#')
    selected=summary[summary.koi_target.isin(['K01240','K00972'])].copy()
    selected=selected.merge(catalog[['kepoi_name','koi_duration','koi_impact','koi_disposition','koi_model_snr']],on='kepoi_name',validate='one_to_one')
    rows=[]; hashes={}
    for _,row in selected.iterrows():
        path=Path(row.raw_result_file)
        hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
        idx=int(row.alderaan_planet_index)
        with fits.open(path) as hdul:
            samples=_sample_frame(hdul['SAMPLES'].data)
            period=paired_period_samples(samples,hdul,idx)
        weight=nested_sample_weights(samples.LN_WT.to_numpy(float),'dynesty')
        t,r,b=[samples[f'{key}_{idx}'].to_numpy(float) for key in ['DUR14','ROR','IMPACT']]
        rho,_=density_and_derivative(t,0,0,period,r,b)
        valid=np.isfinite(rho)&(rho>0)&(t>0)&(b>=0)&(b<1+r)
        out={key:row[key] for key in ['kepoi_name','koi_target','koi_period','alderaan_planet_index','period_relative_difference','rho_true_solar','e50','qc_primary_exclude','qc_reasons','koi_duration','koi_impact','koi_disposition','koi_model_snr']}
        for name,values in [('duration_hours',t*24),('impact',b),('ror',r),('rho_circ_solar',rho)]:
            for label,value in zip(['p16','p50','p84'],quantiles(values[valid],weight[valid])):
                out[name+'_'+label]=float(value)
        ratio=out['rho_circ_solar_p50']/row.rho_true_solar
        out['density_ratio']=ratio
        out['central_density_emin_approx']=abs(ratio**(2/3)-1)/(ratio**(2/3)+1)
        out['duration_ratio_to_catalog']=out['duration_hours_p50']/row.koi_duration
        out['valid_weight']=float(weight[valid].sum())
        rows.append(out)
    result=pd.DataFrame(rows)
    crossing=[]
    for target,sub in result.groupby('koi_target'):
        sub=sub.sort_values('koi_period')
        inner,outer=sub.iloc[0],sub.iloc[-1]
        aratio=(inner.koi_period/outer.koi_period)**(2/3)
        crossing.append(dict(target=target,inner=inner.kepoi_name,outer=outer.kepoi_name,
            inner_to_outer_semimajor_ratio=aratio,
            outer_e_for_radial_overlap_circular_inner=1-aratio,
            outer_central_density_emin_approx=outer.central_density_emin_approx))
    args.output.mkdir(parents=True,exist_ok=True)
    result.to_csv(args.output/'paired_geometry.csv',index=False)
    pd.DataFrame(crossing).to_csv(args.output/'radial_overlap_diagnostic.csv',index=False)
    (args.output/'scope.json').write_text(json.dumps(dict(fits_sha256=hashes,
        input_sha256={key:hashlib.sha256(getattr(args,key).read_bytes()).hexdigest() for key in ['summary','catalog']},
        caveats='Minimum eccentricities use central density ratios and the short-transit approximation; not posterior bounds. Radial overlap is not proof of dynamical instability or a justified exclusion; orientation, masses, resonances and density uncertainties are not modeled.'),indent=2)+'\n')
    print(result[['kepoi_name','koi_duration','duration_hours_p50','koi_impact','impact_p50','density_ratio','central_density_emin_approx']].to_string(index=False))
    print(pd.DataFrame(crossing).to_string(index=False))


if __name__=='__main__':
    main()
