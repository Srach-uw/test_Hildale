"""Check selection conventions with an exactly constant measurement likelihood."""
import argparse
import json
from pathlib import Path
import numpy as np
from audit_likelihood_floor import stable_terms


def null_contrasts():
    e = np.linspace(0, .6, 801)
    omega = np.linspace(0, 2*np.pi, 360, endpoint=False)
    widths = np.gradient(e)
    widths[[0,-1]] *= .5
    posterior = np.broadcast_to(widths[:,None], (len(e),len(omega))).copy()
    posterior /= posterior.sum()
    transit = (1+e[:,None]*np.sin(omega))/(1-e[:,None]**2)
    sigmas = np.linspace(.002,.3,300)
    matrices = {'none':posterior.sum(axis=1),
                'legacy_forward_norm':(posterior*transit).sum(axis=1),
                'manuscript_reciprocal':(posterior/transit).sum(axis=1)}
    return {mode: stable_terms(mass[None,:],e,sigmas,mode)[0]
            for mode,mass in matrices.items()}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    curves=null_contrasts()
    result={mode:dict(log_likelihood_range=float(np.ptp(curve)),
                     hot_minus_cold_per_planet=float(curve[-1]-curve[0]),
                     hot_minus_cold_887_planets=float(887*(curve[-1]-curve[0])))
            for mode,curve in curves.items()}
    result['scope']='Exactly constant likelihood on uniform e/omega interim measure; geometric-transit-only selected generative model. Not a test of ALDERAAN or the paper implementation.'
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
