"""Check a log-duration interim measure against a direct density integral.

Fixed geometry and the circular-duration/g approximation only; not an exact
eccentric transit model or an ALDERAAN population-recovery validation.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from duration_measure_control import duration_from_density, density_and_derivative


def check():
    rows=[]
    reference=quad(lambda rho: np.exp(-.5*((rho-1)/.2)**2), .1, 2.8,
                   epsabs=1e-12)[0]
    for b in [0., .5, .8]:
        for e in [0., .3, .8]:
            for omega in [0., np.pi/2, 3*np.pi/2]:
                low=duration_from_density(2.8,e,omega,b=b)
                high=duration_from_density(.1,e,omega,b=b)
                def integrand(t, corrected):
                    rho, derivative=density_and_derivative(t,e,omega,b=b)
                    density=np.exp(-.5*((rho-1)/.2)**2)
                    # Canceling the 1/T proposal and changing variables to rho.
                    return density*derivative if corrected else density/t
                corrected=quad(lambda t: integrand(t,True),low,high,epsabs=1e-12)[0]
                uncorrected=quad(lambda t: integrand(t,False),low,high,epsabs=1e-12)[0]
                rows.append(dict(b=b,e=e,omega=omega,corrected=corrected,
                                 uncorrected=uncorrected,reference=reference))
    return rows


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    rows=check()
    error=max(abs(r['corrected']/r['reference']-1) for r in rows)
    if error>1e-9:
        raise ValueError('Change-of-variable identity failed')
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(dict(maximum_relative_error=error,rows=rows,
        limitations=__doc__),indent=2))
    print(f'27 geometry cases: maximum relative integration error {error:.3g}')


if __name__=='__main__':
    main()
