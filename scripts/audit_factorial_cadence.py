"""Record actual short-cadence evidence in the archived factorial experiment."""
import argparse
import hashlib
import json
from pathlib import Path


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    runner=a.root/'run_one_target.sh'
    command=next(line for line in runner.read_text().splitlines()
                 if line.startswith('python bin/detrend'))
    rows=[]
    for arm in ['original_lcsc','reference_lcsc']:
        for log in sorted((a.root/'logs'/'targets'/('sagear_validation_'+arm)).glob('*.stdout.log')):
            target=log.name.removesuffix('.stdout.log')
            products=list((a.root/'projects'/arm).rglob(target+'_sc_detrended.fits'))
            rows.append(dict(arm=arm,target=target,
                explicit_no_sc='No short cadence data' in log.read_text(errors='replace'),
                archived_sc_products=[str(x) for x in products],
                log_sha256=hashlib.sha256(log.read_bytes()).hexdigest()))
    result=dict(archived_runner_sha256=hashlib.sha256(runner.read_bytes()).hexdigest(),
        detrend_command=command,records=rows,
        limitation='Absence of an archived product alone does not prove absence at runtime; explicit logs and archived command provide corroboration.')
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,indent=2))
    print(f'{len(rows)} logs; {sum(r["explicit_no_sc"] for r in rows)} explicitly report no SC')


if __name__=='__main__':
    main()
