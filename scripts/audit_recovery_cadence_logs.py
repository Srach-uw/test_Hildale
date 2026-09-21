"""Compare recovery file inventories with explicit processed-cadence logs."""
import argparse
import hashlib
import json
from pathlib import Path
import pandas as pd


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    rows=[]
    for provenance in sorted((a.root/'provenance').rglob('*.tsv')):
        fields={}
        for line in provenance.read_text(errors='replace').splitlines():
            if '\t' in line:
                k,v=line.split('\t',1); fields[k]=v
        if 'target' not in fields or 'short_cadence_files' not in fields:
            continue
        archive=provenance.parent.name
        log=a.root/'logs'/archive/'targets'/fields['run_id']/(fields['target']+'.stdout.log')
        text=log.read_text(errors='replace') if log.exists() else ''
        n=int(fields['short_cadence_files'])
        rows.append(dict(target=fields['target'],archive=archive,available_sc_files=n,
            log_available=log.exists(),explicit_no_sc='No short cadence data' in text,
            inventory_vs_log_flag=n>0 and 'No short cadence data' in text,
            provenance_sha256=hashlib.sha256(provenance.read_bytes()).hexdigest(),
            log_sha256=hashlib.sha256(log.read_bytes()).hexdigest() if log.exists() else '',
            provenance_path=str(provenance),log_path=str(log)))
    frame=pd.DataFrame(rows)
    if frame.empty:
        raise ValueError('No provenance records found')
    a.output.mkdir(parents=True,exist_ok=True)
    frame.to_csv(a.output/'cadence_inventory_vs_logs.csv',index=False)
    scope=dict(records=len(frame),unique_targets=frame.target.nunique(),
        available_sc_records=int(frame.available_sc_files.gt(0).sum()),
        flagged_records=int(frame.inventory_vs_log_flag.sum()),
        missing_logs=int((~frame.log_available).sum()),
        caveat='File inventory is not proof of usable cadence; explicit log messages can follow filtering. Flags require investigating loader and quarter selection, not automatic refits.')
    (a.output/'scope.json').write_text(json.dumps(scope,indent=2))
    print(json.dumps(scope,indent=2))
    print(frame.loc[frame.inventory_vs_log_flag,['target','archive','available_sc_files']].to_string(index=False))


if __name__=='__main__':
    main()
