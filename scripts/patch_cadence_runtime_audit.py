"""Instrument a staged ALDERAAN driver without changing its likelihood."""
import argparse
import hashlib
import json
from pathlib import Path

ANCHOR = '    texp[all_dtype == "long"] = lcit\n'
MARKER = '# Recording selected cadence inputs.'
BLOCK = '''
    # Recording selected cadence inputs.
    import json as _cadence_json
    _cadence_records = []
    for _q in quarters:
        _cadence_records.append(dict(
            quarter=int(_q), cadence=str(all_dtype[_q]),
            selected_points=int(len(all_time[_q])),
            per_planet_mask_points=np.sum(all_mask[_q], axis=1).astype(int).tolist(),
            exposure_days=float(texp[_q]), oversample=int(oversample[_q])))
    _cadence_path = os.path.join(RESULTS_DIR, f"{TARGET}_sampler_cadence.json")
    with open(_cadence_path, "w") as _cadence_file:
        _cadence_json.dump(dict(target=TARGET, run_id=RUN_ID,
            quarters=_cadence_records), _cadence_file, indent=2)
    if os.environ.get("ALDERAAN_REQUIRE_SC") == "1":
        if not any(r["cadence"] == "short" and r["selected_points"] > 0
                   for r in _cadence_records):
            raise RuntimeError("SC validation requires selected short-cadence transit points")
'''


def patched(source):
    if MARKER in source:
        raise ValueError('Runtime audit already present; refusing an ambiguous repatch')
    if source.count(ANCHOR) != 1:
        raise ValueError('Expected exposure assignment not uniquely present')
    result=source.replace(ANCHOR,ANCHOR+BLOCK)
    compile(result,'staged_alderaan_driver','exec')
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if a.output.exists() or a.source.resolve()==a.output.resolve():
        raise ValueError('Output must be a new staged file')
    before=a.source.read_bytes()
    result=patched(before.decode('utf-8').replace('\r\n','\n'))
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(result,encoding='utf-8')
    manifest=dict(source=str(a.source),source_sha256=hashlib.sha256(before).hexdigest(),
        output=str(a.output),output_sha256=hashlib.sha256(a.output.read_bytes()).hexdigest(),
        purpose='Runtime selected cadence record; optional fail-closed SC validation; no likelihood change')
    a.output.with_suffix('.audit.json').write_text(json.dumps(manifest,indent=2))
    print(json.dumps(manifest,indent=2))


if __name__=='__main__':
    main()
