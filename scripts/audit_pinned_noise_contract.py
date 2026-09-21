"""Inspect the pinned public likelihood without claiming runtime provenance."""

import argparse
import ast
import hashlib
import json
from pathlib import Path
from urllib.request import urlopen


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    commit = '7443dff16b7f9092e14a6f0cc1f8948d457c9e0b'
    root = f'https://raw.githubusercontent.com/gjgilbert/alderaan/{commit}/'
    evidence = {'commit': commit, 'files': {}, 'runtime_attestation': False}
    for name in ['alderaan/dynesty_helpers.py', 'bin/fit_transit_shape_simultaneous_nested.py']:
        with urlopen(root+name, timeout=45) as response:
            raw = response.read()
        source = raw.decode()
        tree = ast.parse(source)
        switches = []
        for node in ast.walk(tree):
            if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id == 'USE_GP' for t in node.targets):
                switches.append({'line': node.lineno, 'value': ast.literal_eval(node.value)})
        evidence['files'][name] = {'url': root+name, 'sha256': hashlib.sha256(raw).hexdigest(), 'gp_assignments': switches}
        args.output.mkdir(parents=True, exist_ok=True)
        (args.output/Path(name).name).write_bytes(raw)
    patches = Path(__file__).resolve().parents[1]/'cloud/ld_validation/patch_alderaan_repro.py'
    text = patches.read_text()
    evidence['current_repro_patch'] = {'sha256': hashlib.sha256(patches.read_bytes()).hexdigest(),
        'mentions_gp_switch': 'USE_GP' in text,
        'mentions_likelihood_helper': 'dynesty_helpers' in text}
    evidence['scope'] = 'Pinned upstream code plus current local reproducibility patch; not an execution-time source hash.'
    (args.output/'contract.json').write_text(json.dumps(evidence,indent=2)+'\n')
    print(json.dumps(evidence,indent=2))


if __name__ == '__main__':
    main()
