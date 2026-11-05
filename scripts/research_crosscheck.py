#!/usr/bin/env python3
"""Execute the reviewed mcal analytic function on cube's published symmetric case.

Needs NumPy and a local checkout of the exact mcal commit recorded in sources.json.
Does not install software, run quantum calculations, or change the upstream source.
"""
import argparse
import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mcal-source', type=Path, default=ROOT/'build/research_sources/mcal')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    import numpy as np

    checkout = args.mcal_source.resolve()
    manifest = json.loads((ROOT/'docs/research/sources.json').read_text())
    source = next(x for x in manifest['repositories'] if x['name'] == 'mcal')
    commit = subprocess.check_output(['git', '-C', str(checkout), 'rev-parse', 'HEAD'], text=True).strip()
    if commit != source['commit']:
        raise ValueError('mcal commit differs from the reviewed version')
    relative = 'src/mcal/calculations/hopping_mobility_model.py'
    path = checkout/relative
    file_hash = hashlib.sha256(path.read_bytes()).hexdigest()
    expected = next(x['sha256'] for x in source['files'] if x['path'] == relative)
    if file_hash != expected:
        raise ValueError('mcal analytic source has changed')
    spec = importlib.util.spec_from_file_location('reviewed_mcal_model', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    case = ROOT/'runs/e_006_ether_12-17_0'
    audit = json.loads((case/'input_audit.json').read_text())
    reference = json.loads((case/'summary.json').read_text())
    with (case/'rates_and_geometry.csv').open(newline='') as stream:
        rows = list(csv.DictReader(stream))
    result = {
        'date': '2026-11-05',
        'system': reference['system'],
        'upstream_commit': commit,
        'upstream_function_sha256': file_hash,
        'python': sys.version.split()[0],
        'numpy': np.__version__,
        'scope': 'same imported rates; symmetric static network; no disorder or field; analytic tensor only',
        'not_tested': ['upstream quantum parameter calculation', 'complete upstream CLI',
                       'non-symmetric or finite-field network', 'full upstream test suite'],
        'input_sha256': {name: hashlib.sha256((case/name).read_bytes()).hexdigest()
                         for name in ['input_audit.json', 'summary.json', 'rates_and_geometry.csv']},
        'carriers': {},
    }
    for carrier in ['hole', 'electron']:
        hops = [(int(r['mol_i'])-1, int(r['mol_j'])-1,
                 int(r['image_a']), int(r['image_b']), int(r['image_c']),
                 float(r[f'{carrier}_rate_s_1'])) for r in rows]
        actual = module.diffusion_coefficient_tensor(np.diag(audit['box_nm']), hops)
        target = np.array(reference['carriers'][carrier]['ensemble']['exact_diffusion_tensor_nm2_s'])
        relative_error = float(np.linalg.norm(actual-target)/np.linalg.norm(target))
        if not np.isfinite(relative_error) or relative_error >= 1e-9:
            raise ValueError(f'{carrier}: tensor mismatch {relative_error}')
        result['carriers'][carrier] = {
            'relative_frobenius_error': relative_error,
            'mcal_diffusion_nm2_per_s': actual.tolist(),
            'cube_diffusion_nm2_per_s': target.tolist(),
            'passed': True,
        }
    result['passed'] = True
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: v['relative_frobenius_error'] for k, v in result['carriers'].items()}))


if __name__ == '__main__':
    main()
