#!/usr/bin/env python3
"""Run actual CPU/GPU QC comparisons, not just a CUDA availability probe.

Run in the repository's uv quantum-cuda13 environment. Exit nonzero on failed
acceptance. Output retains every result and the exact environment/code hashes.
Small-system timings are observations, not a general GPU speedup benchmark.
"""
import argparse
from dataclasses import replace
import hashlib
import importlib.metadata
import io
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time
import unittest

import numpy as np

from universal_kmc.quantum.transfer import SCFConfig, calculate_transfer, read_xyz, METHODS


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--case', action='append', choices=('cofacial', 'slipped', 'distorted', 'separated', 'hf', 'multiorbital', 'benzene', 'rigid_transform'))
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    args.output.mkdir(parents=True, exist_ok=True)
    a = read_xyz(root/'examples/quantum/ethylene_a.xyz')

    def moved(atoms, translation, stretch=1.):
        return [(symbol, (float(x*stretch+translation[0]), float(y+translation[1]), float(z+translation[2]))) for symbol, (x, y, z) in atoms]

    benzene = [(symbol, (float(radius*np.cos(angle)), float(radius*np.sin(angle)), 0.))
               for symbol, radius in [('C', 1.397), ('H', 2.479)]
               for angle in np.arange(6)*np.pi/3]

    def rigid(atoms):
        return [(symbol, (float(x+4), float(-z-2), float(y+1))) for symbol, (x, y, z) in atoms]

    cases = {
        'cofacial': (moved(a, [0, 0, 3.5]), SCFConfig(), None),
        'slipped': (moved(a, [1.5, 0, 3.5]), SCFConfig(), None),
        'distorted': (moved(a, [1.2, 0, 3.5], stretch=1.12), SCFConfig(), None),
        'separated': (moved(a, [0, 0, 8.0]), SCFConfig(), None),
        'hf': (moved(a, [0, 0, 3.5]), SCFConfig(xc='hf', basis='sto-3g'), None),
        'multiorbital': (moved(a, [0, 0, 3.5]), SCFConfig(), {'hole': ([6, 7], [6, 7]), 'electron': ([8, 9], [8, 9])}),
        'benzene': (moved(benzene, [.5, 0, 3.5]), SCFConfig(), None),
        'rigid_transform': (rigid(moved(a, [0, 0, 3.5])), SCFConfig(), None),
    }
    summary = {'schema': 'cube-transfer-validation-1', 'host_platform': platform.platform(),
               'python': sys.version, 'gpu': subprocess.check_output(['nvidia-smi', '--query-gpu=name,memory.total,driver_version,compute_cap', '--format=csv,noheader'], text=True).strip(),
               'tolerances': {'total_energy_hartree_atol': 1e-6, 'coupling_ev_atol': 2e-5, 'coupling_rtol': 2e-4},
               'cases': {}, 'all_passed': False,
               'timing_scope': 'synchronized calculation wall time; small QC validation, not a speedup benchmark'}
    hashes = {}
    for path in [root/'pyproject.toml', root/'uv.lock', Path(__file__), *sorted((root/'src/universal_kmc/quantum').glob('*.py')), root/'tests/test_quantum_coupling.py']:
        hashes[str(path.relative_to(root))] = hashlib.sha256(path.read_bytes()).hexdigest()
    summary['file_sha256'] = hashes
    import universal_kmc.quantum as quantum_package
    installed = Path(quantum_package.__file__).resolve().parent
    for source in sorted((root/'src/universal_kmc/quantum').glob('*.py')):
        assert (installed/source.name).read_bytes() == source.read_bytes(), ('installed/source mismatch', source.name)
    summary['installed_implementation_matches_source'] = True
    summary['installed_from_source_tree'] = installed == (root/'src/universal_kmc/quantum').resolve()
    direct_url = importlib.metadata.distribution('universal-kmc').read_text('direct_url.json')
    if direct_url:
        metadata = json.loads(direct_url)
        summary['installed_wheel_hashes'] = metadata.get('archive_info', {}).get('hashes', {})
    started = time.perf_counter()
    try:
        # GPU algebra tests are explicitly enabled, and skipped GPU tests fail this gate.
        os.environ['CUBE_GPU_TESTS'] = '1'
        suite = unittest.defaultTestLoader.discover(str(root/'tests'), pattern='test_quantum_coupling.py')
        log = io.StringIO()
        result = unittest.TextTestRunner(stream=log, verbosity=2).run(suite)
        (args.output/'algebra-tests.log').write_text(log.getvalue())
        summary['algebra_tests'] = {'run': result.testsRun, 'failures': len(result.failures), 'errors': len(result.errors), 'skipped': len(result.skipped)}
        assert result.wasSuccessful() and not result.skipped, log.getvalue()
        for name in args.case or cases:
            b, config, indices = cases[name]
            case_a = benzene if name == 'benzene' else rigid(a) if name == 'rigid_transform' else a
            runs = {}
            for backend in ('cpu', 'gpu'):
                print(f'{name}: {backend} SCF and four coupling methods', flush=True)
                runs[backend] = calculate_transfer(case_a, b, config=replace(config, backend=backend),
                                                   methods=METHODS, assume_equivalent=True,
                                                   orbital_indices=indices)
                (args.output/f'{name}-{backend}.json').write_text(json.dumps(runs[backend], indent=2, allow_nan=False)+'\n')
                assert all(x['converged'] for x in runs[backend]['scf'].values())
                if backend == 'gpu':
                    assert runs[backend]['device'] and all(x['class'].startswith('gpu4pyscf.') for x in runs[backend]['scf'].values())
            differences = []
            for fragment in ('fragment_a', 'fragment_b', 'dimer'):
                differences.append(abs(runs['cpu']['scf'][fragment]['energy_hartree']-runs['gpu']['scf'][fragment]['energy_hartree']))
            assert max(differences) <= summary['tolerances']['total_energy_hartree_atol'], (name, differences)
            compared = []
            rejected = {}
            for carrier in ('hole', 'electron'):
                cpu, gpu = runs['cpu']['results'][carrier], runs['gpu']['results'][carrier]
                assert cpu['orbital_indices_a'] == gpu['orbital_indices_a'] and cpu['orbital_indices_b'] == gpu['orbital_indices_b']
                for method in METHODS:
                    cr, gr = cpu['methods'][method], gpu['methods'][method]
                    assert cr['status'] == gr['status'], (name, carrier, method, cr, gr)
                    if cr['status'] != 'ok':
                        rejected[carrier+'/'+method] = cr['reason']
                        continue
                    # Singular values avoid dependence on arbitrary signs/rotations of
                    # degenerate fragment orbitals. Single-state magnitudes also checked.
                    if 'singular_values_ev' in cr:
                        cv, gv = np.array(cr['singular_values_ev']), np.array(gr['singular_values_ev'])
                    else:
                        cv, gv = np.array([cr['scalar_magnitude_ev']]), np.array([gr['scalar_magnitude_ev']])
                    np.testing.assert_allclose(cv, gv, atol=2e-5, rtol=2e-4, err_msg=f'{name}/{carrier}/{method}')
                    compared.append({'carrier': carrier, 'method': method, 'maximum_difference_ev': float(np.max(abs(cv-gv)))})
                for method in ('dipro', 'frozen_density'):
                    assert cpu['methods'][method]['status'] == 'ok'
                if name in ('cofacial', 'slipped', 'hf', 'rigid_transform'):
                    assert cpu['methods']['esid']['status'] == 'ok', (name, carrier, cpu['methods']['esid'])
                    assert cpu['methods']['orbital_gmh']['status'] == 'ok'
                if name in ('distorted', 'multiorbital', 'benzene'):
                    assert cpu['methods']['esid']['status'] == 'not_applicable'
                if name in ('multiorbital', 'benzene'):
                    assert np.array(cpu['methods']['dipro']['coupling_ev']).shape == (2, 2)
                    assert cpu['methods']['dipro']['scalar_magnitude_ev'] is None
            summary['cases'][name] = {'passed': True, 'maximum_scf_energy_difference_hartree': max(differences),
                                       'comparisons': compared, 'expected_method_rejections': rejected}
            print(f'{name}: PASS; maximum SCF difference {max(differences):.3g} Hartree', flush=True)
        if 'cofacial' in summary['cases'] and 'separated' in summary['cases']:
            for backend in ('cpu', 'gpu'):
                close = json.loads((args.output/f'cofacial-{backend}.json').read_text())
                far = json.loads((args.output/f'separated-{backend}.json').read_text())
                for carrier in ('hole', 'electron'):
                    c = close['results'][carrier]['methods']['dipro']['scalar_magnitude_ev']
                    f = far['results'][carrier]['methods']['dipro']['scalar_magnitude_ev']
                    assert f < c/10, (backend, carrier, c, f)
            summary['distance_decay_check'] = True
        if 'cofacial' in summary['cases'] and 'rigid_transform' in summary['cases']:
            for backend in ('cpu', 'gpu'):
                close = json.loads((args.output/f'cofacial-{backend}.json').read_text())
                transformed = json.loads((args.output/f'rigid_transform-{backend}.json').read_text())
                for carrier in ('hole', 'electron'):
                    for method in METHODS:
                        c = close['results'][carrier]['methods'][method]['scalar_magnitude_ev']
                        r = transformed['results'][carrier]['methods'][method]['scalar_magnitude_ev']
                        np.testing.assert_allclose(c, r, atol=2e-5, rtol=2e-4)
            summary['quantum_rigid_transform_invariance'] = True
        print('GPU frozen-density-only path (no dimer SCF)', flush=True)
        frozen = calculate_transfer(a, cases['cofacial'][0], config=SCFConfig(backend='gpu'), methods=('frozen_density',))
        assert 'dimer' not in frozen['scf']
        (args.output/'frozen-only-gpu.json').write_text(json.dumps(frozen, indent=2, allow_nan=False)+'\n')
        summary['frozen_density_without_dimer_scf'] = True
        failures = {}
        for backend in ('cpu', 'gpu'):
            try:
                calculate_transfer(a, cases['distorted'][0], config=SCFConfig(backend=backend, max_cycle=1))
            except RuntimeError as exc:
                assert 'did not converge' in str(exc)
                failures[backend] = str(exc)
            else:
                raise AssertionError('unconverged SCF was incorrectly accepted')
        summary['unconverged_scf_rejected'] = failures
        command = [sys.executable, '-m', 'universal_kmc.quantum.cli', '--fragment-a', str(root/'examples/quantum/ethylene_a.xyz'),
                   '--fragment-b', str(root/'examples/quantum/ethylene_b.xyz'), '--backend', 'gpu', '--basis', 'sto-3g', '--xc', 'hf']
        for method, expected in (('dipro', 0), ('esid', 2)):
            path = args.output/f'cli-{method}.json'
            proc = subprocess.run(command+['--method', method, '--output', str(path)], text=True, capture_output=True)
            assert proc.returncode == expected, (method, proc.returncode, proc.stderr)
            assert path.exists()
            summary.setdefault('cli', {})[method] = {'exit_code': proc.returncode, 'stderr': proc.stderr}
        failed_path = args.output/'cli-failed.json'
        failed_path.write_text((args.output/'cli-dipro.json').read_text())
        proc = subprocess.run(command+['--max-cycle', '1', '--output', str(failed_path)], text=True, capture_output=True)
        failed = json.loads(failed_path.read_text())
        assert proc.returncode == 1 and failed['status'] == 'failed' and 'results' not in failed
        summary['cli']['failed_run_replaces_stale_success'] = {'exit_code': proc.returncode, 'status': failed['status']}
        summary['all_passed'] = True
    except Exception as exc:
        summary['failure'] = {'type': type(exc).__name__, 'message': str(exc)}
        raise
    finally:
        summary['wall_time_seconds'] = time.perf_counter()-started
        (args.output/'summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n')
    print('All CPU/GPU coupling validation gates passed.', flush=True)


if __name__ == '__main__':
    main()
