#!/usr/bin/env python3
"""Actual PySCF DOS validation; --gpu requires a working GPU and never falls back."""
import argparse
import json
from pathlib import Path
import platform
import time

import numpy as np
import pyscf

from wren.dos.control import normalize
from wren.dos.core import broaden, bloch_spectra, kmesh
from wren.dos.runner import run
from wren.quantum.cli import _write_result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('build/dos-validation'))
    parser.add_argument('--gpu', action='store_true')
    parser.add_argument('--md-xyz', type=Path, help='optional actual monomer geometry, separate from H2 checks')
    args = parser.parse_args()
    root = args.output.resolve()
    if root.exists() and any(root.iterdir()):
        parser.error('validation output must be new/empty')
    root.mkdir(parents=True, exist_ok=True)
    (root/'h2.xyz').write_text('2\nH2 DOS functional validation; Angstrom\nH 0 0 0\nH .74 0 0\n')
    report = {'schema': 'wren-dos-device-validation-1', 'status': 'running',
              'versions': {'python': platform.python_version(), 'numpy': np.__version__, 'pyscf': pyscf.__version__},
              'cases': {}, 'scope': 'small-system implementation validation, not material convergence'}
    output = root/'validation.json'
    _write_result(output, report)
    try:
        settings = {'xyz': 'h2.xyz', 'method': 'hf', 'basis-set': 'sto-3g', 'threads': 4,
                    'sigma': .1, 'output': 'h2-cpu'}
        cpu = run(normalize(settings, root))
        report['cases']['h2-cpu'] = {k: cpu[k] for k in ('scope', 'elapsed_seconds', 'expected_integral', 'finite_grid_integral')}
        if args.gpu:
            import cupy as cp
            import gpu4pyscf
            device = cp.cuda.runtime.getDeviceProperties(0)
            name = device['name']
            report.update(device=name.decode() if isinstance(name, bytes) else str(name))
            report['versions'].update(cupy=cp.__version__, gpu4pyscf=gpu4pyscf.__version__)
            gpu = run(normalize({**settings, 'qc-engine': 'gpu-pyscf', 'dos-backend': 'gpu', 'output': 'h2-gpu'}, root))
            raw = [json.loads((root/f'h2-{kind}/spectra.json').read_text())['frames'][0]['molecules'][0]
                   for kind in ('cpu', 'gpu')]
            error = float(np.max(abs(np.array(raw[0]['orbital_energies_ev'])-raw[1]['orbital_energies_ev'])))
            if error > 5e-5:
                raise RuntimeError(f'CPU/GPU H2 orbital disagreement: {error} eV')
            report['cases']['h2-gpu'] = {**{k: gpu[k] for k in ('elapsed_seconds', 'expected_integral', 'finite_grid_integral')},
                'max_orbital_error_ev': error, 'actual_scf_class': raw[1]['scf']['source_class']}
            # Warm both implementations, then time the same spectrum and grid.
            rng = np.random.default_rng(20261123)
            levels = rng.normal(0., .5, 1000)
            grid = np.linspace(-3., 3., 4001)
            curves, timings = {}, {}
            for backend in ('cpu', 'gpu'):
                broaden(levels[:4], grid[:16], .05, backend=backend)
                if backend == 'gpu':
                    cp.cuda.Stream.null.synchronize()
                start = time.perf_counter()
                curves[backend] = broaden(levels, grid, .05, backend=backend)
                if backend == 'gpu':
                    cp.cuda.Stream.null.synchronize()
                timings[backend] = time.perf_counter()-start
            delta = float(abs(curves['cpu']-curves['gpu']).max())
            if not np.allclose(curves['cpu'], curves['gpu'], rtol=1e-10, atol=1e-10):
                raise RuntimeError('CPU/GPU Gaussian DOS mismatch')
            points, _ = kmesh([8, 1, 1])
            pairs = [(i, i+1, [0, 0, 0], float(rng.normal(0., .1))) for i in range(15)]
            pairs.append((15, 0, [1, 0, 0], .1))
            sites = rng.normal(0., .2, 16)
            ecpu = bloch_spectra(sites, pairs, points)
            egpu = bloch_spectra(sites, pairs, points, backend='gpu')
            if not np.allclose(ecpu, egpu, rtol=1e-10, atol=1e-10):
                raise RuntimeError('CPU/GPU Bloch spectra mismatch')
            report['cases']['gaussian'] = {'max_abs_error': delta, 'wall_seconds_including_transfer': timings,
                'levels': len(levels), 'grid_points': len(grid), 'gpu_warmed': True}
            report['cases']['bloch'] = {'sites': len(sites), 'kpoints': len(points),
                                        'max_abs_error_ev': float(abs(ecpu-egpu).max())}
        if args.md_xyz:
            backend = 'gpu-pyscf' if args.gpu else 'pyscf'
            result = run(normalize({'xyz': str(args.md_xyz.resolve()), 'method': 'hf', 'basis-set': 'sto-3g',
                'auxbasis': 'weigend', 'conv-tol': 1e-7, 'threads': 8, 'qc-engine': backend,
                'dos-backend': 'gpu' if args.gpu else 'cpu', 'output': 'md-monomer'}, root))
            report['cases']['md-monomer'] = {k: result[k] for k in
                ('scope', 'elapsed_seconds', 'expected_integral', 'finite_grid_integral', 'frontier_disorder')}
            report['cases']['md-monomer']['note'] = 'HF/STO-3G with Weigend fitting: functionality only; basis/auxiliary basis not material-converged'
        report['status'] = 'completed'
    except Exception as exc:
        report.update(status='failed', error={'type': type(exc).__name__, 'message': str(exc)})
        _write_result(output, report)
        raise
    _write_result(output, report)
    print(f'wren: DOS validation {report["status"]}: {output}')


if __name__ == '__main__':
    main()
