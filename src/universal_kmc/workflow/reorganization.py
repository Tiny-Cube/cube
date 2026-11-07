"""Vacuum or frozen-MM electrostatic four-point reorganization energies.

lambda = [Eq(R0)-Eq(Rq)] + [E0(Rq)-E0(R0)]. Both R0 and Rq are
optimized; MD input coordinates are only starting geometries. Frozen MM
coordinates/charges are identical at all four points. No polarizable MM,
QM-MM van der Waals, solvent or MM outer-sphere relaxation is included.
Optimization API: https://pyscf.org/user/geomopt.html
"""
from dataclasses import asdict
import json
from pathlib import Path
import time

import numpy as np

from universal_kmc.quantum.cli import _write_result
from universal_kmc.quantum.coupling import host
from universal_kmc.quantum.transfer import _versions
from .control import scf_config
from .dimers import file_hash, fingerprint, xyz_text
from .scf import make_mean_field, run_single_point, scf_record


def four_point_terms(e0_r0, eq_r0, eq_rq, e0_rq):
    """Return signed relaxation terms in caller units; never take abs()."""
    values = np.asarray([e0_r0, eq_r0, eq_rq, e0_rq], dtype=float)
    if not np.isfinite(values).all():
        raise ValueError('four-point energies must be finite')
    charged, neutral = float(eq_r0-eq_rq), float(e0_rq-e0_r0)
    return {'charged_relaxation': charged, 'neutral_relaxation': neutral, 'total': charged+neutral}


def _optimize(atoms, config, charge, directory, maxsteps, *, resume=False, embedding=None):
    from pyscf.geomopt import as_pyscf_method, geometric_solver
    spin = 0 if charge == 0 else 1
    label = 'neutral' if charge == 0 else ('cation' if charge > 0 else 'anion')
    from . import scf, embedding as embedding_module
    key = fingerprint({'atoms': atoms, 'config': asdict(config), 'charge': charge, 'maxsteps': maxsteps,
                       'embedding_fingerprint': embedding.fingerprint if embedding else None,
                       'implementation_sha256': {p: file_hash(p) for p in (__file__, scf.__file__, embedding_module.__file__)}})
    metadata = directory/f'{label}.optimization.json'
    if resume and metadata.exists():
        old = json.loads(metadata.read_text())
        if old.get('fingerprint') != key:
            raise ValueError('optimization inputs changed; choose a different output directory')
        if old.get('status') == 'completed':
            for relative, expected in old['output_sha256'].items():
                if not (directory/relative).is_file() or file_hash(directory/relative) != expected:
                    raise ValueError('cannot resume changed/missing optimization output')
            return old['atoms'], old['optimization']
    _write_result(metadata, {'fingerprint': key, 'status': 'running'})
    mf = None
    trace = []

    def evaluate(mol):
        if embedding is not None:
            embedding.validate_qm([(mol.atom_symbol(i), r) for i, r in enumerate(mol.atom_coords(unit='Angstrom'))])
            embedding.assert_attached(scanner.base)
        energy, gradient = scanner(mol)
        record = scf_record(scanner.base, config.backend)
        if embedding is not None:
            embedding.assert_attached(scanner.base)
            record.update(embedding_fingerprint=embedding.fingerprint, mm_frozen=True)
        gradient = np.asarray(host(gradient), dtype=float)
        if gradient.shape != (mol.natm, 3) or not np.isfinite(gradient).all():
            raise RuntimeError('optimizer received an invalid energy gradient')
        trace.append({'evaluation': len(trace)+1, 'energy_hartree': float(energy),
                      'coordinates_angstrom': mol.atom_coords(unit='Angstrom').tolist(),
                      'gradient_hartree_per_bohr': gradient.tolist(),
                      'gradient_rms': float(np.sqrt(np.mean(gradient**2))),
                      'gradient_max': float(abs(gradient).max()), 'scf': record})
        _write_result(directory/f'{label}.trace.json', trace)
        return float(energy), gradient

    try:
        mf = make_mean_field(atoms, config, charge=charge, spin=spin, log=directory/f'{label}.optimization.log', embedding=embedding)
        scanner = mf.nuc_grad_method().as_scanner()
        gradient_class = type(scanner).__module__+'.'+type(scanner).__name__
        if config.backend == 'gpu':
            from gpu4pyscf.grad.rhf import GradientsBase
            if not isinstance(scanner, GradientsBase):
                raise RuntimeError('GPU optimization did not produce GPU4PySCF gradients')
        if embedding is not None:
            if config.backend == 'gpu':
                from gpu4pyscf.qmmm.itrf import QMMMGrad
            else:
                from pyscf.qmmm.itrf import QMMMGrad
            if not isinstance(scanner, QMMMGrad):
                raise RuntimeError('QM/MM optimizer is missing MM contributions to the QM gradient')
        # The optimizer is on CPU; evaluate() explicitly executes the requested
        # CPU/GPU energy+gradient and transfers only the optimizer's gradient.
        adapter = as_pyscf_method(mf.mol, evaluate)
        converged, optimized = geometric_solver.kernel(adapter, maxsteps=maxsteps, assert_convergence=True)
        if not converged or not trace:
            raise RuntimeError(f'{label} geometry optimization did not converge in {maxsteps} steps')
        xyz = optimized.atom_coords(unit='Angstrom')
        result = [(optimized.atom_symbol(i), tuple(map(float, xyz[i]))) for i in range(optimized.natm)]
        environment_name = 'vacuum' if embedding is None else 'frozen-MM electrostatic QM/MM'
        (directory/f'{label}.optimized.xyz').write_text(xyz_text([a[0] for a in result], xyz, f'{label} optimized in {environment_name}'))
        details = {'converged': True, 'evaluations': len(trace), 'gradient_class': gradient_class,
                         'final_gradient_rms': trace[-1]['gradient_rms'],
                         'final_gradient_max': trace[-1]['gradient_max'], 'charge': charge, 'spin_2s': spin,
                         'gradient_mro': [c.__module__+'.'+c.__name__ for c in type(scanner).__mro__]}
        if embedding is not None:
            details.update(mm_frozen=True, embedding_fingerprint=embedding.fingerprint,
                           optimized_region='QM_atoms_only', qm_mm_gradient_verified=True)
        _write_result(metadata, {'fingerprint': key, 'status': 'completed', 'atoms': result, 'optimization': details,
            'output_sha256': {f'{label}{suffix}': file_hash(directory/f'{label}{suffix}')
                              for suffix in ('.trace.json', '.optimized.xyz')}})
        return result, details
    except Exception as exc:
        _write_result(metadata, {'fingerprint': key, 'status': 'failed', 'error': str(exc)})
        raise
    finally:
        if mf is not None:
            mf.mol.stdout.close()


def calculate_reorganization(atoms, settings, directory, *, resume=False, embedding=None):
    from pyscf.data.nist import HARTREE2EV
    if settings['environment'] not in ('vacuum', 'qmmm') or settings['re-method'] != '4-points':
        raise ValueError('use vacuum or frozen-MM qmmm four-point reorganization')
    if (settings['environment'] == 'qmmm') != (embedding is not None):
        raise ValueError('QM/MM requires a frozen MM environment; vacuum must not have MM charges')
    if embedding is not None:
        if not settings.get('mm-frozen', True):
            raise ValueError('only frozen MM coordinates and charges are supported')
        embedding.validate_qm(atoms)
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    config = scf_config(settings)
    from . import scf, embedding as embedding_module
    request = {'atoms': atoms, 'settings': asdict(config), 'carriers': settings['carriers'],
               'environment': settings['environment'], 'opt_maxsteps': settings['opt-maxsteps'],
               'versions': _versions(), 'implementation_sha256': file_hash(__file__)}
    request['support_implementation_sha256'] = {p: file_hash(p) for p in (scf.__file__, embedding_module.__file__)}
    if embedding is not None:
        request.update(embedding_fingerprint=embedding.fingerprint, mm_frozen=True)
    key = fingerprint(request)
    output = directory/'reorganization.json'
    if output.exists():
        old = json.loads(output.read_text())
        if not resume or old.get('fingerprint') != key:
            raise ValueError(f'existing reorganization task {directory}; use --resume with identical inputs')
        if old['status'] == 'completed':
            for relative, expected in old['output_sha256'].items():
                if not (directory/relative).is_file() or file_hash(directory/relative) != expected:
                    raise ValueError(f'cannot resume: changed/missing reorganization output {relative}')
            return {**old, 'reused': True}
    record = {'schema': 'cube-reorganization-1', 'status': 'running', 'fingerprint': key,
              'input': request, 'environment': settings['environment'], 'method': '4-points',
              'scope': ('intramolecular_vacuum_no_MM_no_outer_sphere' if embedding is None else
                        'QM_relaxation_in_frozen_point_charge_MM_no_QM_MM_vdw_no_outer_sphere_relaxation'), 'results': {}}
    if embedding is not None:
        from .embedding import snapshot
        snapshot(embedding, directory/'mm-environment.json')
        record['embedding'] = {k: v for k, v in embedding.payload().items()
                               if k not in ('coordinates_angstrom', 'charges_e', 'provenance')}
        record['embedding'].update(fingerprint=embedding.fingerprint, snapshot='mm-environment.json')
    _write_result(output, record)
    started = time.perf_counter()
    try:
        r0, opt0 = _optimize(atoms, config, 0, directory, settings['opt-maxsteps'], resume=resume, embedding=embedding)
        e0r0 = run_single_point(r0, settings, directory/'E0_R0', resume=resume, embedding=embedding)
        record['neutral_optimization'] = opt0
        for carrier in settings['carriers']:
            charge = 1 if carrier == 'hole' else -1
            rq, optq = _optimize(r0, config, charge, directory, settings['opt-maxsteps'], resume=resume, embedding=embedding)
            eqr0 = run_single_point(r0, settings, directory/f'{carrier}_Eq_R0', resume=resume, charge=charge, spin=1, embedding=embedding)
            eqrq = run_single_point(rq, settings, directory/f'{carrier}_Eq_Rq', resume=resume, charge=charge, spin=1, embedding=embedding)
            e0rq = run_single_point(rq, settings, directory/f'{carrier}_E0_Rq', resume=resume, embedding=embedding)
            energies = {'E0_R0': e0r0['energy_hartree'], 'Eq_R0': eqr0['energy_hartree'],
                        'Eq_Rq': eqrq['energy_hartree'], 'E0_Rq': e0rq['energy_hartree']}
            terms = four_point_terms(*(energies[k] for k in ('E0_R0','Eq_R0','Eq_Rq','E0_Rq')))
            if min(terms['charged_relaxation'], terms['neutral_relaxation'])*HARTREE2EV < -1e-5:
                raise RuntimeError('negative relaxation energy exceeds tolerance; inspect optimized states and SCF solutions')
            record['results'][carrier] = {'energies_hartree': energies, 'relaxation_hartree': terms,
                'relaxation_ev': {k: v*HARTREE2EV for k,v in terms.items()},
                'reorganization_energy_ev': terms['total']*HARTREE2EV, 'charged_optimization': optq,
                'single_points': {'E0_R0': e0r0, 'Eq_R0': eqr0, 'Eq_Rq': eqrq, 'E0_Rq': e0rq}}
            _write_result(output, record)
        record.update(status='completed', elapsed_seconds=time.perf_counter()-started, reused=False,
                      output_sha256={str(p.relative_to(directory)): file_hash(p) for p in directory.iterdir()
                                     if p.is_file() and p != output})
        _write_result(output, record)
        return record
    except Exception as exc:
        record.update(status='failed', error={'type': type(exc).__name__, 'message': str(exc)})
        _write_result(output, record)
        raise
