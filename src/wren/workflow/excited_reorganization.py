"""Restricted singlet SCF/TDA four-point energies in vacuum.

Use the lowest singlet TDA root throughout. Neighboring-root gaps and the
cross-geometry overlaps of transition amplitudes detect ambiguous steps.
These diagnostics neither establish diabatic identity nor verify minima.
PySCF conventions: https://pyscf.org/user/tddft.html
"""
from dataclasses import asdict
from importlib.metadata import version
import json
import math
from pathlib import Path
import time

import numpy as np

from wren import reorganization as core
from wren.quantum.cli import _write_result
from wren.quantum.geometry import normalize_fragment
from . import scf
from .control import scf_config
from .dimers import file_hash, fingerprint, xyz_text


def _screen_roots(energies_ev, overlaps, minimum_gap, minimum_overlap):
    if not np.isfinite(energies_ev).all() or np.min(np.diff(energies_ev), initial=0) < 0:
        raise RuntimeError('TDA returned non-finite or unordered roots')
    gap = min(energies_ev[0], energies_ev[1] - energies_ev[0]) if len(energies_ev) > 1 else energies_ev[0]
    if gap <= minimum_gap:
        raise RuntimeError(f'S1 is unstable or nearly degenerate: nearest energy gap {gap:.6g} eV')
    if overlaps is not None:
        if not np.isfinite(overlaps).all() or overlaps[0] < minimum_overlap:
            raise RuntimeError(f'low S1 transition-amplitude overlap: {overlaps[0]:.6g}')
        if len(overlaps) > 1 and overlaps[0] <= max(overlaps[1:]) + 1e-8:
            raise RuntimeError(f'ambiguous S1 root character, candidate overlaps {overlaps}')
    return float(gap)


class _SingletRoot:
    def __init__(self, minimum_gap, minimum_overlap):
        self.minimum_gap = minimum_gap
        self.minimum_overlap = minimum_overlap
        self.previous = None

    def evaluate(self, mf):
        from pyscf import gto
        from pyscf.data.nist import HARTREE2EV

        occ = mf.mo_coeff[:, mf.mo_occ == 2]
        vir = mf.mo_coeff[:, mf.mo_occ == 0]
        count = min(3, occ.shape[1] * vir.shape[1])
        if not count or occ.shape[1] + vir.shape[1] != mf.mo_coeff.shape[1]:
            raise ValueError('singlet TDA requires closed-shell occupied and virtual orbitals')
        td = mf.TDA()
        td.singlet = True
        td.max_cycle = mf.max_cycle
        # Include nonpositive TDA modes for instability detection. The PySCF
        # default discards low eigenvalues before returning the root list.
        td.positive_eig_threshold = -np.inf
        td.conv_tol = min(mf.conv_tol ** .5, 1e-7)
        excitation, xy = td.kernel(nstates=count)
        if len(excitation) != count or len(xy) != count or not np.all(td.converged):
            raise RuntimeError('singlet TDA roots did not converge')
        amplitudes = []
        for x, _ in xy:
            norm = np.linalg.norm(x)
            if not np.isfinite(norm) or norm <= 0:
                raise RuntimeError('invalid singlet TDA transition amplitude')
            amplitudes.append(x / norm)
        overlaps = None
        if self.previous is not None:
            prev_mol, prev_occ, prev_vir, prev_x = self.previous
            cross = gto.intor_cross('int1e_ovlp', prev_mol, mf.mol)
            occ_overlap = prev_occ.T @ cross @ occ
            vir_overlap = prev_vir.T @ cross @ vir
            # Frobenius product of the normalized occupied-to-virtual TDA
            # amplitudes after cross-AO projection. Absolute value removes
            # the arbitrary sign of each root; no orbital-index matching.
            projected = occ_overlap.T @ prev_x @ vir_overlap
            overlaps = [float(abs(np.vdot(projected, x))) for x in amplitudes]
        roots_ev = excitation * HARTREE2EV
        gap = _screen_roots(roots_ev, overlaps, self.minimum_gap, self.minimum_overlap)
        self.previous = (mf.mol.copy(), occ, vir, amplitudes[0])
        return td, {'root': 1, 'excitation_energies_ev': roots_ev.tolist(),
                    'nearest_gap_ev': gap, 'transition_amplitude_overlaps': overlaps}


def _optimize(mol, evaluator, directory, label, maxsteps):
    from pyscf.geomopt import as_pyscf_method, geometric_solver

    trace = []

    def energy_gradient(geom):
        e0, estate, grad, details = evaluator(geom, gradient=True)
        if grad.shape != (geom.natm, 3) or not np.isfinite(grad).all() or not math.isfinite(estate):
            raise RuntimeError(f'{label} optimization received invalid energy or gradient')
        trace.append({'ground_energy_hartree': e0, 'state_energy_hartree': estate,
                      'gradient_rms_hartree_per_bohr': float(np.sqrt(np.mean(grad**2))),
                      'coords_angstrom': geom.atom_coords(unit='Angstrom').tolist(), **details})
        _write_result(directory / f'{label}.trace.json', trace)
        return estate, grad

    method = as_pyscf_method(mol, energy_gradient)
    converged, geom = geometric_solver.kernel(method, maxsteps=maxsteps, assert_convergence=True)
    if not converged or not trace:
        raise RuntimeError(f'{label} geometry optimization failed to converge in {maxsteps} steps')
    symbols = [geom.atom_symbol(i) for i in range(geom.natm)]
    (directory / f'{label}.xyz').write_text(xyz_text(symbols, geom.atom_coords(unit='Angstrom'),
                                                    f'{label} optimized in vacuum'))
    return geom


def calculate_excited_reorganization(atoms, settings, directory, *, resume=False, state_a=None, state_b=None):
    """Calculate four ground/TDA-singlet total energies. XYZ coordinates are in Å."""
    from pyscf.data.nist import HARTREE2EV

    state_a = state_a or core.ElectronicState('S0', 0, 1, 0, 'scf')
    state_b = state_b or core.ElectronicState('S1', 0, 1, 1, 'tda')
    if (state_a.solver, state_a.multiplicity, state_a.root) != ('scf', 1, 0) or (
            state_b.solver, state_b.multiplicity, state_b.root) != ('tda', 1, 1):
        raise ValueError('automatic singlet solver requires SCF root 0 and singlet TDA root 1')
    if state_a.charge != state_b.charge or state_a.name == state_b.name:
        raise ValueError('SCF and TDA states require matching QM charge and distinct names')
    if settings.get('qc-engine') != 'pyscf' or settings.get('environment') != 'vacuum':
        raise ValueError('automatic singlet TDA four-point calculation requires CPU PySCF and vacuum')
    if settings.get('re-method') != '4-points' or settings['auxbasis'].lower() != 'none':
        raise ValueError('use re-method=4-points and auxbasis=none for conventional singlet TDA gradients')
    config = scf_config(settings)
    config.validate()
    minimum_gap = float(settings.get('root-gap-ev', .05))
    minimum_overlap = float(settings.get('root-min-overlap', .7))
    maxsteps = settings['opt-maxsteps']
    if (not math.isfinite(minimum_gap) or minimum_gap <= 0 or not 0 < minimum_overlap <= 1
            or type(maxsteps) is not int or maxsteps < 1):
        raise ValueError('invalid root-gap-ev, root-min-overlap, or opt-maxsteps')

    atoms = normalize_fragment(atoms)
    request = {'atoms_angstrom': atoms, 'scf': asdict(config), 'state_a': asdict(state_a),
               'state_b': asdict(state_b), 'opt_maxsteps': maxsteps, 'root_gap_ev': minimum_gap,
               'root_min_overlap': minimum_overlap,
               'versions': {pkg: version(pkg) for pkg in ('pyscf', 'geometric', 'numpy', 'scipy')},
               'implementation_sha256': {Path(p).name: file_hash(p)
                                         for p in (__file__, scf.__file__, core.__file__)}}
    key = fingerprint(request)
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    output = directory / 'reorganization.json'
    if output.exists():
        old = json.loads(output.read_text())
        if not resume or old.get('fingerprint') != key or old.get('status') != 'completed':
            raise ValueError('existing singlet task; resume only a completed calculation with identical inputs')
        for name, expected in old['output_sha256'].items():
            if not (directory / name).is_file() or file_hash(directory / name) != expected:
                raise ValueError(f'cannot resume changed or missing singlet output {name}')
        return {**old, 'reused': True}
    if any(directory.iterdir()):
        raise ValueError(f'use an empty output directory for singlet reorganization: {directory}')
    record = {'schema': 'wren-excited-reorganization-1', 'status': 'running', 'input': request,
              'fingerprint': key, 'scope': 'vacuum_restricted_scf_singlet_tda_conventional_integrals',
              'environment': 'vacuum', 'method': '4-points', 'reused': False}
    # Claim the calculation before writing any log, trace, or optimized geometry.
    with output.open('x') as initial:
        json.dump(record, initial)
    start = time.perf_counter()
    mf = None
    try:
        mf = scf.make_mean_field(atoms, config, charge=state_a.charge, spin=0, density_fit=False,
                                 log=directory / 'scf-tda.log')
        scanner = mf.as_scanner()
        root = _SingletRoot(minimum_gap, minimum_overlap)

        def ground(geom, *, gradient):
            scanner(geom)
            energy = scf.scf_record(scanner, 'cpu')['energy_hartree']
            grad = scanner.nuc_grad_method().kernel() if gradient else None
            return energy, energy, grad, {}

        def excited(geom, *, gradient):
            e0, _, _, _ = ground(geom, gradient=False)
            td, details = root.evaluate(scanner)
            energy = e0 + float(td.e[0])
            grad = td.nuc_grad_method().kernel(state=1) if gradient else None
            return e0, energy, grad, details

        ra = _optimize(mf.mol, ground, directory, 'RA', maxsteps)
        aa, ba, _, info_ra = excited(ra, gradient=False)
        rb = _optimize(ra, excited, directory, 'RB', maxsteps)
        ab, bb, _, info_rb = excited(rb, gradient=False)
        analysis = core.analyze_four_point_pair(state_a, state_b, core.FourPointEnergies(
            aa, ba, bb, ab, 'hartree'))
        relaxation_ev = {k: v * HARTREE2EV for k, v in analysis['relaxation'].items()}
        if min(relaxation_ev['on_a'], relaxation_ev['on_b']) < -1e-5:
            raise RuntimeError(f'negative relaxation exceeds tolerance: {relaxation_ev} eV')
        # Flush PySCF output before including it in the verified artifact list.
        mf.mol.stdout.flush()
        record.update(status='completed', elapsed_seconds=time.perf_counter()-start, analysis=analysis,
                      relaxation_ev=relaxation_ev, vertical_root_diagnostics={'RA': info_ra, 'RB': info_rb},
                      geometry_optimization_converged={'RA': True, 'RB': True},
                      state_policy='energy_ordered_S1_with_gap_and_adjacent_evaluation_character_checks',
                      output_sha256={p.name: file_hash(p) for p in directory.iterdir() if p != output})
        _write_result(output, record)
        return record
    except Exception as exc:
        record.update(status='failed', error={'type': type(exc).__name__, 'message': str(exc)})
        _write_result(output, record)
        raise
    finally:
        if mf is not None:
            mf.mol.stdout.close()
