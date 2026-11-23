"""Compare a vacuum H2 four-point calculation with the vertical-gradient approximation.

From a Wren environment with PySCF, NumPy, SciPy and geomeTRIC, run:
    python runs/wren_state_pairs_2026-11-23/validate_h2.py build/h2-four-point-vg

Use a new output directory. The full quantum/optimization traces are created
there; only the input, harmonic modes and compact result accompany this script.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from pyscf import gto, scf
from pyscf.data.nist import HARTREE2EV, BOHR
from scipy import constants

from wren.workflow.excited_reorganization import calculate_excited_reorganization
from wren.workflow.vertical_gradient import read_vertical_gradient


def run(output):
    settings = {'qc-engine': 'pyscf', 'environment': 'vacuum', 're-method': '4-points',
                'method': 'hf', 'basis-set': '6-31g', 'auxbasis': 'none', 'threads': 1,
                'grid-level': 3, 'conv-tol': 1e-10, 'max-cycle': 100,
                'max-memory-mb': 1200, 'opt-maxsteps': 55}
    xyz = Path(__file__).with_name('h2.xyz')
    output.mkdir(parents=True, exist_ok=False)
    qc = output / 'four-point'
    result = calculate_excited_reorganization(xyz, settings, qc)

    mol = gto.M(atom=str(qc / 'RA.xyz'), basis=settings['basis-set'], verbose=0)
    mf = scf.RHF(mol).run(conv_tol=1e-12)
    td = mf.TDA().run(nstates=3, conv_tol=1e-9)
    if not mf.converged or not np.all(td.converged):
        raise RuntimeError('independent RHF/TDA calculation did not converge')
    grad = td.nuc_grad_method().kernel(state=1) * HARTREE2EV / BOHR
    hess = mf.Hessian().kernel().transpose(0, 2, 1, 3).reshape(6, 6) * HARTREE2EV / BOHR**2
    masses = np.full(2, 1.00782503223)  # isotopic mass (1H) in u
    mass_xyz = np.repeat(masses, 3)
    curvature, vectors = np.linalg.eigh(hess / np.sqrt(np.outer(mass_xyz, mass_xyz)))
    stretch = vectors[:, -1]  # diatomic bond stretch
    # Choose the arbitrary mode sign so H2 moves along +z as the bond elongates.
    if stretch[-1] < 0:
        stretch *= -1
    angular_frequency = np.sqrt(curvature[-1] * constants.electron_volt / (constants.atomic_mass * 1e-20))
    wavenumber = angular_frequency / (2 * np.pi * constants.c * 100)
    modes = {'schema': 'wren-vertical-gradient-1',
             'method': 'RHF/6-31G analytical reference Hessian and singlet TDA target gradient',
             'environment': 'vacuum', 'states': result['analysis']['states'],
             'reference': 'a', 'reference_stationary': True, 'mode_subspace': 'H2 bond stretch',
             'cartesian': {'atom_ids': ['H1', 'H2'], 'masses_u': masses.tolist(),
                           'gradient_ev_per_angstrom': grad.tolist()},
             'modes': [{'id': 'stretch', 'wavenumber_cm1': float(wavenumber),
                        'mass_weighted_direction': stretch.reshape(2, 3).tolist()}]}
    vg_input = output / 'h2-vg.json'
    vg_input.write_text(json.dumps(modes, indent=2, allow_nan=False) + '\n')
    vg = read_vertical_gradient(vg_input)
    approx = vg['relaxation_ev']
    direct = .5 * (stretch @ (grad.ravel() / np.sqrt(mass_xyz)))**2 / curvature[-1]
    if not np.isclose(approx, direct, rtol=1e-8, atol=1e-10):
        raise AssertionError(f'harmonic analytical check failed: {approx} vs {direct} eV')

    geometries = {}
    for label in ('RA', 'RB'):
        trace = json.loads((qc / f'{label}.trace.json').read_text())
        coords = trace[-1]['coords_angstrom']
        geometries[label] = {'coordinates_angstrom': coords, 'bond_angstrom': float(np.linalg.norm(
            np.subtract(coords[0], coords[1]))), 'evaluations': len(trace),
            'final_gradient_rms_hartree_per_bohr': trace[-1]['gradient_rms_hartree_per_bohr']}
    rb_trace = json.loads((qc / 'RB.trace.json').read_text())
    summary = {'system': 'H2, neutral closed shell, vacuum; numerical validation of these methods',
               'settings': settings, 'versions': result['input']['versions'],
               'four_point': {'states': result['analysis']['states'],
                              'energies_hartree': result['analysis']['energies'],
                              'relaxation_ev': result['relaxation_ev'], 'geometries': geometries,
                              'converged': result['geometry_optimization_converged'],
                              'minimum_tested_gap_ev': min(s['nearest_gap_ev'] for s in rb_trace),
                              'minimum_adjacent_S1_amplitude_overlap': min(
                                  s['transition_amplitude_overlaps'][0] for s in rb_trace)},
               'vertical_gradient': {'reference': 'S0', 'gradient_state': 'S1',
                                     'reference_max_abs_gradient_hartree_per_bohr': float(
                                         np.max(np.abs(mf.nuc_grad_method().kernel()))),
                                     'input': vg_input.name, 'modes': vg['modes'],
                                     'relaxation_on_S1_ev': approx,
                                     'direct_inverse_curvature_ev': float(direct)},
               'fully_optimized_on_S1_minus_harmonic_ev': result['relaxation_ev']['on_b'] - approx}
    (output / 'h2-validation.json').write_text(json.dumps(summary, indent=2, allow_nan=False) + '\n')
    return summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    summary = run(parser.parse_args().output)
    print('H2 vertical-gradient S1 relaxation (eV)', summary['vertical_gradient']['relaxation_on_S1_ev'])
    print('H2 four-point S1 relaxation (eV)', summary['four_point']['relaxation_ev']['on_b'])
