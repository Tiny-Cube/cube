"""Native PySCF/GPU4PySCF H2 singlet benchmark; uses public package APIs.

Run with python validate_gpu_tda.py result.json in a GPU4PySCF environment.
Recorded hardware: NVIDIA GeForce RTX 5090, driver 615.71.09. Versions:
PySCF 2.15.0, GPU4PySCF 1.9.1, geomeTRIC 1.1.1, CuPy 13.6.0, NumPy 2.5.4.

This measures the backend's capability separately from Wren's CPU singlet
workflow. Energies in each list follow E_A_RA, E_B_RA, E_B_RB, E_A_RB.
"""
import argparse
import json
from pathlib import Path

import cupy
import gpu4pyscf
import numpy as np
from pyscf import gto, lib
from pyscf.data.nist import BOHR
from pyscf.geomopt.geometric_solver import kernel


def hydrogen(distance_angstrom):
    return gto.M(atom=[('H', (0, 0, 0)), ('H', (0, 0, distance_angstrom))],
                 basis='6-31g', verbose=0)


def solve(mol, method, gpu):
    mf = mol.RHF() if method == 'hf' else mol.RKS(xc=method)
    mf.conv_tol = 1e-12
    if gpu:
        mf = mf.to_gpu()
    mf.kernel()
    if not mf.converged:
        raise RuntimeError('SCF failed')
    td = mf.TDA().run(nstates=3, conv_tol=1e-9)
    if not all(td.converged):
        raise RuntimeError('TDA failed')
    return mf, td


def independent_gradient_check(method):
    mf, td = solve(hydrogen(.79), method, True)
    grad = td.nuc_grad_method()
    gpu_grad = float(grad.kernel(state=1)[1, 2])
    if not isinstance(td.xy[0][0], cupy.ndarray) or not type(grad).__module__.startswith('gpu4pyscf.'):
        raise RuntimeError('GPU TDA amplitude or analytical gradient provider could not be verified')

    # Independent CPU reference diagonalizes the full singlet TDA A matrix.
    reference = []
    for distance in (.789, .790, .791):
        cpu_mf, cpu_td = solve(hydrogen(distance), method, False)
        a, _ = cpu_td.get_ab()
        n = a.shape[0] * a.shape[1]
        reference.append(float(cpu_mf.e_tot + np.linalg.eigvalsh(a.reshape(n, n))[0]))
    derivative = BOHR * (reference[2] - reference[0]) / .002
    if abs(derivative - gpu_grad) >= 5e-5:
        raise AssertionError('GPU S1 gradient disagrees with finite differences')
    return {'gradient_provider': type(grad).__module__, 'amplitudes_array': type(td.xy[0][0]).__module__,
            'gpu_gradient_hartree_per_bohr': gpu_grad, 'cpu_finite_difference_hartree_per_bohr': derivative,
            'cpu_explicit_S1_hartree': reference[1], 'gpu_S1_hartree': float(mf.e_tot + td.e[0])}


def four_points(method, gpu):
    initial, _ = solve(hydrogen(.78), method, gpu)
    ok, ra = kernel(initial, maxsteps=55)
    if not ok:
        raise RuntimeError('S0 geometry optimization failed')
    ma, ta = solve(ra, method, gpu)
    ok, rb = kernel(ta.nuc_grad_method().as_scanner(state=1), maxsteps=55)
    if not ok:
        raise RuntimeError('S1 geometry optimization failed')
    mb, tb = solve(rb, method, gpu)
    return {'energies_hartree': [float(ma.e_tot), float(ma.e_tot + ta.e[0]),
                                 float(mb.e_tot + tb.e[0]), float(mb.e_tot)],
            'bond_angstrom': [float(np.linalg.norm(np.diff(m.atom_coords(unit='Angstrom'), axis=0)))
                              for m in (ra, rb)]}


def run():
    lib.num_threads(1)
    report = {}
    for method in ('hf', 'b3lyp'):
        row = independent_gradient_check(method)
        row['cpu'] = four_points(method, False)
        row['gpu'] = four_points(method, True)
        diff = np.subtract(row['gpu']['energies_hartree'], row['cpu']['energies_hartree'])
        row['max_energy_difference_hartree'] = float(np.max(np.abs(diff)))
        if row['max_energy_difference_hartree'] >= 1e-5:
            raise AssertionError('four-point total energies disagree')
        report[method] = row
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    with args.output.open('x') as file:
        json.dump(run(), file, indent=2, allow_nan=False)
        file.write('\n')
