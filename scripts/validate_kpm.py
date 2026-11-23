#!/usr/bin/env python3
"""Run sparse DOS on a seeded disordered lattice and optionally an MD PySCF checkpoint."""
import argparse
import json
from pathlib import Path
import platform
import time

import numpy as np
import scipy

from wren.dos.control import normalize
from wren.dos.runner import run
from wren.workflow.dimers import file_hash


def solve(root, label, sites, pairs, controls, backend):
    network = root/f'{label}.json'
    network.write_text(json.dumps({'schema': 'wren-dos-tb-1',
        'coupling_convention': 'signed_real_consistent_orbital_phases',
        'site_energies_ev': sites.tolist(), 'pairs': pairs}, separators=(',', ':')))
    settings = normalize({'input': network, 'output': label, 'mode': 'tight-binding',
                          'solver': 'kpm', 'dos-backend': backend, **controls}, root)
    result = run(settings)
    moments = np.genfromtxt(root/label/'moments.csv', delimiter=',', names=True)
    summary = {key: result[key] for key in ('expected_integral', 'finite_grid_integral', 'integral_relative_error',
                                             'elapsed_seconds', 'kpm', 'broadening')}
    summary['network_sha256'] = file_hash(network)
    if result['integral_relative_error'] > 1e-10:
        raise RuntimeError('full-support integrated DOS disagrees with state count')
    return summary, moments


def lattice(size):
    # Each undirected bond in a square periodic layer occurs once, including wrapped bonds.
    n = size**2
    sites = np.random.default_rng(20261123).normal(0., .05, n)
    tx, ty = -.10, -.035
    pairs = []
    for y in range(size):
        for x in range(size):
            i = y*size+x
            pairs.append([i, y*size+(x+1)%size, [int(x+1 == size), 0, 0], tx])
            pairs.append([i, ((y+1)%size)*size+x, [0, int(y+1 == size), 0], ty])
    return sites, pairs


def lattice_case(root, size, backend):
    sites, pairs = lattice(size)
    n = len(sites)
    tx, ty = pairs[0][-1], pairs[1][-1]
    summary, rows = solve(root, 'disordered-layer', sites, pairs,
                          {'kpm-moments': 128, 'kpm-vectors': 16, 'kpm-seed': 314159}, backend)
    a, b = summary['kpm']['scale_ev'], summary['kpm']['center_ev']
    # At Gamma, Tr(H-bI)^2 is determined by the diagonal and the two bond families.
    first = np.sum(sites-b)/a
    second = 2*(np.sum((sites-b)**2)+2*n*(tx**2+ty**2))/a**2-n
    errors = rows['trace_moment'][1:3] - [first, second]
    se = rows['trace_standard_error'][1:3]
    if np.any(abs(errors) > 6*se+1e-8):
        raise RuntimeError('stochastic trace moments disagree with direct lattice identities')
    summary.update(sites=n, bonds=len(pairs), standard_errors_for_m1_m2=se.tolist(),
                   errors_in_standard_errors_for_m1_m2=(errors/se).tolist(), onsite_sd_ev=.05,
                   hopping_x_ev=tx, hopping_y_ev=ty, structure='synthetic periodic square layer, Gamma')
    return summary


def convergence_case(root, backend):
    """Separate finite polynomial width from probe noise on one unchanged network."""
    from wren.dos.kpm import bin_densities
    sites, pairs = lattice(16)
    h = np.diag(sites)
    for i, j, _, value in pairs:
        h[i, j] += value
        h[j, i] += value
    energies = np.linalg.eigvalsh(h)
    reports = []
    for order, vectors in ((64, 8), (64, 64), (128, 8), (128, 64)):
        label = f'convergence-m{order}-r{vectors}'
        summary, _ = solve(root, label, sites, pairs, {'kpm-moments': order,
            'kpm-vectors': vectors, 'kpm-seed': 314159, 'points': 1024}, backend)
        bins = np.genfromtxt(root/label/'dos.csv', names=True, delimiter=',')
        edges = np.r_[bins['left_ev'], bins['right_ev'][-1]]
        a, b = summary['kpm']['scale_ev'], summary['kpm']['center_ev']
        exact = np.cos(np.outer(np.arange(order), np.arccos((energies-b)/a))).sum(axis=1)
        ref = np.empty(len(bins))
        for sl, density, _ in bin_densities(exact[None, :], edges, b, a):
            ref[sl] = density
        widths = np.diff(edges)
        # Maximum CDF differences at the right edges, in fractions of all states.
        ref_cdf = np.cumsum(ref*widths)/len(energies)
        stochastic_cdf = np.cumsum(bins['dos']*widths)/len(energies)
        exact_cdf = np.searchsorted(energies, edges[1:], side='right')/len(energies)
        reports.append({'moments': order, 'vectors': vectors,
                        'max_cdf_difference_stochastic_vs_exact_jackson': float(np.max(abs(stochastic_cdf-ref_cdf))),
                        'max_cdf_difference_exact_jackson_vs_discrete': float(np.max(abs(ref_cdf-exact_cdf))),
                        'midband_resolution_scale_ev': np.pi*a/(order+1)})
    return {'sites': len(sites), 'mesh': 'Gamma', 'comparisons': reports}


def checkpoint_case(root, path, backend):
    from pyscf.data.nist import HARTREE2EV
    from pyscf.scf import chkfile
    from scipy.linalg import sqrtm
    mol, scf = chkfile.load_scf(str(path))
    c, eig = np.asarray(scf['mo_coeff']), np.asarray(scf['mo_energy'])*HARTREE2EV
    if c.shape != (mol.nao_nr(), mol.nao_nr()) or not np.isrealobj(c) or eig.shape != (len(c),):
        raise ValueError('validation requires a full real restricted molecular orbital basis')
    u = sqrtm(mol.intor_symmetric('int1e_ovlp')) @ c
    if not np.allclose(u.T @ u, np.eye(len(c)), rtol=0, atol=1e-8):
        raise ValueError('checkpoint orbitals fail the overlap metric check')
    h = (u*eig) @ u.T
    pairs = [[i, j, [0, 0, 0], h[i, j]] for i in range(len(c)) for j in range(i+1, len(c)) if h[i, j]]
    summary, rows = solve(root, 'md-monomer', h.diagonal(), pairs,
                          {'kpm-moments': 128, 'kpm-trace': 'exact', 'spin-degeneracy': 2}, backend)
    a, b = summary['kpm']['scale_ev'], summary['kpm']['center_ev']
    angles = np.arccos((eig-b)/a)
    exact = 2*np.cos(np.outer(rows['order'], angles)).sum(axis=1)
    error = float(np.max(np.abs(rows['trace_moment']-exact)))
    if error > 1e-7:
        raise RuntimeError('KPM moments disagree with checkpoint orbital eigenvalues')
    summary.update(atoms=mol.natm, spatial_orbitals=len(c), checkpoint_sha256=file_hash(path),
                   max_absolute_moment_error=error, scope='complete individual molecular Fock spectrum in symmetric Lowdin AOs')
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--size', type=int, default=256, help='square lattice edge in sites')
    parser.add_argument('--checkpoint', type=Path)
    parser.add_argument('--backend', choices=['cpu', 'gpu'], default='cpu')
    args = parser.parse_args()
    root = args.output.resolve()
    if args.size < 3 or root.exists() and any(root.iterdir()):
        parser.error('provide a new/empty output directory and a lattice edge of at least 3')
    root.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    report = {'backend': args.backend, 'platform': platform.platform(), 'numpy': np.__version__,
              'scipy': scipy.__version__}
    if args.backend == 'gpu':
        import cupy as cp
        device = cp.cuda.runtime.getDeviceProperties(cp.cuda.runtime.getDevice())['name']
        report.update(device=device.decode() if isinstance(device, bytes) else str(device),
                      cupy=cp.__version__, cuda_runtime=cp.cuda.runtime.runtimeGetVersion())
    report['disordered_layer'] = lattice_case(root, args.size, args.backend)
    report['convergence'] = convergence_case(root, args.backend)
    if args.checkpoint:
        report['md_checkpoint'] = checkpoint_case(root, args.checkpoint, args.backend)
    report['wall_seconds'] = time.perf_counter()-started
    (root/'validation.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(root/'validation.json')


if __name__ == '__main__':
    main()
