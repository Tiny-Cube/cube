"""Reproduce the 100-site e_006 frontier bands with the installed wren-dos CLI."""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys

import numpy as np


HERE = Path(__file__).resolve().parent
ARCHIVE = HERE.parent/'e_006_ether_12-17_0'
NSITES = 100
SIGMA = .005  # eV
MESHES = (1, 2, 4)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def dump(path, data):
    path.write_text(json.dumps(data, indent=2, allow_nan=False)+'\n')


def onsite_values(csv_path, reference):
    if csv_path is None:
        if reference is not None:
            raise ValueError('--site-reference requires --sites')
        return np.zeros((NSITES, 2)), {'model': 'uniform_zero_both_channels',
                                       'reference': 'independent arbitrary zero for each frontier manifold'}
    if not reference or not reference.strip():
        raise ValueError('declare a shared electron-energy reference in --site-reference')
    result = {}
    with csv_path.open(newline='', encoding='utf-8-sig') as stream:
        rows = csv.DictReader(stream)
        if rows.fieldnames != ['mol', 'homo_ev', 'lumo_ev']:
            raise ValueError('site table must have header mol,homo_ev,lumo_ev with electron energies in eV')
        for row in rows:
            if None in row or None in row.values():
                raise ValueError('incorrect number of site columns')
            mol = int(row['mol'])
            if not 1 <= mol <= NSITES or mol in result:
                raise ValueError(f'invalid or duplicate molecule {mol}')
            values = float(row['homo_ev']), float(row['lumo_ev'])
            if not np.isfinite(values).all() or values[0] > values[1]:
                raise ValueError(f'nonfinite or inverted HOMO/LUMO energies on molecule {mol}')
            result[mol] = values
    missing = sorted(set(range(1, NSITES+1)) - result.keys())
    if missing:
        raise ValueError(f'missing site energies for molecules {missing}')
    return np.array([result[i] for i in range(1, NSITES+1)]), {
        'model': 'complete_user_supplied_electron_energies', 'reference': reference.strip(),
        'units': 'eV', 'source_sha256': sha(csv_path), 'reference_supplied_by_user': True}


def material():
    audit = json.loads((ARCHIVE/'input_audit.json').read_text())
    box = np.asarray(audit['box_nm'], dtype=float)
    if audit['molecules'] != NSITES or box.shape != (3,) or not np.all(box > 0):
        raise ValueError('unexpected e_006 box or molecule count')
    with (ARCHIVE/'sites.csv').open(newline='') as stream:
        centers = {int(r['mol_id']): np.array([float(r[f'core_{c}_nm']) for c in 'xyz'])
                   for r in csv.DictReader(stream)}
    if set(centers) != set(range(1, NSITES+1)):
        raise ValueError('archived molecule IDs differ from 1..100')
    pairs = {'homo': [], 'lumo': []}
    seen = set()
    with (ARCHIVE/'rates_and_geometry.csv').open(newline='') as stream:
        for r in csv.DictReader(stream):
            i, j = int(r['mol_i']), int(r['mol_j'])
            image = [int(r['image_'+axis]) for axis in 'abc']
            delta = [float(r[f'd{c}_nm']) for c in 'xyz']
            if (i, j) in seen or not 1 <= i < j <= NSITES:
                raise ValueError('unexpected duplicate or reversed e_006 pair')
            seen.add((i, j))
            if not np.allclose(centers[j]-centers[i]+image*box, delta, rtol=0, atol=1e-10):
                raise ValueError(f'periodic image and measured displacement disagree for {(i, j)}')
            for name, column in [('homo', 'VH_meV'), ('lumo', 'VL_meV')]:
                value = float(r[column])/1000.
                if not math.isfinite(value):
                    raise ValueError(f'invalid coupling on {(i, j)}')
                pairs[name].append([i-1, j-1, image, value])
    if len(seen) != audit['pairs'] or len(seen) != 550:
        raise ValueError('incomplete e_006 coupling archive')
    return pairs, np.diag(10*box).tolist()


def network(energies, pairs, convention, lattice, level_info):
    terms = []
    for i, j, image, value in pairs:
        if convention == 'flipped':
            # One fixed choice per orbital, including all its periodic images.
            value *= (-1 if i % 3 == 0 else 1)*(-1 if j % 3 == 0 else 1)
        elif convention == 'absolute':
            value = abs(value)
        terms.append([i,j,image,value])
    return {'schema': 'wren-dos-tb-1', 'coupling_convention': 'signed_real_consistent_orbital_phases',
            'site_energies_ev': energies.tolist(), 'pairs': terms,
            'lattice_angstrom': lattice, 'energy_reference': level_info['reference'],
            'model_scope': {'signed': 'assume common site gauge in archived MOMAP calculation',
                            'flipped': 'gauge-equivalent demonstration using fixed site signs',
                            'absolute': 'counterfactual different Hamiltonian with all positive couplings'}[convention]}


def run_case(root, channel, convention, mesh, grid):
    case = root/f'{channel}-{convention}-k{mesh}'
    case.mkdir()
    emin, emax, points = grid
    control = case/'dos.toml'
    control.write_text(f'''[dos]
mode = "tight-binding"
input = "../{channel}-{convention}.json"
output = "results"
kmesh = [{mesh}, {mesh}, {mesh}]
spin-degeneracy = 2
sigma = {SIGMA}
points = {points}
emin = {emin:.17g}
emax = {emax:.17g}
''')
    executable = Path(sys.executable).with_name('wren-dos')
    cmd = str(executable) if executable.is_file() else 'wren-dos'
    subprocess.run([cmd, str(control)], check=True,
                   env={**os.environ, 'OPENBLAS_NUM_THREADS':'1', 'OMP_NUM_THREADS':'1', 'MKL_NUM_THREADS':'1'})
    out = case/'results'
    result = json.loads((out/'result.json').read_text())
    for name in ('spectra.json','dos.csv'):
        if sha(out/name) != result['output_sha256'][name]:
            raise RuntimeError(f'Wren output integrity failed: {out/name}')
    frame = json.loads((out/'spectra.json').read_text())['frames'][0]
    levels = np.asarray(frame['eigenvalues_ev'])
    if (result['status'] != 'completed' or levels.shape != (mesh**3, NSITES)
            or not np.allclose(frame['kpoint_weights'], 1./mesh**3, rtol=0, atol=1e-14)):
        raise RuntimeError('incomplete or incorrectly weighted frontier spectrum')
    axis, density = np.loadtxt(out/'dos.csv', delimiter=',', skiprows=1, usecols=(0,1), unpack=True)
    density /= NSITES
    integral = float(np.trapezoid(density, axis))
    if abs(integral-2.) > 1e-8 or result['expected_integral'] != 2*NSITES:
        raise RuntimeError('frontier state count failed')
    metrics = {'kpoints': mesh**3, 'min_ev': float(levels.min()), 'max_ev': float(levels.max()),
               'bandwidth_ev': float(np.ptp(levels)), 'integral_states_per_molecule': integral}
    return axis, density, levels, metrics


def compare(axis, a, b):
    return float(.25*np.trapezoid(np.abs(a-b), axis))  # TV between two unit-normalized densities


def channel_report(root, channel, energies, pairs, lattice, level_info):
    # All variants have the same Gershgorin bounds. Gaussian tails beyond ±9σ
    # contribute negligibly; all meshes are compared on one energy axis.
    radius = np.zeros(NSITES)
    for i, j, _, v in pairs:
        radius[[i,j]] += abs(v)
    lo, hi = float(min(energies-radius) - 9*SIGMA), float(max(energies+radius) + 9*SIGMA)
    points = math.ceil((hi-lo)/(.2*SIGMA))+1
    if points > 20000:
        raise ValueError('site-energy range exceeds 20,000 grid points; check declared eV units')
    grid = lo, hi, points
    raw = root/'work'
    raw.mkdir(exist_ok=True)
    axes = []
    curves = {}
    spectra = {}
    mesh_report = {}
    for convention in ('signed', 'flipped', 'absolute'):
        dump(raw/f'{channel}-{convention}.json', network(energies, pairs, convention, lattice, level_info))
        for mesh in MESHES if convention == 'signed' else (MESHES[-1],):
            key = f'{convention}_k{mesh}'
            x, curves[key], spectra[key], metrics = run_case(raw, channel, convention, mesh, grid)
            axes.append(x)
            mesh_report[key] = metrics
    x = axes[0]
    if not all(np.array_equal(x, y) for y in axes):
        raise RuntimeError('inconsistent energy grids')
    # There is exactly one edge per unordered pair. The first two spectral
    # moments at each k can therefore be checked without using an eigensolver.
    m1 = float(np.mean(energies))
    m2 = float(np.mean(energies**2) + 2*sum(p[3]**2 for p in pairs)/NSITES)
    for case, levels in spectra.items():
        if not np.allclose(np.mean(levels, axis=1), m1, rtol=1e-12, atol=1e-12):
            raise RuntimeError(f'{case} fails the analytic first spectral moment')
        if not np.allclose(np.mean(levels**2, axis=1), m2, rtol=1e-12, atol=1e-12):
            raise RuntimeError(f'{case} fails the analytic second spectral moment')
    finest = 'signed_k4'
    signed = spectra[finest]
    gauge_error = float(np.max(np.abs(signed-spectra['flipped_k4'])))
    gauge_curve_error = float(np.max(np.abs(curves[finest]-curves['flipped_k4'])))
    if gauge_error > 1e-9 or gauge_curve_error > 1e-8:
        raise RuntimeError('site gauge transformation changed the spectrum')
    keys = ['signed_k1', 'signed_k2', finest, 'absolute_k4']
    with (root/f'{channel}.csv').open('w', newline='') as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(['energy_ev', *keys])
        writer.writerows(zip(x, *(curves[k] for k in keys), strict=True))
    return {'spectra': mesh_report, 'mesh_tv_k1_to_k2': compare(x, curves['signed_k1'], curves['signed_k2']),
            'mesh_tv_k2_to_k4': compare(x, curves['signed_k2'], curves[finest]),
            'site_flip_max_band_change_ev': gauge_error, 'site_flip_max_dos_change': gauge_curve_error,
            'absolute_max_band_change_ev': float(np.max(np.abs(signed-spectra['absolute_k4']))),
            'absolute_dos_tv': compare(x, curves[finest], curves['absolute_k4']),
            'first_moment_ev': m1, 'second_moment_ev2': m2,
            'sha256_csv': sha(root/f'{channel}.csv'), 'grid_spacing_ev': float(x[1]-x[0])}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sites', type=Path, help='complete mol,homo_ev,lumo_ev table')
    parser.add_argument('--site-reference', help='common electronic energy zero and level of theory for all 100 sites')
    parser.add_argument('--output', type=Path, default=HERE.parents[1]/'build/e006-frontier-dos',
                        help='new directory for CSVs, summary, Wren controls and raw outputs')
    args = parser.parse_args()
    energies, source = onsite_values(args.sites, args.site_reference)
    terms, lattice = material()
    output = args.output.resolve()
    if output.exists():
        parser.error(f'output exists; choose a new directory: {output}')
    output.mkdir(parents=True)
    report = {'material': 'e_006_ether_12-17_0', 'frame': 2000, 'molecules': NSITES, 'pairs':550,
              'model': 'one orthonormal frontier orbital per molecule per channel, fixed periodic morphology',
              'source_sha256':{p: sha(ARCHIVE/p) for p in ('sites.csv','rates_and_geometry.csv','input_audit.json')},
              'script_sha256': sha(Path(__file__)), 'onsite':source, 'lattice_angstrom':lattice, 'sigma_ev':SIGMA, 'spin_degeneracy':2,
              'energy_axis':'electron energies in eV', 'csv_density_units':'states / (eV molecule)',
              'gauge_from_original_checkpoints_independently_verified':False, 'channels':{}}
    for idx, name in enumerate(('homo','lumo')):
        report['channels'][name] = channel_report(output,name,energies[:,idx],terms[name],lattice,source)
    dump(output/'report.json',report)
    print(output/'report.json')


if __name__ == '__main__':
    main()
