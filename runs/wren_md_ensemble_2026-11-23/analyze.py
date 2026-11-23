"""Reproduce the finite-trajectory DOS and blocking diagnostics from saved MO levels."""
import argparse
import csv
import json
from pathlib import Path

import numpy as np
from wren.dos.core import broaden
from wren.dos.sampling import block_curve_statistics


def analyze(path, destination):
    data = json.loads(path.read_text())
    records = data['frames']
    frames = [row['frame'] for row in records]
    times = [row['time_ps'] for row in records]
    energies = np.array([row['orbital_energies_ev'] for row in records])
    occupations = np.array([row['occupations'] for row in records])
    if not (len(frames) == 8 and energies.shape == (8, 302) and occupations.shape == energies.shape):
        raise ValueError('expected eight full 302-spatial-orbital spectra')
    if not (np.isfinite(energies).all() and np.all((occupations == 0) | (occupations == 2))):
        raise ValueError('expected finite restricted closed-shell spectra')

    homo = np.max(np.where(occupations > 0, energies, -np.inf), axis=1)
    lumo = np.min(np.where(occupations == 0, energies, np.inf), axis=1)
    descriptors = np.column_stack([homo, lumo, lumo - homo])
    grid = np.linspace(-8., 6., 701)
    sigma = 0.1
    weights = np.full(energies.shape[1], 2.)
    curves = np.stack([broaden(levels, grid, sigma, weights=weights) for levels in energies])
    blocks = (1, 2, 4)
    names = ('homo_ev', 'lumo_ev', 'molecular_orbital_gap_ev')
    errors, curve_se, summaries = {}, {}, {}
    for size in blocks:
        params = {'frame_indices': frames, 'times_ps': times, 'block_size': size}
        result = block_curve_statistics(descriptors, **params)
        curve_result = block_curve_statistics(curves, **params)
        direct_se = np.std(descriptors.reshape(-1, size, 3).mean(axis=1),
                           axis=0, ddof=1) / np.sqrt(len(frames)//size)
        np.testing.assert_allclose(direct_se, result['standard_error'], rtol=1e-12, atol=1e-15)
        curve_se[size] = curve_result['standard_error']
        errors[size] = dict(zip(names, result['standard_error'].tolist()))
        summaries[size] = result['metadata']
    mean = curves.mean(axis=0)

    # Check spectral normalization independently of the finite plotted window.
    # Integrate every occupied/empty spatial MO with spin weight two on the
    # same full-spectrum grid; the plotted frontier window is not rescaled.
    full = np.arange(float(energies.min())-8*sigma, float(energies.max())+9*sigma, sigma/4)
    integrals = [float(np.trapezoid(broaden(row, full, sigma, weights=weights), full)) for row in energies]
    if not np.allclose(integrals, weights.sum(), rtol=1e-9, atol=0):
        raise ArithmeticError('full spin-state spectral sum rule failed')

    report = {
        'schema': 'wren-md-monomer-sampling-1',
        'scope': 'finite_time_series_one_complete_isolated_vacuum_monomer',
        'energy_reference': 'unshifted_molecular_orbital_energies_vacuum_boundary',
        'frame_count': len(frames), 'frame_indices': frames, 'times_ps': times,
        'molecule_id': 1, 'qm_atom_count': 86, 'spatial_orbitals_per_frame': int(energies.shape[1]),
        'spin_degeneracy': 2, 'gaussian_sigma_ev': sigma, 'plot_window_ev': [float(grid[0]), float(grid[-1])],
        'total_spin_states_expected': float(weights.sum()), 'all_orbitals_integral_by_frame': integrals,
        'dos_units': 'spin states per molecule per eV',
        'observed': {name: {'mean': float(v.mean()), 'min': float(v.min()), 'max': float(v.max()),
                            'population_standard_deviation': float(v.std(ddof=0)),
                            'sample_standard_deviation': float(v.std(ddof=1)),
                            'first_half_mean': float(v[:4].mean()), 'second_half_mean': float(v[4:].mean())}
                     for name, v in zip(names, descriptors.T)},
        'block_standard_errors_ev': errors, 'block_metadata': summaries,
        'pointwise_dos_standard_error_csv': 'frontier-dos.csv',
        'equilibration_demonstrated': False, 'block_independence_demonstrated': False,
        'material_dos_convergence_demonstrated': False,
    }
    destination.mkdir(parents=True, exist_ok=True)
    (destination/'ensemble-report.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    with (destination/'frontier-orbitals.csv').open('w', newline='') as out:
        writer = csv.writer(out, lineterminator="\n")
        writer.writerow(['frame', 'time_ps', *names])
        writer.writerows((frame, time, *row) for frame, time, row in zip(frames, times, descriptors))
    with (destination/'frontier-dos.csv').open('w', newline='') as out:
        writer = csv.writer(out, lineterminator="\n")
        writer.writerow(['energy_ev', 'mean_dos', *[f'se_block_{size}_frames' for size in blocks]])
        writer.writerows(zip(grid, mean, *(curve_se[size] for size in blocks)))
    return report, grid, curves, curve_se[2], descriptors


def plot(report, grid, curves, uncertainty, descriptors, destination):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, (spectrum, levels) = plt.subplots(1, 2, figsize=(9.6, 3.5), layout='constrained')
    mean = curves.mean(axis=0)
    for row in curves:
        spectrum.plot(grid, row, color='#83949d', alpha=.24, lw=.8)
    spectrum.fill_between(grid, mean - uncertainty, mean + uncertainty, color='#2e728e', alpha=.25,
                          label='Mean ± SE (four blocks)')
    spectrum.plot(grid, mean, color='#185570', lw=1.6, label='8-frame mean')
    spectrum.set(xlabel='Orbital energy (eV)', ylabel='Spin states / molecule / eV', xlim=(grid[0], grid[-1]))
    spectrum.legend(fontsize=8)
    for column, title, color, marker in [(0, 'HOMO', '#166486', 'o'), (1, 'LUMO', '#b45325', 's')]:
        values = descriptors[:, column]
        levels.plot(np.asarray(report['times_ps'])/1000, values-values.mean(), label=title,
                    color=color, marker=marker, lw=1.3, ms=4.5)
    levels.axhline(0, color='#a1abae', lw=.8, zorder=0)
    levels.set(xlabel='MD time (ns)', ylabel='Level minus its series mean (eV)')
    levels.legend(fontsize=8)
    fig.savefig(destination/'ensemble-dos.png', dpi=180)
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('spectra', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--plot', action='store_true')
    args = parser.parse_args()
    results = analyze(args.spectra, args.output)
    if args.plot:
        plot(*results, args.output)
