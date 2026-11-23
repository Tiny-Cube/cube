#!/usr/bin/env python3
"""Compare two completed Wren periodic DOS runs of the same model and k mesh."""
import argparse
import json
from pathlib import Path

import numpy as np

from wren.workflow.dimers import file_hash


def load(directory):
    result = json.loads((directory/'result.json').read_text())
    if result['status'] != 'completed' or result['mode'] != 'periodic':
        raise ValueError(f'{directory} requires a completed periodic DOS calculation')
    for name, digest in result['output_sha256'].items():
        if file_hash(directory/name) != digest:
            raise ValueError(f'changed output: {directory/name}')
    frame, = json.loads((directory/result['frames']).read_text())['frames']
    return result, frame


def curves(directory, name):
    data = np.genfromtxt(directory/name, names=True, delimiter=',')
    return data['energy_ev'], {key: data[key] for key in data.dtype.names[1:] if not key.startswith('frame_')}


def compare(first, second):
    results, frames = zip(load(first), load(second))
    models = []
    for result in results:
        request = result['request']
        settings = {k: v for k, v in request['settings'].items()
                    if k not in ('input', 'output', 'qc-engine', 'dos-backend', 'threads', 'max-memory-mb')}
        models.append((settings, request['source_sha256'], request['datasets']))
    if models[0] != models[1]:
        raise ValueError('controls or crystal structures differ')
    for key in ('kpoints_fractional', 'kpoint_weights', 'occupations'):
        if frames[0][key] != frames[1][key]:
            raise ValueError(f'unmatched {key}')
    delta_bands = max(float(np.max(np.abs(np.subtract(a, b)))) for a, b in
                      zip(*(frame['orbital_energies_ev'] for frame in frames), strict=True))
    files = ['dos.csv']
    if 'element_projections' in results[0]:
        files.append('pdos.csv')
    residuals = {}
    for name in files:
        grid, left = curves(first, name)
        other_grid, right = curves(second, name)
        if left.keys() != right.keys():
            raise ValueError(f'{name} has unmatched columns')
        # Defaults may place endpoints at slightly different energies on different devices.
        samples = np.union1d(grid, other_grid)
        for key in left:
            a = np.interp(samples, grid, left[key], left=0., right=0.)
            b = np.interp(samples, other_grid, right[key], left=0., right=0.)
            residuals[f'{name}/{key}'] = float(np.max(abs(a-b)))
    return {
        'controls': models[0][0], 'input_sha256': models[0][1],
        'runs': [{'fingerprint': result['fingerprint'], 'versions': result['request']['versions'],
                  'scf': frame['scf'], 'elapsed_seconds': result['elapsed_seconds'],
                  'total_energy_hartree': frame['total_energy_hartree'], 'band_edges': frame['band_edges'],
                  'integrated_states_per_cell': result['finite_grid_integral']} for result, frame in zip(results, frames)],
        'abs_total_energy_difference_hartree': abs(frames[0]['total_energy_hartree']-frames[1]['total_energy_hartree']),
        'max_abs_orbital_energy_difference_ev': delta_bands,
        'max_abs_interpolated_curve_difference_per_ev_per_cell': residuals}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('first', type=Path)
    parser.add_argument('second', type=Path)
    args = parser.parse_args()
    print(json.dumps(compare(args.first, args.second), indent=2, allow_nan=False))
