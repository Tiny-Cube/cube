"""Rebuild the ensemble and complete-block standard error from raw SCF bands."""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
from wren.workflow.dimers import file_hash


def check(directory, output):
    result = json.loads((directory/'result.json').read_text())
    for name, digest in result['output_sha256'].items():
        if file_hash(directory/name) != digest:
            raise ValueError(f'missing or changed {name}')
    frames = json.loads((directory/'spectra.json').read_text())['frames']
    dos = np.genfromtxt(directory/'dos.csv', delimiter=',', names=True)
    x = dos['energy_ev']
    sigma = result['broadening']['sigma_ev']
    curves, rows = [], []
    for index, frame in enumerate(frames):
        ref = frame['energy_alignment']['reference_ev']
        spectrum = np.zeros_like(x)
        for bands, w in zip(frame['orbital_energies_ev'], frame['kpoint_weights'], strict=True):
            energy = np.asarray(bands)-ref
            spectrum += 2*w*np.exp(-.5*((x[:, None]-energy)/sigma)**2).sum(axis=1)/(sigma*np.sqrt(2*np.pi))
        np.testing.assert_allclose(spectrum, dos[f'frame_{index}'], rtol=1e-12, atol=1e-12)
        curves.append(spectrum)
        rows.append({'frame': frame['frame'], 'time_ps': frame['time_ps'], 'reference_ev': ref,
                     'sampled_gap_ev': frame['band_edges']['sampled_gap_ev'],
                     'aligned_vbm_ev': frame['energy_alignment']['aligned_vbm_ev'],
                     'aligned_cbm_ev': frame['energy_alignment']['aligned_cbm_ev'],
                     'occupied_electrons': frame['scf']['occupation_integral_per_cell'],
                     'scf_seconds': frame['scf']['elapsed_seconds']})
    curves = np.asarray(curves)
    mean = curves.mean(axis=0)
    pdos = np.genfromtxt(directory/'pdos.csv', delimiter=',', names=True)
    np.testing.assert_array_equal(x, pdos['energy_ev'])
    np.testing.assert_allclose(pdos['H']+pdos['C'], mean, rtol=1e-12, atol=1e-12)
    archive = directory.with_name(directory.name+'-blocks')
    study = json.loads((archive/'block-study.json').read_text())
    errors = {}
    for item in study['studies']:
        size = item['block_size_frames']
        n = item['full_blocks']
        block = curves[:size*n].reshape(n, size, len(x)).mean(axis=1)
        se = np.sqrt(((block-block.mean(axis=0))**2).sum(axis=0)/(n*(n-1)))
        external = np.genfromtxt(archive/item['csv'], delimiter=',', names=True)
        np.testing.assert_allclose(se, external['standard_error'], rtol=1e-11, atol=1e-12)
        errors[size] = se
    summary = {'schema': 'wren-periodic-md-validation-1', 'software_commit': '2d561b7849321ba9038a25b6fe94603384b1ca11',
               'input_sha256': result['request']['source_sha256'], 'versions': result['request']['versions'],
               'device': frames[0]['scf']['device'], 'cell_volume_bohr3': frames[0]['cell_volume_bohr3'],
               'reference': result['energy_reference'], 'sigma_ev': sigma, 'frames': rows,
               'integrated_spin_states_per_cell': float(np.trapezoid(mean, x)),
               'max_total_dos_difference_from_raw_bands': float(abs(mean-dos['dos']).max()),
               'max_element_sum_residual': float(abs(pdos['H']+pdos['C']-mean).max()),
               'elapsed_seconds': result['elapsed_seconds'], 'peak_pointwise_se': {b: float(v.max()) for b, v in errors.items()}}
    output.mkdir(parents=True, exist_ok=True)
    (output/'validation.json').write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n')
    with (output/'dos.csv').open('w', newline='') as stream:
        writer = csv.writer(stream, lineterminator='\n')
        writer.writerow(['energy_ev', 'mean', 'C', 'H', *[f'se_block_{b}' for b in errors]])
        writer.writerows(zip(x, mean, pdos['C'], pdos['H'], *errors.values(), strict=True))
    return summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    print(json.dumps(check(args.directory, args.output), indent=2))
