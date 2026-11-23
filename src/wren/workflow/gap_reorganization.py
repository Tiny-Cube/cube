"""Classical two-ensemble linear response from paired, vertical total energies."""
import csv
import hashlib
import io
import math
from pathlib import Path
from statistics import fmean, stdev, variance

KB_EV_K = 8.617333262145179e-5


def _ensemble(path, block_size, temperature):
    path = Path(path)
    raw = path.read_bytes()
    reader = csv.DictReader(io.StringIO(raw.decode('utf-8-sig')))
    if reader.fieldnames != ['frame', 'E_A_eV', 'E_B_eV']:
        raise ValueError(f'{path}: expected columns frame,E_A_eV,E_B_eV')
    gaps = []
    last = None
    stride = None
    for row in reader:
        try:
            if None in row:
                raise ValueError
            frame = int(row['frame'])
            a, b = float(row['E_A_eV']), float(row['E_B_eV'])
            gap = b - a
            if (frame < 0 or not math.isfinite(gap) or
                    last is not None and (frame <= last or stride is not None and frame - last != stride)):
                raise ValueError
        except (ValueError, TypeError) as exc:
            raise ValueError(f'{path}:{reader.line_num}: expected uniformly spaced ordered frame IDs '
                             'and two finite total energies in eV') from exc
        if last is not None:
            stride = frame - last
        last = frame
        gaps.append(gap)
    if len(gaps) < 2:
        raise ValueError(f'{path}: at least two equilibrium samples are required')
    mean = fmean(gaps)
    var = variance(gaps)
    # Use all samples for point estimates. Error estimate uses complete,
    # equal-size blocks and needs two or more blocks in each independent run.
    nblocks = len(gaps) // block_size if block_size else 0
    block_means = [fmean(gaps[i*block_size:(i+1)*block_size]) for i in range(nblocks)]
    if block_size and nblocks < 2:
        raise ValueError(f'{path}: block-size must leave at least two complete blocks')
    return {'path': str(path), 'sha256': hashlib.sha256(raw).hexdigest(), 'n': len(gaps),
            'frame_stride': stride, 'mean_gap_ev': mean, 'sample_variance_ev2': var,
            'lambda_variance_ev': var / (2 * KB_EV_K * temperature),
            'full_blocks': nblocks, 'samples_outside_full_blocks': len(gaps) % block_size if block_size else None,
            'block_mean_se_ev': stdev(block_means)/math.sqrt(nblocks) if nblocks else None}


def analyze_gap_ensembles(settings):
    a, b = (_ensemble(settings['ensemble-' + key], settings.get('block-size'), settings['temperature'])
            for key in ('a', 'b'))
    if Path(a['path']).samefile(b['path']):
        raise ValueError('linear response requires separate equilibrium ensembles for A and B')
    stokes = (a['mean_gap_ev'] - b['mean_gap_ev']) / 2
    free_energy = (a['mean_gap_ev'] + b['mean_gap_ev']) / 2
    error = None
    if a['block_mean_se_ev'] is not None:
        error = math.hypot(a['block_mean_se_ev'], b['block_mean_se_ev']) / 2
    return {'schema': 'wren-gap-reorganization-1', 'status': 'completed', 'method': 'linear-response',
            'states': {k: settings['state-'+k] for k in ('a', 'b')}, 'temperature_k': settings['temperature'],
            'environment': settings['environment'], 'gap_definition': 'E_B(R) - E_A(R)',
            'ensembles': {'a': a, 'b': b}, 'lambda_stokes_ev': stokes, 'delta_g_lr_ev': free_energy,
            'lambda_stokes_block_se_ev': error, 'block_size_samples': settings.get('block-size'),
            'diagnostics': {'negative_lambda_stokes': stokes < 0,
                            'variance_minus_stokes_ev': {k: v['lambda_variance_ev'] - stokes
                                                         for k, v in [('a', a), ('b', b)]},
                            'linear_response_verified': False,
                            'equilibration_and_block_independence_verified': False}}
