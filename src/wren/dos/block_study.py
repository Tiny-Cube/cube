"""Reblock archived frame-resolved DOS without loading a quantum backend."""
import csv
import json
import os
from pathlib import Path
import tempfile

import numpy as np

from wren.workflow.dimers import file_hash, fingerprint
from . import sampling


def _source_data(source, expected_mode):
    result_path = source/'result.json'
    result = json.loads(result_path.read_text())
    if result.get('schema') != 'wren-dos-1' or result.get('status') != 'completed':
        raise ValueError('block study requires a completed wren-dos-1 result')
    mode = result.get('mode')
    if mode not in ('periodic', 'monomer-orbitals') or (expected_mode is not None and expected_mode != mode):
        raise ValueError('block study requires a matching monomer-orbitals or periodic DOS result')
    request = result['request']
    if fingerprint(request) != result['fingerprint'] or request['settings']['mode'] != mode:
        raise ValueError('saved DOS request fingerprint or mode does not match')
    if result.get('frames') != 'spectra.json' or result.get('energy_unit') != 'eV':
        raise ValueError('expected frame-resolved spectra.json with energies in eV')
    hashes = {'result.json': file_hash(result_path)}
    for name in ('spectra.json', 'dos.csv'):
        hashes[name] = file_hash(source/name)
        if result.get('output_sha256', {}).get(name) != hashes[name]:
            raise ValueError(f'{name}: missing or mismatched output hash')

    spectra = json.loads((source/'spectra.json').read_text())
    frames = spectra['frames']
    count = len(frames)
    if (spectra.get('schema') != 'wren-dos-spectra-1' or count != result.get('frame_count')
            or count != len(request['datasets']) or count < 2):
        raise ValueError('block study requires consistent metadata for at least two MD frames')
    indices = [f.get('frame') for f in frames]
    times = [f.get('time_ps') for f in frames]
    if (indices != [d.get('frame') for d in request['datasets']]
            or times != [d.get('time_ps') for d in request['datasets']]):
        raise ValueError('stored DOS frames or times disagree with source request order')
    if mode == 'periodic':
        reference = request['settings'].get('periodic-reference')
        if reference not in ('external', 'frame-vbm', 'cell-mean-local') or any(
                f.get('energy_alignment', {}).get('mode') != reference
                or not np.isfinite(f['energy_alignment'].get('reference_ev', float('nan'))) for f in frames):
            raise ValueError('periodic block study requires aligned per-frame energy references')
    times = None if all(t is None for t in times) else times
    sampling.validate_sampling(indices, 1, times_ps=times)

    with (source/'dos.csv').open(newline='') as stream:
        header = next(csv.reader([stream.readline()]))
        if header != ['energy_ev', 'dos', *[f'frame_{i}' for i in range(count)]]:
            raise ValueError('dos.csv needs energy_ev,dos and ordered frame_0,... columns')
        data = np.loadtxt(stream, delimiter=',', ndmin=2)
    if data.shape[1] != count+2 or len(data) < 2 or not np.isfinite(data).all():
        raise ValueError('dos.csv requires finite values on a common energy grid')
    grid, mean, curves = data[:, 0], data[:, 1], data[:, 2:].T
    spacing = result.get('grid_spacing_ev')
    if (spacing is None or not np.isfinite(spacing) or spacing <= 0 or np.any(np.diff(grid) <= 0)
            or not np.allclose(np.diff(grid), spacing, rtol=1e-8, atol=1e-12)):
        raise ValueError('dos.csv energy axis must increase with the reported spacing')
    if not np.allclose(curves.mean(axis=0), mean, rtol=1e-10, atol=1e-12):
        raise ValueError('dos.csv mean DOS disagrees with the frame columns')
    info = {'mode': mode, 'source_directory': str(source), 'source_fingerprint': result['fingerprint'],
            'source_sha256': hashes, 'frame_indices': indices, 'times_ps': times,
            'normalization': result['normalization'], 'energy_reference': result['energy_reference'],
            'broadening': result['broadening']}
    return info, grid, curves


def run_block_study(source_dir, block_sizes, *, resume=False, check_only=False, expected_mode=None):
    """Verify and analyze archived DOS. A resumed new selection keeps previous curves."""
    source = Path(source_dir).expanduser().resolve()
    sizes = list(block_sizes)
    if not sizes or any(type(n) is not int or n <= 0 for n in sizes) or len(set(sizes)) != len(sizes):
        raise ValueError('block sizes must be distinct positive integers')
    info, grid, curves = _source_data(source, expected_mode)
    indices, times = info['frame_indices'], info['times_ps']
    for size in sizes:
        sampling.validate_sampling(indices, size, times_ps=times)
    destination = source.with_name(source.name+'-blocks')
    output = destination/'block-study.json'
    if check_only:
        return {'status': 'validated', 'source': str(source), 'output': str(output),
                'frame_count': len(indices), 'block_sizes': sizes}

    implementation = {Path(p).name: file_hash(p) for p in (__file__, sampling.__file__)}
    key = fingerprint({'source_sha256': info['source_sha256'], 'implementation_sha256': implementation,
                       'numpy_version': np.__version__})
    updating = output.is_file()
    if updating:
        previous = json.loads(output.read_text())
        if (not resume or previous.get('schema') != 'wren-dos-block-study-1' or previous.get('status') != 'completed'
                or previous.get('source_and_method_fingerprint') != key):
            raise ValueError('existing block study; --resume requires the same source and analysis implementation')
        for size in previous['block_sizes']:
            if type(size) is not int or size <= 0:
                raise ValueError('invalid stored block size')
            name = f'block-{size}.csv'
            if previous['output_sha256'].get(name) != file_hash(destination/name):
                raise ValueError(f'{name}: block result missing or changed')
        if previous['block_sizes'] == sizes:
            return {**previous, 'output': str(output), 'reused': True}
    elif destination.exists():
        raise ValueError('untracked block study directory; expected a completed block-study.json')

    report = {'schema': 'wren-dos-block-study-1', 'status': 'completed', **info, 'block_sizes': sizes,
              'source_and_method_fingerprint': key, 'implementation_sha256': implementation,
              'numpy_version': np.__version__, 'studies': [], 'output_sha256': {}}
    with tempfile.TemporaryDirectory(dir=destination.parent, prefix='.wren-block-study-') as tmp:
        stage = Path(tmp)/'study'
        stage.mkdir()
        for size in sizes:
            stats = sampling.block_curve_statistics(curves, indices, size, times_ps=times)
            se = stats['standard_error']
            peak = int(np.argmax(se))
            name = f'block-{size}.csv'
            with (stage/name).open('w', newline='') as stream:
                writer = csv.writer(stream)
                writer.writerow(['energy_ev', 'retained_mean', 'standard_error'])
                writer.writerows(zip(grid, stats['retained_mean'], se))
            report['studies'].append({**stats['metadata'], 'csv': name,
                                      'peak_standard_error': float(se[peak]), 'peak_se_energy_ev': float(grid[peak])})
            report['output_sha256'][name] = file_hash(stage/name)
        (stage/output.name).write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
        if not updating:
            stage.rename(destination)
        else:
            # Keep older per-size files. Each published manifest refers only to
            # complete, immutable CSVs from a single verified source.
            for name, digest in report['output_sha256'].items():
                try:
                    os.link(stage/name, destination/name)
                except FileExistsError:
                    if file_hash(destination/name) != digest:
                        raise ValueError(f'{name}: conflicting existing block curve')
            (stage/output.name).replace(output)
    return {**report, 'output': str(output), 'reused': False}
