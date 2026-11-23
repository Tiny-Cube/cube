"""Compare periodic DOS on successively refined meshes of the same crystal."""
import json
from pathlib import Path

import numpy as np

from wren.quantum.cli import _write_result
from .core import broaden


def lines(frame):
    bands = [np.asarray(row) for row in frame['orbital_energies_ev']]
    weights = np.concatenate([np.full(len(row), 2*w) for row, w in
                              zip(bands, frame['kpoint_weights'], strict=True)])
    return np.concatenate(bands), weights


def shape_change(previous, current, sigma):
    """Total variation of normalized Gaussian DOS on a common unshifted energy grid."""
    e1, w1 = lines(previous)
    e2, w2 = lines(current)
    if not np.isclose(w1.sum(), w2.sum(), atol=1e-8, rtol=1e-10):
        raise ValueError('mesh study needs the same number of spin states per cell')
    lo, hi = min(e1.min(), e2.min())-8*sigma, max(e1.max(), e2.max())+8*sigma
    count = int(np.ceil((hi-lo)/(sigma/5)))+1
    if count > 1_000_000:
        raise ValueError('mesh study needs over a million energy samples; increase sigma')
    grid = np.linspace(lo, hi, count)
    d1 = broaden(e1, grid, sigma, weights=w1)
    d2 = broaden(e2, grid, sigma, weights=w2)
    return float(np.trapezoid(abs(d1-d2), grid)/(2*w1.sum()))


def run_study(settings, meshes, *, resume=False):
    from .runner import run
    from pyscf.data.nist import HARTREE2EV
    single = Path(settings['output'])
    root = single.with_name(single.name+'-kmesh')
    output = root/'kmesh-study.json'
    if (root/'result.json').exists():
        raise ValueError('choose a study output separate from a single DOS run')
    if output.exists() and not resume:
        raise ValueError('existing mesh study; use --resume or choose another output')
    rows, previous = [], None
    for mesh in meshes:
        folder = root/('mesh-'+'x'.join(map(str, mesh)))
        result = run(dict(settings, kmesh=mesh, output=str(folder)), resume=resume)
        frame, = json.loads((folder/result['frames']).read_text())['frames']
        row = {'kmesh': mesh, 'result': str(folder.relative_to(root)/'result.json'),
               'fingerprint': result['fingerprint'], 'total_energy_hartree': frame['total_energy_hartree'],
               **frame['band_edges'], 'finite_grid_integral': result['finite_grid_integral'],
               'delta_total_ev_per_cell': None, 'delta_sampled_gap_ev': None, 'dos_total_variation': None}
        if previous is not None:
            row['delta_total_ev_per_cell'] = (row['total_energy_hartree']-previous['total_energy_hartree'])*HARTREE2EV
            old_gap = previous['band_edges']['sampled_gap_ev']
            if old_gap is not None and row['sampled_gap_ev'] is not None:
                row['delta_sampled_gap_ev'] = row['sampled_gap_ev']-old_gap
            row['dos_total_variation'] = shape_change(previous, frame, settings['sigma'])
        rows.append(row)
        previous = frame
    _write_result(output, {'schema': 'wren-periodic-mesh-study-1', 'sigma_ev': settings['sigma'],
                           'energy_reference': 'same_cell_unshifted', 'meshes': rows})
    return output
