"""Jackson kernel polynomial DOS for orthonormal sparse Bloch Hamiltonians.

Recurrence, stochastic trace and damping follow Weiße et al., Rev. Mod. Phys.
78, 275 (2006), doi:10.1103/RevModPhys.78.275. Output bins integrate the
Chebyshev series analytically, including bins intersecting the support edges.
"""
import csv

import numpy as np
from scipy import sparse

from .core import _arrays, _check_pairs, kmesh


def _scaling(sites, left, right, couplings, bounds):
    # Triangle inequality before inserting Bloch phases bounds every k at once.
    radii = np.zeros(len(sites))
    np.add.at(radii, left, np.abs(couplings))
    np.add.at(radii, right, np.abs(couplings))
    low, high = float(np.min(sites-radii)), float(np.max(sites+radii))
    if not np.isfinite([low, high]).all():
        raise ValueError('nonfinite Gershgorin energy bounds')
    lower, upper = (low, high) if bounds is None else bounds
    if lower > low or upper < high:
        raise ValueError(f'kpm-bounds-ev must enclose the conservative bounds [{low:g}, {high:g}]')
    if lower == upper:
        raise ValueError('flat spectrum: supply kpm-bounds-ev to set a finite kernel energy scale')
    center = .5*lower + .5*upper
    scale = 1.01*(.5*upper - .5*lower)
    if not np.isfinite([center-scale, center+scale]).all() or scale <= 0:
        raise ValueError('KPM scaling exceeds floating-point energy range')
    return center, scale, [low, high]


def _probe_moments(matrix, order, settings, xp):
    """Return probe-by-order traces, or one summed deterministic trace."""
    n = matrix.shape[0]
    exact = settings['kpm-trace'] == 'exact'
    count = n if exact else settings['kpm-vectors']
    rng = np.random.default_rng(settings.get('kpm-seed', 0))
    sums = np.zeros((1 if exact else count, order))
    for start in range(0, count, 8):
        batch = min(8, count-start)
        if exact:
            probe = xp.zeros((n, batch), dtype=matrix.dtype)
            index = xp.arange(batch)
            probe[start+index, index] = 1.
        else:
            signs = 2*rng.integers(0, 2, (batch, n), dtype=np.int8)-1
            probe = xp.asarray(signs.T, dtype=matrix.dtype, order='C')
        previous, current = probe, matrix @ probe
        moments = xp.empty((batch, order), dtype=xp.float64)
        moments[:, 0] = 1 if exact else n
        # All probes are real; only the real part contributes to their expectation.
        for m in range(1, order):
            if exact:
                moments[:, m] = current[start+index, index].real
            else:
                moments[:, m] = xp.sum(probe.real * current.real, axis=0)
            if m+1 < order:
                previous, current = current, 2*(matrix @ current)-previous
        values = moments if xp is np else xp.asnumpy(moments)
        if exact:
            sums[0] += values.sum(axis=0)
        else:
            sums[start:start+batch] = values
    return sums


def bin_densities(moments, edges, center, scale):
    """Yield slices of bin means and probe standard errors in states/eV."""
    order = moments.shape[1]
    n = np.arange(1, order)
    phi = np.pi/(order+1)
    jackson = ((order-n+1)*np.cos(n*phi) + np.sin(n*phi)/np.tan(phi))/(order+1)
    x = np.clip((edges-center)/scale, -1., 1.)
    x[edges <= center-scale] = -1.
    x[edges >= center+scale] = 1.
    theta = np.arccos(x)
    # Bound temporary kernel and per-probe bin arrays independently of grid size.
    block = max(1, 262144//max(order, len(moments)))
    for start in range(0, len(edges)-1, block):
        stop = min(len(edges)-1, start+block)
        left, right = theta[start:stop], theta[start+1:stop+1]
        delta = left-right
        mid = .5*(left+right)
        # Integral of T_n(x)/sqrt(1-x²) is (sin(n*left)-sin(n*right))/n.
        kernel = (4*jackson/n)*np.cos(mid[:, None]*n)*np.sin(.5*delta[:, None]*n)
        values = moments[:, :1]*delta + moments[:, 1:] @ kernel.T
        values /= np.pi*np.diff(edges[start:stop+1])
        sem = values.std(axis=0, ddof=1)/np.sqrt(len(values)) if len(values)>1 else np.zeros(len(delta))
        yield slice(start, stop), values.mean(axis=0), sem


def calculate(directory, network, pairs, settings):
    sites = np.asarray(network['site_energies_ev'], dtype=float)
    if sites.ndim != 1 or sites.size == 0 or not np.isfinite(sites).all():
        raise ValueError('KPM requires finite site energies')
    checked = _check_pairs(pairs, len(sites))
    left = np.asarray([p[0] for p in checked], dtype=np.intp)
    right = np.asarray([p[1] for p in checked], dtype=np.intp)
    images = np.asarray([p[2] for p in checked], dtype=float).reshape(-1, 3)
    couplings = np.asarray([p[3] for p in checked], dtype=float)
    center, scale, bounds = _scaling(sites, left, right, couplings, settings.get('kpm-bounds-ev'))
    order = settings['kpm-moments']
    bins = settings.get('points', max(512, order*4))
    lower, upper = settings.get('emin', center-scale), settings.get('emax', center+scale)
    if lower >= upper or bins > 1000000:
        raise ValueError('KPM output needs increasing energy bounds and at most 1000000 bins')
    edges = np.linspace(lower, upper, bins+1)
    if not np.all(np.diff(edges)>0):
        raise ValueError('KPM bin widths are too small for the energy scale')

    xp = _arrays(settings['dos-backend'])
    if xp is not np:
        from cupyx.scipy import sparse as gpu_sparse
    points, kw = kmesh(settings['kmesh'])
    indices = np.arange(len(sites))
    row, col = np.r_[indices, left, right], np.r_[indices, right, left]
    total = None
    max_nnz = 0
    for k, weight in zip(points, kw, strict=True):
        hopping = couplings * np.exp(2j*np.pi*(images @ k))
        if not hopping.imag.any():
            hopping = hopping.real
        values = np.r_[(sites-center)/scale, hopping/scale, hopping.conj()/scale]
        matrix = sparse.csr_matrix((values, (row, col)), shape=(len(sites), len(sites)))
        matrix.eliminate_zeros()
        max_nnz = max(max_nnz, matrix.nnz)
        if xp is not np:
            matrix = gpu_sparse.csr_matrix(matrix)
        # Reuse the same probes across k, then estimate sampling variance across probes.
        traces = _probe_moments(matrix, order, settings, xp) * weight
        if total is None:
            total = traces
        else:
            total += traces
    total *= settings['spin-degeneracy']
    if not np.isfinite(total).all():
        raise ValueError('nonfinite Chebyshev trace')
    mean = total.mean(axis=0)
    se = total.std(axis=0, ddof=1)/np.sqrt(len(total)) if len(total)>1 else np.zeros(order)
    with (directory/'moments.csv').open('w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['order', 'trace_moment', 'trace_standard_error'])
        writer.writerows(zip(range(order), mean, se, strict=True))

    integral = 0.
    with (directory/'dos.csv').open('w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['energy_ev', 'dos', 'trace_standard_error', 'left_ev', 'right_ev'])
        for sl, density, stderr in bin_densities(total, edges, center, scale):
            a, b = edges[:-1][sl], edges[1:][sl]
            integral += float(np.dot(density, b-a))
            writer.writerows(zip(a+(b-a)/2, density, stderr, a, b, strict=True))
    expected = len(sites)*settings['spin-degeneracy']
    return {'scope': 'orthonormal_effective_Hamiltonian', 'normalization': 'states_per_model_cell_per_eV',
            'solver': 'kpm', 'broadening': {'kernel': 'Jackson', 'moments': order},
            'frames': [{'moments': 'moments.csv', 'kpoints_fractional': points.tolist(), 'kpoint_weights': kw.tolist()}],
            'frame_count': 1, 'expected_integral': expected, 'finite_grid_integral': integral,
            'integral_relative_error': abs(integral-expected)/expected,
            'kpm': {'gershgorin_bounds_ev': bounds, 'center_ev': center, 'scale_ev': scale,
                    'trace': settings['kpm-trace'], 'random_vectors': settings.get('kpm-vectors', 0),
                    'seed': settings.get('kpm-seed'), 'max_scaled_csr_nnz': max_nnz,
                    'grid': 'bin_averages', 'bin_count': bins, 'trace_error_estimated': len(total)>1},
            'grid_spacing_ev': float(edges[1]-edges[0])}
