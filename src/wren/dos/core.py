"""Spectral sums; Gaussian width is sigma, not physical energetic disorder.

Gaussian-broadened DOS: sum_n w_n exp(-(E-e_n)^2/(2 sigma^2))/(sqrt(2pi)sigma).
Bloch phases use fractional k coordinates and integer lattice translations.
"""
import math
import numpy as np


def _arrays(backend):
    if backend == 'cpu':
        return np
    if backend == 'gpu':
        import cupy
        # Force initialization, so an unavailable device never falls back to CPU.
        cupy.cuda.runtime.getDeviceCount()
        return cupy
    raise ValueError('DOS backend must be cpu or gpu')


def broaden(energies, grid, sigma, *, weights=None, backend='cpu', chunk=2048):
    """Return states/eV on an ascending eV grid; integral tends to sum(weights).

    No finite-window renormalization or occupation weighting is implicit.
    Supply weight 2 for spin-degenerate spatial orbitals, 1 for each spin orbital.
    """
    levels, grid = np.asarray(energies, dtype=float), np.asarray(grid, dtype=float)
    if levels.ndim != 1 or not len(levels) or not np.isfinite(levels).all():
        raise ValueError('energies need a nonempty finite 1D eV array')
    if grid.ndim != 1 or len(grid) < 2 or not np.isfinite(grid).all() or np.any(np.diff(grid) <= 0):
        raise ValueError('grid must be finite, strictly ascending and contain at least two points')
    if isinstance(sigma, bool) or not math.isfinite(sigma) or sigma <= 0:
        raise ValueError('Gaussian sigma must be finite and positive in eV')
    if type(chunk) is not int or chunk < 1:
        raise ValueError('chunk must be a positive integer')
    weights = np.ones_like(levels) if weights is None else np.asarray(weights, dtype=float)
    if weights.shape != levels.shape or not np.isfinite(weights).all() or np.any(weights < 0):
        raise ValueError('spectral weights must match levels and be finite/nonnegative')
    xp = _arrays(backend)
    x = xp.asarray(grid)
    out = xp.zeros_like(x)
    # Bound both axes to avoid allocating a levels-by-grid dense matrix.
    for j in range(0, len(grid), chunk):
        here = x[j:j+chunk]
        acc = xp.zeros_like(here)
        for i in range(0, len(levels), chunk):
            delta = (here[None, :]-xp.asarray(levels[i:i+chunk])[:, None])/sigma
            acc += xp.sum(xp.asarray(weights[i:i+chunk])[:, None]*xp.exp(-.5*delta*delta), axis=0)
        out[j:j+chunk] = acc/(math.sqrt(2*math.pi)*sigma)
    return np.asarray(out) if xp is np else xp.asnumpy(out)


def kmesh(mesh):
    """Unshifted full fractional k mesh, equal weights summing to one."""
    if len(mesh) != 3 or any(type(n) is not int or n < 1 for n in mesh):
        raise ValueError('kmesh requires three positive integers')
    axes = [np.arange(n)/n for n in mesh]
    points = np.stack(np.meshgrid(*axes, indexing='ij'), axis=-1).reshape(-1, 3)
    return points, np.full(len(points), 1/len(points))


def bloch_spectra(site_energies, pairs, kpoints, *, backend='cpu'):
    """Eigenvalues of an orthonormal one-orbital-per-site tight-binding model.

    Pairs are (i,j,image[3],signed_real_coupling_eV), each undirected channel ONCE.
    Do not pass coupling magnitudes: they cannot recover interference/band DOS.
    Self-image channels produce 2V cos(k.R), never a constant site shift.
    """
    sites = np.asarray(site_energies, dtype=float)
    ks = np.asarray(kpoints, dtype=float)
    if sites.ndim != 1 or not len(sites) or not np.isfinite(sites).all():
        raise ValueError('finite site energies required')
    if ks.ndim != 2 or ks.shape[1] != 3 or not len(ks) or not np.isfinite(ks).all():
        raise ValueError('finite fractional Nx3 kpoints required')
    checked, seen = [], set()
    for i, j, image, coupling in pairs:
        if type(i) is not int or type(j) is not int or not (0 <= i < len(sites) and 0 <= j < len(sites)):
            raise ValueError('pair site index out of range')
        image = tuple(image)
        if len(image) != 3 or any(type(n) is not int for n in image):
            raise ValueError('pair image must contain three integers')
        if i == j and image == (0, 0, 0):
            raise ValueError('onsite energy belongs in site_energies')
        key = min((i, j, *image), (j, i, *(-n for n in image)))
        if key in seen:
            raise ValueError('duplicate/reversed pair would double count a coupling')
        seen.add(key)
        if isinstance(coupling, bool) or not math.isfinite(coupling):
            raise ValueError('finite signed real coupling required')
        checked.append((i, j, np.asarray(image), coupling))
    xp = _arrays(backend)
    values = []
    for k in ks:
        h = xp.diag(xp.asarray(sites, dtype=complex))
        for i, j, image, coupling in checked:
            value = coupling*np.exp(2j*np.pi*np.dot(k, image))
            h[i, j] += value
            h[j, i] += value.conjugate()
        eigenvalues = xp.linalg.eigvalsh(h)
        values.append(np.asarray(eigenvalues) if xp is np else xp.asnumpy(eigenvalues))
    return np.asarray(values)


def disorder_statistics(levels):
    """Exact population-variance decomposition of a balanced time-by-site table.

    static = var_i(mean_t e_ti); dynamic = mean_i(var_t e_ti).
    This finite-sample descriptor is not an independence test or a confidence interval.
    """
    values = np.asarray(levels, dtype=float)
    if values.ndim != 2 or not values.size or not np.isfinite(values).all():
        raise ValueError('balanced finite frames-by-sites energy table required')
    static = float(np.var(np.mean(values, axis=0)))
    dynamic = float(np.mean(np.var(values, axis=0)))
    return {'mean_ev': float(np.mean(values)), 'variance_total_ev2': float(np.var(values)),
            'variance_static_ev2': static, 'variance_dynamic_ev2': dynamic,
            'frames': len(values), 'sites': values.shape[1],
            'statistical_uncertainty_estimated': False}
