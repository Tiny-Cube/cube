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
    For projected DOS pass a (levels, channels) weight matrix. The result then
    has shape (grid, channels), sharing one Gaussian evaluation across channels.
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
    if (weights.ndim not in (1, 2) or weights.shape[0] != len(levels) or not weights.size
            or not np.isfinite(weights).all() or np.any(weights < 0)):
        raise ValueError('spectral weights must match levels and be finite/nonnegative')
    xp = _arrays(backend)
    x = xp.asarray(grid)
    out = xp.zeros((len(x), *weights.shape[1:]), dtype=x.dtype)
    # Bound both axes to avoid allocating a levels-by-grid dense matrix.
    for j in range(0, len(grid), chunk):
        here = x[j:j+chunk]
        acc = xp.zeros_like(out[j:j+chunk])
        for i in range(0, len(levels), chunk):
            delta = (here[None, :]-xp.asarray(levels[i:i+chunk])[:, None])/sigma
            kernel = xp.exp(-.5*delta*delta)
            w = xp.asarray(weights[i:i+chunk])
            if weights.ndim == 1:
                acc += xp.sum(w[:, None]*kernel, axis=0)
            else:
                acc += kernel.T @ w
        out[j:j+chunk] = acc/(math.sqrt(2*math.pi)*sigma)
    return np.asarray(out) if xp is np else xp.asnumpy(out)


def kmesh(mesh):
    """Unshifted full fractional k mesh, equal weights summing to one."""
    if len(mesh) != 3 or any(type(n) is not int or n < 1 for n in mesh):
        raise ValueError('kmesh requires three positive integers')
    axes = [np.arange(n)/n for n in mesh]
    points = np.stack(np.meshgrid(*axes, indexing='ij'), axis=-1).reshape(-1, 3)
    return points, np.full(len(points), 1/len(points))


def bloch_spectra(site_energies, pairs, kpoints, *, backend='cpu', overlap_pairs=None):
    """Eigenvalues of a normalized one-orbital-per-site tight-binding model.

    Pairs are (i,j,image[3],signed_real_coupling_eV), each undirected channel ONCE.
    Do not pass coupling magnitudes: they cannot recover interference/band DOS.
    Self-image channels produce 2V cos(k.R), never a constant site shift.
    """
    return np.asarray([values for values, _ in bloch_eigensystems(
        site_energies, pairs, kpoints, backend=backend, overlap_pairs=overlap_pairs)])


def _check_pairs(pairs, nsites):
    checked, seen = [], set()
    for i, j, image, coupling in pairs:
        if type(i) is not int or type(j) is not int or not (0 <= i < nsites and 0 <= j < nsites):
            raise ValueError('pair site index out of range')
        image = tuple(image)
        if len(image) != 3 or any(type(n) is not int for n in image):
            raise ValueError('pair image must contain three integers')
        if i == j and image == (0, 0, 0):
            raise ValueError('diagonal zero-image entries belong in the onsite diagonal')
        key = min((i, j, *image), (j, i, *(-n for n in image)))
        if key in seen:
            raise ValueError('duplicate/reversed pair would double count a coupling')
        seen.add(key)
        if isinstance(coupling, bool) or not math.isfinite(coupling):
            raise ValueError('finite signed real coupling required')
        checked.append((i, j, np.asarray(image), coupling))
    return checked


def _bloch_matrix(diagonal, pairs, k, xp):
    matrix = xp.diag(xp.asarray(diagonal, dtype=complex))
    for i, j, image, coefficient in pairs:
        value = coefficient*np.exp(2j*np.pi*np.dot(k, image))
        matrix[i, j] += value
        matrix[j, i] += value.conjugate()
    return matrix


def _overlap_power(overlap, exponent, xp):
    """Symmetric matrix power of a well-conditioned Hermitian overlap."""
    if (overlap.ndim != 2 or not overlap.size or overlap.shape[0] != overlap.shape[1]
            or not xp.isfinite(overlap).all()
            or not xp.allclose(overlap, overlap.conj().T, rtol=1e-10, atol=1e-12)):
        raise ValueError('finite Hermitian square overlap matrix required')
    values, vectors = xp.linalg.eigh(overlap)
    if float(values[0]) <= 1e-10 * max(float(values[-1]), 0.):
        raise ValueError(f'overlap must be positive definite and well conditioned (min={float(values[0]):g}, max={float(values[-1]):g})')
    return (vectors * values**exponent) @ vectors.conj().T


def element_probabilities(mol, coefficients, overlap, *, backend='cpu'):
    """Return element labels and (bands, elements) symmetric Löwdin populations."""
    xp = _arrays(backend)
    c = xp.asarray(coefficients)
    s = xp.asarray(overlap)
    if c.ndim != 2 or c.shape[0] != mol.nao_nr() or s.shape != (c.shape[0], c.shape[0]):
        raise ValueError('AO coefficient and overlap dimensions must match the basis')
    u = _overlap_power(s, .5, xp) @ c
    if not xp.allclose(u.conj().T @ u, xp.eye(c.shape[1]), atol=1e-8, rtol=0):
        raise ValueError('orbitals are not orthonormal in the supplied AO overlap metric')
    symbols = [mol.atom_pure_symbol(i) for i in range(mol.natm)]
    labels = sorted(set(symbols))
    weights = xp.abs(u)**2
    groups = xp.zeros((c.shape[1], len(labels)))
    for symbol, (_, _, start, end) in zip(symbols, mol.aoslice_by_atom(), strict=True):
        groups[:, labels.index(symbol)] += weights[start:end].sum(axis=0)
    return labels, groups if xp is np else xp.asnumpy(groups)


def bloch_eigensystems(site_energies, pairs, kpoints, *, backend='cpu', project_sites=False,
                       overlap_pairs=None):
    """Yield energies and optional symmetric Löwdin (bands, sites) probabilities."""
    sites = np.asarray(site_energies, dtype=float)
    ks = np.asarray(kpoints, dtype=float)
    if sites.ndim != 1 or not len(sites) or not np.isfinite(sites).all():
        raise ValueError('finite site energies required')
    if ks.ndim != 2 or ks.shape[1] != 3 or not len(ks) or not np.isfinite(ks).all():
        raise ValueError('finite fractional Nx3 kpoints required')
    checked = _check_pairs(pairs, len(sites))
    overlaps = None if overlap_pairs is None else _check_pairs(overlap_pairs, len(sites))
    xp = _arrays(backend)
    for k in ks:
        h = _bloch_matrix(sites, checked, k, xp)
        if overlaps:
            s = _bloch_matrix(np.ones(len(sites)), overlaps, k, xp)
            try:
                x = _overlap_power(s, -.5, xp)
            except ValueError as exc:
                raise ValueError(f'k={k.tolist()}: {exc}') from exc
            h = x @ h @ x
            h = (h + h.conj().T) * .5
        probabilities = None
        if project_sites:
            eigenvalues, vectors = xp.linalg.eigh(h)
            probabilities = (xp.abs(vectors)**2).T
            if xp is not np:
                probabilities = xp.asnumpy(probabilities)
        else:
            eigenvalues = xp.linalg.eigvalsh(h)
        yield (eigenvalues if xp is np else xp.asnumpy(eigenvalues)), probabilities


def participation_groups(energies, probabilities, tolerance_ev):
    """Participation in model-cell sites for contiguous unresolved energy groups.

    One row represents one group at a single k point. For a rank-r group the
    normalized diagonal projector is sum_band(|C_site,band|²)/r. Its second
    moment is rotation invariant; individual IPRs are returned only at rank one.
    """
    energies = np.asarray(energies, dtype=float)
    probabilities = np.asarray(probabilities, dtype=float)
    if energies.ndim != 1 or not energies.size or not np.isfinite(energies).all() or np.any(np.diff(energies) < 0):
        raise ValueError('participation analysis needs sorted finite energies')
    if (probabilities.shape != (len(energies), len(energies))
            or not np.isfinite(probabilities).all() or np.any(probabilities < 0)
            or not np.allclose(probabilities.sum(axis=0), 1., rtol=0, atol=1e-10)
            or not np.allclose(probabilities.sum(axis=1), 1., rtol=0, atol=1e-10)):
        raise ValueError('complete normalized band-by-site probabilities required')
    if isinstance(tolerance_ev, bool) or not math.isfinite(tolerance_ev) or tolerance_ev < 0:
        raise ValueError('degeneracy tolerance must be finite and nonnegative')
    boundaries = [0, *(np.flatnonzero(np.diff(energies) > tolerance_ev) + 1), len(energies)]
    for start, stop in zip(boundaries, boundaries[1:]):
        rank = int(stop-start)
        profile = probabilities[start:stop].mean(axis=0)
        profile_ipr = float(profile @ profile)
        yield {'first_band': int(start), 'last_band': int(stop-1), 'rank': rank,
               'min_ev': float(energies[start]), 'max_ev': float(energies[stop-1]),
               'ipr': profile_ipr if rank == 1 else None,
               'participation_number': 1/profile_ipr if rank == 1 else None,
               'subspace_profile_ipr': profile_ipr}


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
