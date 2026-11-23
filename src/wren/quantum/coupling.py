"""Electronic-coupling algebra, in the energy units supplied by the caller.

References are placed beside each implementation. Only short source excerpts
are quoted; formulas are independently implemented, not copied upstream code.
Real, orthonormal adiabatic states or explicitly supplied overlap metrics only.
"""
import numpy as np


class ApplicabilityError(ValueError):
    """A requested approximation is not justified by the supplied states."""


def array_module(value):
    if hasattr(value, '__cuda_array_interface__'):
        import cupy
        return cupy
    return np


def host(value):
    return value.get() if hasattr(value, 'get') else np.asarray(value)


def _real(value, xp):
    if xp.iscomplexobj(value):
        raise ValueError('complex orbitals are not supported by this real-state implementation')
    value = xp.asarray(value, dtype=xp.float64)
    if not bool(xp.all(xp.isfinite(value))):
        raise ValueError('matrix entries must be finite')
    return value


def lowdin(hamiltonian, overlap, relative_cutoff=1e-8):
    """Return S^(-1/2) H S^(-1/2), S^(-1/2), and overlap eigenvalues.

    Baumeier, Kirkpatrick & Andrienko, PCCP 12, 11103 (2010), Sec. II,
    Eq. (7), DOI:10.1039/C002337J. Original abstract excerpt:
    "forming the expectation value of a dimer Fock operator"
    Author record: https://research.tue.nl/en/publications/
    density-functional-based-determination-of-intermolecular-charge-t/
    The two-state off-diagonal consequence is independently checked against
    V=(J-Sab*(ea+eb)/2)/(1-Sab**2), also in the FHI-aims manual Eq. (3.259).
    """
    xp = array_module(hamiltonian)
    h, s = _real(hamiltonian, xp), _real(overlap, xp)
    if h.ndim != 2 or h.shape[0] != h.shape[1] or h.shape != s.shape or not h.size:
        raise ValueError('H and S must be nonempty square matrices of equal shape')
    if not 0 < relative_cutoff < 1:
        raise ValueError('overlap cutoff must be between zero and one')
    if not bool(xp.allclose(h, h.T, atol=1e-11, rtol=1e-10)) or not bool(xp.allclose(s, s.T, atol=1e-11, rtol=1e-10)):
        raise ValueError('H and S must be symmetric')
    eigenvalues, vectors = xp.linalg.eigh(s)
    if float(eigenvalues[-1]) <= 0 or float(eigenvalues[0]) <= relative_cutoff*float(eigenvalues[-1]):
        raise ApplicabilityError('fragment overlap is singular or ill-conditioned; no inverse-square-root truncation is applied')
    inverse_sqrt = (vectors / xp.sqrt(eigenvalues)[None, :]) @ vectors.T
    return inverse_sqrt @ h @ inverse_sqrt, inverse_sqrt, eigenvalues


def energy_splitting(energies, site_energies, *, equivalent=False, tolerance=1e-4):
    """Magnitude |V|=|E2-E1|/2 ONLY in the equivalent-site two-state model.

    Cave & Newton, JCP 106, 9213 (1997), Sec. I,
    DOI:10.1063/1.474023; author copy https://scholarship.claremont.edu/hmc_fac_pub/8/.
    The original treatment restricts half splitting to symmetry-equivalent
    donor/acceptor groups and an adequate two-state model. Valeev et al.,
    JACS 128, 9882 (2006), DOI:10.1021/ja061827h, warn in the abstract:
    "neglect of electronic polarization leads to qualitatively incorrect values"
    This function requires explicit equivalence AND a site-energy check;
    the QC workflow additionally checks the selected dimer MO subspace.
    """
    e, sites = _real(energies, np), _real(site_energies, np)
    if e.shape != (2,) or sites.shape != (2,):
        raise ValueError('two adiabatic and two orthogonalized site energies are required')
    if not np.isfinite(tolerance) or tolerance < 0:
        raise ValueError('site-energy tolerance must be finite and nonnegative')
    if not equivalent:
        raise ApplicabilityError('ESID requires an explicit symmetry-equivalent-fragment assumption')
    if abs(sites[1]-sites[0]) > tolerance:
        raise ApplicabilityError('ESID rejected: projected site energies are not equivalent')
    return float(abs(e[1]-e[0])/2)


def gmh(energies, dipoles, *, axis=None, dipole_tolerance=1e-8):
    """Two-state GMH magnitude, for orthonormal adiabatic input states.

    Cave & Newton, Chem. Phys. Lett. 249, 15 (1996),
    DOI:10.1016/0009-2614(95)01310-5. Original abstract excerpt:
    "quantities obtained solely in terms of the adiabatic states"
    Author record: https://scholarship.claremont.edu/hmc_fac_pub/816/.
    V=|deltaE*mu12|/sqrt((mu22-mu11)**2+4*mu12**2).
    Project vector dipoles on the difference-dipole direction. When that
    vector vanishes by symmetry, use the transition-dipole direction instead.
    Applying this algebra to KS orbitals is an ORBITAL approximation; it does
    not generate many-electron adiabatic states or their transition densities.
    """
    e, mu = _real(energies, np), _real(dipoles, np)
    if e.shape != (2,) or mu.shape != (3, 2, 2):
        raise ValueError('energies (2,) and dipoles (3,2,2) are required')
    if not np.allclose(mu, mu.transpose(0, 2, 1), atol=1e-11, rtol=1e-10):
        raise ValueError('dipole matrices must be symmetric')
    if not np.isfinite(dipole_tolerance) or dipole_tolerance <= 0:
        raise ValueError('dipole tolerance must be positive and finite')
    if axis is None:
        axis = mu[:, 1, 1]-mu[:, 0, 0]
        if np.linalg.norm(axis) <= dipole_tolerance:
            axis = mu[:, 0, 1]
    direction = _real(axis, np)
    if direction.shape != (3,) or np.linalg.norm(direction) <= dipole_tolerance:
        raise ApplicabilityError('GMH undefined: no nonzero dipole direction')
    direction = direction/np.linalg.norm(direction)
    projected = np.einsum('x,xij->ij', direction, mu)
    delta, transition = projected[1, 1]-projected[0, 0], projected[0, 1]
    denominator = np.hypot(delta, 2*transition)
    if denominator <= dipole_tolerance:
        raise ApplicabilityError('GMH undefined: zero dipole contrast and transition moment')
    return {'magnitude': float(abs((e[1]-e[0])*transition)/denominator),
            'axis': direction.tolist(), 'projected_dipoles': projected.tolist()}
