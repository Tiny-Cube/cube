"""Independent-trajectory FSSH with explicit forces and NAC momentum adjustment.

Atomic units throughout; real fixed orthonormal diabatic basis only. Nuclear
velocity Verlet + unitary exponential-midpoint electronic substeps. Integrated
population flux determines at most one hop at the end of a nuclear step. This
is a finite-step approximation: both nuclear/electronic convergence are needed.
"""
from dataclasses import dataclass
import numpy as np


@dataclass
class AdiabaticFrame:
    energies: np.ndarray
    vectors: np.ndarray
    gradients: np.ndarray
    nac: np.ndarray


def adiabatic_frame(provider, q, previous=None, gap_tolerance=1e-10):
    h, dh = provider.evaluate(q)
    if np.iscomplexobj(h) or np.iscomplexobj(dh):
        raise ValueError('complex Hamiltonians/SOC are not supported by this real-basis engine')
    h, dh = np.asarray(h, dtype=float), np.asarray(dh, dtype=float)
    n = h.shape[0] if h.ndim == 2 else 0
    if n < 2 or h.shape != (n, n) or dh.shape != (n, n, len(q)):
        raise ValueError('provider requires H[n,n], dH[n,n,dof], n>=2')
    if not np.isfinite(h).all() or not np.isfinite(dh).all():
        raise ValueError('non-finite Hamiltonian or derivative')
    if not np.allclose(h, h.T, atol=1e-12, rtol=0) or not np.allclose(dh, dh.swapaxes(0, 1), atol=1e-12, rtol=0):
        raise ValueError('real Hamiltonian and derivatives must be symmetric')
    energies, vectors = np.linalg.eigh(h)
    gaps = energies[None, :]-energies[:, None]
    off_diagonal = ~np.eye(n, dtype=bool)
    if np.any(abs(gaps[off_diagonal]) <= gap_tolerance):
        raise ValueError('degenerate/near-degenerate states require a local-diabatization provider')
    if previous is not None:
        if previous.shape != vectors.shape:
            raise ValueError('provider changed electronic state count')
        overlaps = np.sum(previous*vectors, axis=0)
        if np.min(abs(overlaps)) < .5:
            raise RuntimeError('adiabatic state overlap too small; reduce nuclear timestep')
        vectors *= np.where(overlaps < 0, -1., 1.)
    else:
        pivots = vectors[np.argmax(abs(vectors), axis=0), np.arange(n)]
        vectors *= np.where(pivots < 0, -1., 1.)
    transformed = np.einsum('ai,abk,bj->ijk', vectors, dh, vectors)
    nac = np.zeros_like(transformed)
    nac[off_diagonal] = transformed[off_diagonal]/gaps[off_diagonal, None]
    return AdiabaticFrame(energies, vectors, transformed[np.arange(n), np.arange(n)], nac)


def electronic_step(c, left, right, velocity, active, dt, substeps):
    """Return propagated amplitudes and FSSH probabilities (signed flux integral).

    i dc/dt = [diag(E) - i v.d]c. For a -> j the outward flux is
    2 Re[c_a* c_j (v.d_aj)]. Divide its integral by initial |c_a|².
    Negative integrated flux is discarded only after electronic subcycling.
    """
    if type(substeps) is not int or substeps < 1 or not np.isfinite(dt) or dt <= 0:
        raise ValueError('positive dt and integer substeps required')
    c = np.asarray(c, dtype=complex).copy()
    if c.shape != left.energies.shape or not 0 <= active < len(c) or not np.isfinite(c).all():
        raise ValueError('invalid amplitudes or active state')
    if not np.isclose(np.vdot(c, c).real, 1., atol=1e-10, rtol=0):
        raise ValueError('electronic amplitudes must be normalized')
    population = abs(c[active])**2
    if population < 1e-14:
        raise RuntimeError('active-state population vanished; reduce nuclear timestep')
    flux = np.zeros(len(c))
    small_dt = dt/substeps
    for index in range(substeps):
        fraction = (index+.5)/substeps
        energies = (1-fraction)*left.energies+fraction*right.energies
        nac = (1-fraction)*left.nac+fraction*right.nac
        tau = np.einsum('ijk,k->ij', nac, velocity)
        h = np.diag(energies)-1j*tau
        # One diagonalization per substep, using midpoint amplitudes for flux.
        e, u = np.linalg.eigh(h)
        spectral = u.conj().T @ c
        phases = np.exp(-.5j*e*small_dt)
        midpoint = u @ (phases*spectral)
        flux += 2*np.real(midpoint[active].conjugate()*midpoint*tau[active])*small_dt
        c = u @ (phases*phases*spectral)
    probability = np.maximum(flux/population, 0.)
    probability[active] = 0.
    if probability.sum() > 1.+1e-12:
        raise RuntimeError('FSSH probability exceeds one; reduce nuclear timestep')
    return c, probability


def rescale_momentum(momentum, masses, direction, delta_energy):
    """NAC-directed, minimal-change energy-conserving momentum jump.

    Solve A gamma² + B gamma + deltaE = 0 for p' = p + gamma d.
    A frustrated hop is rejected without velocity reversal in this version.
    """
    p, m, d = (np.asarray(v, dtype=float) for v in (momentum, masses, direction))
    if p.shape != m.shape or p.shape != d.shape or p.ndim != 1 or not np.isfinite([*p, *m, *d, delta_energy]).all() or np.any(m <= 0):
        raise ValueError('finite vectors and positive coordinate masses required')
    # Normalize the direction for scale independence and stable discriminants.
    scale = np.max(abs(d))
    if scale == 0:
        return p.copy(), False
    d = d/scale
    a = .5*np.sum(d*d/m)
    b = np.sum(p*d/m)
    discriminant = b*b-4*a*delta_energy
    if discriminant < 0:
        return p.copy(), False
    root = np.sqrt(discriminant)
    # Stable expression for the smaller-magnitude root (avoids cancellation).
    denominator = b+np.copysign(root, b)
    gamma = -2*delta_energy/denominator if denominator else 0.
    return p+gamma*d, True


def energy_decoherence(c, energies, active, kinetic, dt, constant=.1):
    """Granucci–Persico energy-based correction, Eq.17 (2007), atomic units."""
    if kinetic < 0 or constant <= 0:
        raise ValueError('EDC requires nonnegative kinetic energy and positive C')
    c = np.asarray(c, dtype=complex).copy()
    rates = abs(energies-energies[active])/(1+constant/kinetic) if kinetic > 0 else np.zeros(len(c))
    c *= np.exp(-dt*rates)
    others = np.arange(len(c)) != active
    phase = c[active]/abs(c[active]) if abs(c[active]) else 1.+0j
    c[active] = phase*np.sqrt(max(0., 1.-np.sum(abs(c[others])**2)))
    return c


@dataclass
class TrajectoryState:
    q: np.ndarray
    p: np.ndarray
    c: np.ndarray
    active: int
    frame: AdiabaticFrame
    time: float = 0.


class FSSH:
    """Provider-independent real-Hamiltonian engine, arbitrary states and DOFs."""

    def __init__(self, provider, masses, *, dt=.5, electronic_substeps=4, decoherence='none', edc_constant=.1):
        self.provider = provider
        self.masses = np.array(masses, dtype=float)
        if self.masses.ndim != 1 or not len(self.masses) or not np.isfinite(self.masses).all() or np.any(self.masses <= 0):
            raise ValueError('masses must be positive finite coordinate masses in electron masses')
        if not np.isfinite(dt) or dt <= 0 or type(electronic_substeps) is not int or electronic_substeps < 1:
            raise ValueError('dt must be positive; electronic_substeps must be a positive integer')
        if decoherence not in ('none', 'edc') or not np.isfinite(edc_constant) or edc_constant <= 0:
            raise ValueError('decoherence must be none or edc, with positive finite edc_constant')
        self.dt, self.substeps = dt, electronic_substeps
        self.decoherence, self.edc_constant = decoherence, edc_constant

    def initialize(self, coordinates, momenta, active=0, amplitudes=None):
        q, p = np.array(coordinates, dtype=float), np.array(momenta, dtype=float)
        if q.shape != self.masses.shape or p.shape != q.shape or not np.isfinite(q).all() or not np.isfinite(p).all():
            raise ValueError('finite coordinates/momenta must match coordinate masses')
        frame = adiabatic_frame(self.provider, q)
        if type(active) is not int or not 0 <= active < len(frame.energies):
            raise ValueError('active adiabatic state is zero based and out of range')
        c = np.zeros(len(frame.energies), dtype=complex)
        c[active] = 1.
        if amplitudes is not None:
            c = np.array(amplitudes, dtype=complex)
        if c.shape != frame.energies.shape or not np.isfinite(c).all() or not np.isclose(np.vdot(c, c).real, 1., atol=1e-12, rtol=0):
            raise ValueError('initial electronic amplitudes must be finite and normalized')
        if abs(c[active])**2 < 1e-14:
            raise ValueError('initial active state must have nonzero amplitude')
        return TrajectoryState(q, p, c, active, frame)

    def energy(self, state):
        return float(np.sum(state.p**2/(2*self.masses))+state.frame.energies[state.active])

    def step(self, state, rng):
        old_active = state.active
        p_half = state.p-.5*self.dt*state.frame.gradients[old_active]
        velocity = p_half/self.masses
        q = state.q+self.dt*velocity
        frame = adiabatic_frame(self.provider, q, state.frame.vectors)
        c, probability = electronic_step(state.c, state.frame, frame, velocity, old_active, self.dt, self.substeps)
        p = p_half-.5*self.dt*frame.gradients[old_active]
        target = int(np.searchsorted(np.cumsum(probability), rng.random(), side='right'))
        event, active, jump_error = 'none', old_active, 0.
        if target < len(probability):
            delta = frame.energies[target]-frame.energies[old_active]
            adjusted, accepted = rescale_momentum(p, self.masses, frame.nac[old_active, target], delta)
            if accepted:
                jump_error = float(np.sum((adjusted**2-p**2)/(2*self.masses))+delta)
                p, active, event = adjusted, target, 'accepted'
            else:
                event = 'frustrated'
        if self.decoherence == 'edc':
            c = energy_decoherence(c, frame.energies, active, float(np.sum(p*p/(2*self.masses))), self.dt, self.edc_constant)
        new = TrajectoryState(q, p, c, active, frame, state.time+self.dt)
        return new, {'event': event, 'from': old_active, 'to': target if target < len(probability) else None,
                     'probability_sum': float(probability.sum()), 'jump_energy_error': jump_error}
