"""Real, orthonormal, coordinate-independent diabatic Hamiltonian providers.

All quantities use atomic units. evaluate(q) returns H[n,n], dH[n,n,dof].
An atomistic provider must include all energy and force terms in this contract;
a moving/nonorthogonal fragment basis needs additional basis-derivative terms.
"""
from typing import Protocol
import math
import numpy as np


class HamiltonianProvider(Protocol):
    def evaluate(self, coordinates: np.ndarray) -> tuple[np.ndarray, np.ndarray]: ...


class FragmentSitesModel:
    """Signed radial transfer, site-energy shifts and common pair springs.

    Each finite-cluster site has one fixed orthonormal electronic basis state
    and three Cartesian nuclear coordinates. All quantities are in atomic units.
    """

    def __init__(self, parameters):
        self.parameters = parameters

    def evaluate(self, coordinates):
        n = len(self.parameters.energies)
        q = np.asarray(coordinates)
        if np.iscomplexobj(q) or q.shape != (3*n,) or not np.isfinite(q).all():
            raise ValueError('fragment model requires three finite real coordinates per center')
        centers = np.asarray(q, dtype=float).reshape(n, 3)
        h = np.diag(self.parameters.energies)
        dh = np.zeros((n, n, n, 3))
        elastic_energy, elastic_gradient = 0., np.zeros((n, 3))
        for pair in self.parameters.contacts:
            i, j = pair.i, pair.j
            displacement = centers[i] - centers[j]
            r = math.hypot(*displacement)
            if not math.isfinite(r) or r == 0:
                raise ValueError('fragment contact has undefined radial derivative at this separation')
            direction = displacement / r
            dr = r - pair.reference
            try:
                transfer = pair.coupling * math.exp(-pair.decay * dr) if pair.coupling else 0.
            except OverflowError as exc:
                raise ValueError('fragment transfer overflows: check separation, decay and units') from exc
            h[i, j] = h[j, i] = transfer
            h[i, i] += pair.slope_i * dr
            h[j, j] += pair.slope_j * dr
            radial_gradient = -pair.decay * transfer * direction
            dh[i, j, i] = dh[j, i, i] = radial_gradient
            dh[i, j, j] = dh[j, i, j] = -radial_gradient
            for site, slope in ((i, pair.slope_i), (j, pair.slope_j)):
                dh[site, site, i] += slope * direction
                dh[site, site, j] -= slope * direction
            elastic_energy += .5 * pair.spring * dr * dr
            spring_gradient = pair.spring * dr * direction
            elastic_gradient[i] += spring_gradient
            elastic_gradient[j] -= spring_gradient
        diagonal = np.arange(n)
        h[diagonal, diagonal] += elastic_energy
        dh[diagonal, diagonal] += elastic_gradient
        if not np.isfinite(h).all() or not np.isfinite(dh).all():
            raise ValueError('nonfinite fragment energy/gradient: check geometry, parameters and units')
        return h, dh.reshape(n, n, 3*n)


class LinearVibronicModel:
    """H(q)=H0 + sum_k G_k q_k + 1/2 sum_k K_k q_k² I.

    A model interface for site/exciton Hamiltonians in a fixed orthonormal basis.
    Kinetic energy uses the coordinate masses supplied to FSSH. It is not an
    atomistic force field or an automatic fit to a static coupling CSV.
    """

    def __init__(self, h0, linear, harmonic):
        if any(np.iscomplexobj(v) for v in (h0, linear, harmonic)):
            raise ValueError('LVC parameters must be real')
        self.h0, self.linear, self.harmonic = (np.array(v, dtype=float) for v in (h0, linear, harmonic))
        n = len(self.h0)
        if self.h0.shape != (n, n) or n < 2 or self.harmonic.ndim != 1 or not len(self.harmonic) or self.linear.shape != (n, n, len(self.harmonic)):
            raise ValueError('LVC requires H0[n,n], linear[n,n,dof], harmonic[dof]')
        if not all(np.isfinite(v).all() for v in (self.h0, self.linear, self.harmonic)) or np.any(self.harmonic < 0):
            raise ValueError('finite LVC parameters and nonnegative harmonic constants required')
        if not np.allclose(self.h0, self.h0.T, atol=1e-12, rtol=0) or not np.allclose(self.linear, self.linear.swapaxes(0, 1), atol=1e-12, rtol=0):
            raise ValueError('LVC matrices must be symmetric')

    def evaluate(self, coordinates):
        q = np.asarray(coordinates, dtype=float)
        if q.shape != self.harmonic.shape or not np.isfinite(q).all():
            raise ValueError('LVC coordinate dimension/values invalid')
        identity = np.eye(len(self.h0))
        h = self.h0+np.einsum('ijk,k->ij', self.linear, q)+.5*np.sum(self.harmonic*q*q)*identity
        dh = self.linear+identity[..., None]*(self.harmonic*q)
        return h, dh


class TullyModel:
    """Models 1–3 and original parameters from Tully, JCP 93, 1061 (1990)."""

    def __init__(self, number=1):
        if number not in (1, 2, 3):
            raise ValueError('Tully model number must be 1, 2 or 3')
        self.number = number

    def evaluate(self, coordinates):
        q = np.asarray(coordinates, dtype=float)
        if q.shape != (1,) or not np.isfinite(q).all():
            raise ValueError('Tully models require one finite coordinate in Bohr')
        x = float(q[0])
        if self.number == 1:
            a, b, c, d = .01, 1.6, .005, 1.
            v = np.sign(x)*a*(-np.expm1(-b*abs(x)))
            dv = a*b*np.exp(-b*abs(x))
            w = c*np.exp(-d*x*x)
            dw = -2*d*x*w
            h = [[v, w], [w, -v]]
            dh = [[dv, dw], [dw, -dv]]
        elif self.number == 2:
            a, b, c, d, e = .10, .28, .015, .06, .05
            v = -a*np.exp(-b*x*x)+e
            dv = 2*a*b*x*np.exp(-b*x*x)
            w = c*np.exp(-d*x*x)
            dw = -2*d*x*w
            h = [[0., w], [w, v]]
            dh = [[0., dw], [dw, dv]]
        else:
            a, b, c = .0006, .10, .90
            exponential = np.exp(-c*abs(x))
            w = b*exponential if x < 0 else b*(2-exponential)
            dw = b*c*exponential
            h = [[a, w], [w, -a]]
            dh = [[0., dw], [dw, 0.]]
        return np.array(h), np.array(dh)[..., None]
