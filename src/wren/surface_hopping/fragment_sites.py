"""Input for a finite, fixed-basis fragment-center Hamiltonian (standard library)."""
from dataclasses import dataclass
import hashlib
import json
import math
from pathlib import Path


@dataclass(frozen=True)
class Contact:
    i: int
    j: int
    reference: float
    coupling: float
    decay: float
    spring: float
    slope_i: float
    slope_j: float


@dataclass(frozen=True)
class FragmentParameters:
    names: tuple
    coordinates: tuple
    masses: tuple
    energies: tuple
    contacts: tuple


def _number(value, name):
    if type(value) not in (int, float) or not math.isfinite(value):
        raise ValueError(f'fragment model: {name} must be a finite real number')
    return float(value)


def _array(value, name, length):
    if not isinstance(value, list) or len(value) != length:
        raise ValueError(f'fragment model: {name} must contain {length} numbers')
    return tuple(_number(item, name) for item in value)


def read_fragment_model(filename):
    raw = Path(filename).read_bytes()
    data = json.loads(raw)
    keys = {'schema', 'basis', 'coupling_convention', 'sites', 'pairs'}
    if not isinstance(data, dict) or data.keys() != keys or data['schema'] != 'wren-fragment-sites-1':
        raise ValueError('fragment model requires wren-fragment-sites-1, basis, coupling_convention, sites and pairs')
    if data['basis'] != 'fixed_orthonormal' or data['coupling_convention'] != 'signed_real_consistent_orbital_phases':
        raise ValueError('fragment model requires a fixed orthonormal basis and signed phase-consistent couplings')
    if not isinstance(data['sites'], list) or len(data['sites']) < 2:
        raise ValueError('fragment model requires at least two sites')
    indices, centers, masses, energies = {}, [], [], []
    for site in data['sites']:
        if not isinstance(site, dict) or site.keys() != {'id', 'center_bohr', 'mass_me', 'energy_hartree'}:
            raise ValueError('fragment site requires id, center_bohr, mass_me and energy_hartree')
        name = site['id']
        if not isinstance(name, str) or not name.strip() or name in indices:
            raise ValueError('fragment site IDs must be unique nonempty strings')
        indices[name] = len(indices)
        centers.append(_array(site['center_bohr'], 'center_bohr', 3))
        mass = _number(site['mass_me'], 'mass_me')
        if mass <= 0:
            raise ValueError('fragment mass_me must be positive')
        masses.extend([mass] * 3)
        energies.append(_number(site['energy_hartree'], 'energy_hartree'))

    if not isinstance(data['pairs'], list) or not data['pairs']:
        raise ValueError('fragment model requires a nonempty contact list')
    required = {'sites', 'coupling_hartree', 'decay_per_bohr', 'spring_hartree_per_bohr2'}
    contacts, seen = [], set()
    for pair in data['pairs']:
        if not isinstance(pair, dict) or required - pair.keys() or pair.keys() - (required | {'site_slopes_hartree_per_bohr'}):
            raise ValueError('fragment pair requires sites, coupling_hartree, decay_per_bohr and spring_hartree_per_bohr2')
        ends = pair['sites']
        if not isinstance(ends, list) or len(ends) != 2 or any(not isinstance(s, str) or s not in indices for s in ends):
            raise ValueError('fragment pair sites must be two known IDs')
        i, j = (indices[name] for name in ends)
        edge = tuple(sorted((i, j)))
        if i == j or edge in seen:
            raise ValueError('duplicate or self fragment pair')
        seen.add(edge)
        r0 = math.dist(centers[i], centers[j])
        if not math.isfinite(r0) or r0 == 0:
            raise ValueError('fragment pair requires finite nonzero reference separation')
        coupling = _number(pair['coupling_hartree'], 'coupling_hartree')
        decay = _number(pair['decay_per_bohr'], 'decay_per_bohr')
        spring = _number(pair['spring_hartree_per_bohr2'], 'spring_hartree_per_bohr2')
        if min(decay, spring) < 0:
            raise ValueError('fragment decay and spring constants must be nonnegative')
        slopes = _array(pair.get('site_slopes_hartree_per_bohr', [0., 0.]), 'site_slopes_hartree_per_bohr', 2)
        contacts.append(Contact(i, j, r0, coupling, decay, spring, *slopes))
    coordinates = tuple(x for center in centers for x in center)
    parameters = FragmentParameters(tuple(indices), coordinates, tuple(masses), tuple(energies), tuple(contacts))
    return parameters, hashlib.sha256(raw).hexdigest()
