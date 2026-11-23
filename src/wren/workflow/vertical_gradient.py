"""Harmonic vertical-gradient estimate in an explicitly supplied mode subspace.

At the declared minimum of one state, use its positive normal curvatures as
an approximation to the curvature of the other state. In mass-weighted
coordinates, a target-state gradient g gives dQ = -K⁺g and λ = gᵀK⁺g/2.
This is one target-surface relaxation. No geometry optimization, Hessian
calculation, cutoff, rigid-body projection or mode matching is performed.
"""
from dataclasses import asdict
import hashlib
import json
import math
from pathlib import Path

from wren.reorganization import ElectronicState
from .external_reorganization import _fields, _text, _unique_fields

# c, h, e are exact SI constants; the atomic mass constant is CODATA 2022.
C_M_S = 299792458.0
H_J_S = 6.62607015e-34
EV_J = 1.602176634e-19
U_KG = 1.66053906892e-27
GRADIENT = 'gradient_ev_per_angstrom_sqrt_u'


def _number(value, label):
    if type(value) not in (int, float) or not math.isfinite(value):
        raise ValueError(f'{label} must be finite and numeric')
    return float(value)


def _cartesian_projection(data, modes):
    import numpy as np

    _fields(data, 'cartesian', ('atom_ids', 'masses_u', 'gradient_ev_per_angstrom'))
    atom_ids = data['atom_ids']
    if (not isinstance(atom_ids, list) or not atom_ids or
            any(not isinstance(a, str) or not a.strip() for a in atom_ids) or len(set(atom_ids)) != len(atom_ids)):
        raise ValueError('provide distinct, ordered string atom_ids for the mobile Cartesian input')
    n = len(atom_ids)

    def array(value, shape, label):
        arr = np.asarray(value, dtype=object)
        if arr.shape != shape or any(type(x) not in (int, float) for x in arr.flat):
            raise ValueError(f'{label} must contain numbers with shape {shape}')
        arr = arr.astype(float)
        if not np.isfinite(arr).all():
            raise ValueError(f'{label} must be finite')
        return arr

    masses = array(data['masses_u'], (n,), 'masses_u')
    if np.any(masses <= 0):
        raise ValueError('atomic masses in u must be positive')
    gradient = array(data['gradient_ev_per_angstrom'], (n, 3), 'Cartesian gradient')
    basis = array([m['mass_weighted_direction'] for m in modes], (len(modes), n, 3),
                  'mass_weighted_direction').reshape(len(modes), -1)
    if len(modes) > 3*n or not np.allclose(basis @ basis.T, np.eye(len(modes)), atol=1e-6, rtol=0):
        raise ValueError('mass-weighted normal vectors must be orthonormal in the supplied Cartesian subspace')
    g = (gradient / np.sqrt(masses[:, None])).ravel()
    projected = basis @ g
    residual = float(np.linalg.norm(g - basis.T @ projected))
    if not math.isfinite(residual):
        raise ValueError('unrepresented gradient norm overflowed')
    return projected.tolist(), {'input': 'cartesian', 'atom_ids': atom_ids, 'cartesian_dimension': 3*n,
                                'unrepresented_gradient_norm_ev_per_angstrom_sqrt_u': residual}


def read_vertical_gradient(path):
    """Read wren-vertical-gradient-1 JSON. Mode gradients are dE/dQ, in eV/(Å√u)."""
    path = Path(path).resolve()
    raw = path.read_bytes()
    data = json.loads(raw, object_pairs_hook=_unique_fields)
    _fields(data, 'VG input', ('schema', 'method', 'environment', 'states', 'reference',
                               'reference_stationary', 'mode_subspace', 'modes'), ('cartesian',))
    if data['schema'] != 'wren-vertical-gradient-1':
        raise ValueError('expected wren-vertical-gradient-1')
    method = _text(data['method'], 'method')
    environment = _text(data['environment'], 'environment')
    subspace = _text(data['mode_subspace'], 'mode_subspace')
    if data['reference_stationary'] is not True or data['reference'] not in ('a', 'b'):
        raise ValueError('declare a stationary reference geometry with reference "a" or "b"')
    _fields(data['states'], 'states', ('a', 'b'))
    states = {}
    for side, value in data['states'].items():
        _fields(value, f'state {side}', ('name', 'charge', 'multiplicity', 'root', 'solver'))
        states[side] = asdict(ElectronicState(**value))
    if states['a']['name'] == states['b']['name']:
        raise ValueError('states must have distinct names')
    modes = data['modes']
    if not isinstance(modes, list) or not modes:
        raise ValueError('provide at least one selected positive-frequency mode')
    components = 'mass_weighted_direction' if 'cartesian' in data else GRADIENT
    for mode in modes:
        _fields(mode, 'mode', ('id', 'wavenumber_cm1', components))
    ids = [_text(mode['id'], 'mode id') for mode in modes]
    if len(set(ids)) != len(ids):
        raise ValueError('mode IDs must be distinct')
    if 'cartesian' in data:
        gradients, projection = _cartesian_projection(data['cartesian'], modes)
    else:
        gradients = [_number(mode[GRADIENT], GRADIENT) for mode in modes]
        projection = {'input': 'projected', 'unrepresented_gradient_norm_ev_per_angstrom_sqrt_u': None}

    contributions = []
    for mode_id, mode, g in zip(ids, modes, gradients):
        wavenumber = _number(mode['wavenumber_cm1'], 'wavenumber_cm1')
        if wavenumber <= 0:
            raise ValueError('select positive wavenumbers in cm^-1; resolve imaginary or zero modes explicitly')
        frequency_hz = C_M_S * 100 * wavenumber
        angular_frequency = 2 * math.pi * frequency_hz
        curvature = angular_frequency * angular_frequency * U_KG * 1e-20 / EV_J
        quantum = H_J_S * frequency_hz / EV_J
        if not math.isfinite(curvature) or min(curvature, quantum) <= 0:
            raise ValueError('wavenumber outside numerical range')
        displacement = -g / curvature
        relaxation = -.5 * g * displacement
        hr = relaxation / quantum
        if not all(map(math.isfinite, (displacement, relaxation, hr))):
            raise ValueError('mode relaxation exceeds numerical range')
        contributions.append({'id': mode_id, 'wavenumber_cm1': wavenumber, GRADIENT: g,
                              'curvature_ev_per_angstrom2_u': curvature, 'quantum_ev': quantum,
                              'displacement_angstrom_sqrt_u': displacement, 'relaxation_ev': relaxation,
                              'huang_rhys': hr})
    total = math.fsum(mode['relaxation_ev'] for mode in contributions)
    return {'schema': 'wren-vertical-gradient-result-1', 'status': 'completed', 'method': 'vertical-gradient',
            'input': str(path), 'source_sha256': {str(path): hashlib.sha256(raw).hexdigest()},
            'electronic_method': method, 'environment': environment, 'states': states,
            'reference': data['reference'], 'target': 'b' if data['reference'] == 'a' else 'a',
            'scope': 'one_target_surface_relaxation_in_supplied_reference_mode_subspace',
            'mode_subspace': subspace, 'projection': projection, 'modes': contributions,
            'relaxation_ev': total, 'diagnostics': {'reference_stationary_declared': True,
                'reference_minimum_verified': False, 'target_reference_curvatures_equal_verified': False}}
