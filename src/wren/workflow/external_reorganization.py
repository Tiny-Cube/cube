"""Analyze state-pair total energies transcribed from any electronic-structure code."""
import hashlib
import json
import math
from pathlib import Path

from wren.reorganization import ElectronicState, FourPointEnergies, analyze_four_point_pair


POINTS = ('E_A_RA', 'E_B_RA', 'E_B_RB', 'E_A_RB')
# Use the same Hartree/eV factor as the existing PySCF four-point adapter.
HARTREE_TO_EV = 27.21138602


def _unique_fields(items):
    record = {}
    for key, value in items:
        if key in record:
            raise ValueError(f'duplicate JSON key {key!r}')
        record[key] = value
    return record


def _fields(value, label, required, optional=()):
    if not isinstance(value, dict):
        raise ValueError(f'{label} must be a JSON object')
    missing = set(required) - value.keys()
    extra = value.keys() - set(required) - set(optional)
    if missing or extra:
        raise ValueError(f'{label}: missing fields {sorted(missing)}, unknown fields {sorted(extra)}')
    return value


def _text(value, label):
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f'{label} must be a nonempty string')
    return value.strip()


def _analyze_pair(pair, index, unit, base, hashes):
    label = f'pairs[{index}]'
    _fields(pair, label, ('id', 'states', 'energies'), ('sources', 'converged'))
    identity = _text(pair['id'], f'{label}.id')
    state_records = _fields(pair['states'], f'{identity}.states', ('a', 'b'))
    states = []
    for side in ('a', 'b'):
        record = _fields(state_records[side], f'{identity}.states.{side}',
                         ('name', 'charge', 'multiplicity', 'root', 'solver'))
        states.append(ElectronicState(**record))
    energies = _fields(pair['energies'], f'{identity}.energies', POINTS)
    try:
        values = FourPointEnergies(*(energies[key] for key in POINTS), unit=unit)
        result = analyze_four_point_pair(*states, values)
    except (OverflowError, ValueError) as exc:
        raise ValueError(f'{identity}: invalid four-point energies or states: {exc}') from exc
    factor = HARTREE_TO_EV if unit == 'hartree' else 1.0
    for field in ('relaxation', 'gaps'):
        converted = {key: value * factor for key, value in result[field].items()}
        if not all(map(math.isfinite, converted.values())):
            raise ValueError(f'{identity}: {field} overflows in eV')
        result[field + '_ev'] = converted

    result['id'] = identity
    result['declared_convergence'] = {}
    if 'converged' in pair:
        flags = _fields(pair['converged'], f'{identity}.converged', (), POINTS)
        for key, flag in flags.items():
            if flag is not True:
                raise ValueError(f'{identity}.{key}: converged must be true when supplied')
        result['declared_convergence'] = flags
    result['sources'] = {}
    if 'sources' in pair:
        sources = _fields(pair['sources'], f'{identity}.sources', (), POINTS)
        for key, filename in sources.items():
            source = (base / _text(filename, f'{identity}.sources.{key}')).resolve()
            if not source.is_file():
                raise ValueError(f'{identity}.{key}: source does not exist: {source}')
            if str(source) not in hashes:
                with source.open('rb') as stream:
                    hashes[str(source)] = hashlib.file_digest(stream, 'sha256').hexdigest()
            result['sources'][key] = str(source)
    return result


def read_external_four_points(path):
    """Read a wren-four-point-input-1 dataset and return JSON-serializable results.

    Sources are paths relative to the input JSON. Hashing a supplied source
    records its identity. Numerical values and convergence flags are the
    author's declarations and are not extracted or verified from source logs.
    """
    path = Path(path).resolve()
    raw = path.read_bytes()
    data = json.loads(raw, object_pairs_hook=_unique_fields)
    _fields(data, 'four-point dataset',
            ('schema', 'energy_kind', 'energy_unit', 'method', 'environment', 'pairs'))
    if data['schema'] != 'wren-four-point-input-1':
        raise ValueError('expected schema wren-four-point-input-1')
    if data['energy_kind'] != 'total':
        raise ValueError('energy_kind must be total: supply state total energies at each geometry')
    unit = data['energy_unit']
    if unit not in ('hartree', 'eV'):
        raise ValueError('energy_unit must be hartree or eV for all four points')
    method = _text(data['method'], 'method')
    environment = _text(data['environment'], 'environment')
    if not isinstance(data['pairs'], list) or not data['pairs']:
        raise ValueError('pairs must be a nonempty array')

    hashes = {str(path): hashlib.sha256(raw).hexdigest()}
    pairs = [_analyze_pair(pair, index, unit, path.parent, hashes)
             for index, pair in enumerate(data['pairs'])]
    if len({pair['id'] for pair in pairs}) != len(pairs):
        raise ValueError('four-point pair IDs must be unique')
    return {'schema': 'wren-external-four-point-result-1', 'status': 'completed',
            'method': '4-points', 'electronic_method': method, 'environment': environment,
            'energy_unit': unit, 'energy_kind': 'total', 'hartree_to_ev': HARTREE_TO_EV,
            'input': str(path), 'source_sha256': hashes, 'pairs': pairs}
