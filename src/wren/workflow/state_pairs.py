"""Explicit electronic states used by workflow controls (standard library only)."""
from dataclasses import asdict
import tomllib

from wren.reorganization import ElectronicState


def parse_state(text):
    """Accept common singlet shortcuts or a TOML inline state definition."""
    if text.upper() in ('S0', 'S1'):
        excited = text.upper() == 'S1'
        return ElectronicState(text.upper(), 0, 1, int(excited), 'tda' if excited else 'scf')
    try:
        data = tomllib.loads('state = ' + text)['state']
        if not isinstance(data, dict) or not {'name', 'charge', 'multiplicity', 'solver', 'root'} == data.keys():
            raise ValueError
        return ElectronicState(**data)
    except (ValueError, TypeError, KeyError) as exc:
        raise ValueError('state must be S0, S1, or an inline table with name, charge, multiplicity, solver and root') from exc


def explicit_pair(settings):
    if not {'state-a', 'state-b'} <= settings.keys():
        raise ValueError('both state-a and state-b are required')
    a, b = (parse_state(settings[key]) for key in ('state-a', 'state-b'))
    if a.name == b.name:
        raise ValueError('state names must differ')
    return {'state-a': asdict(a), 'state-b': asdict(b)}


def quantum_pair_adapter(pair):
    a, b = (ElectronicState(**pair[key]) for key in ('state-a', 'state-b'))
    if all(s.solver == 'scf' and s.root == 0 for s in (a, b)):
        if (a.charge, a.multiplicity) == (b.charge, b.multiplicity):
            raise ValueError('SCF alone cannot distinguish these states; supply external state-specific energies')
        return 'scf'
    if (a.solver, a.root, b.solver, b.root) == ('scf', 0, 'tda', 1) and (
            a.multiplicity == b.multiplicity == 1 and a.charge == b.charge):
        return 'tda'
    raise ValueError('automatic four-point supports two SCF sectors or singlet SCF→TDA root 1; '
                     'use energy-source=external for other defined states')
