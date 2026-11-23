"""Quantum-model input; coordinates are read from single-frame XYZ in Å."""
import hashlib
import json
import math
from pathlib import Path

from wren.quantum.geometry import read_xyz
from .control import _integer, _real, _vector, sha256

# Conversion between chemical input units and the atomic units of the propagator.
BOHR_ANGSTROM = 0.52917721092
U_TO_ELECTRON_MASS = 1822.888486209


def read_casscf_model(filename):
    path = Path(filename).resolve()
    raw = path.read_bytes()
    data = json.loads(raw)
    required = {'schema', 'xyz', 'masses_u', 'basis', 'charge', 'active_electrons', 'active_orbitals'}
    optional = {'states', 'initial_active_mo_1based', 'threads', 'max_memory_mb', 'max_cycles',
                'min_state_overlap', 'min_gap_hartree', 'use_etfs', 'orbital_gradient_tolerance',
                'multiplicity', 'fragments'}
    if not isinstance(data, dict) or required - data.keys() or data.keys() - required - optional:
        raise ValueError('invalid fields in SA-CASSCF model')
    if data['schema'] != 'wren-sa-casscf-1' or type(data['basis']) is not str or not data['basis'].strip():
        raise ValueError('use wren-sa-casscf-1 and a PySCF orbital basis name')
    if not isinstance(data['xyz'], str) or not data['xyz'].strip() or type(data['charge']) is not int:
        raise ValueError('provide an XYZ path and an integer total charge')
    xyz_path = (path.parent / data['xyz']).resolve()
    xyz_hash = sha256(xyz_path)
    atoms = read_xyz(xyz_path)
    if sha256(xyz_path) != xyz_hash:
        raise ValueError('XYZ changed during read')
    masses = _vector(data['masses_u'], 'masses_u', positive=True)
    if len(atoms) != len(masses):
        raise ValueError('provide one isotope or atomic mass in u per XYZ atom')
    ncas = _integer(data['active_orbitals'], 'active_orbitals', 2)
    ne = _integer(data['active_electrons'], 'active_electrons')
    multiplicity = _integer(data.get('multiplicity', 1), 'multiplicity')
    spin = multiplicity-1
    if multiplicity not in (1, 2) or (ne-spin) % 2 or ne >= 2*ncas:
        raise ValueError('active space must be compatible with singlet/doublet multiplicity and contain holes')
    states = _integer(data.get('states', 2), 'states', 2)
    na, nb = (ne+spin)//2, (ne-spin)//2
    dimension = math.comb(ncas, na)*math.comb(ncas, nb)
    if na < ncas and nb:
        dimension -= math.comb(ncas, na+1)*math.comb(ncas, nb-1)
    if states > dimension:
        raise ValueError('more states than independent CAS functions for this spin')
    if 'fragments' in data:
        fragments = data['fragments']
        if (not isinstance(fragments, dict) or len(fragments) < 2 or
                any(not name.strip() or not isinstance(group, list) or not group or
                    any(type(i) is not int for i in group) for name, group in fragments.items())):
            raise ValueError('fragments must map distinct names to nonempty arrays of 1-based atom numbers')
        if sorted(i for group in fragments.values() for i in group) != list(range(1,len(atoms)+1)):
            raise ValueError('fragments must partition every atom exactly once')
    mo_indices = data.get('initial_active_mo_1based')
    if mo_indices is not None and (not isinstance(mo_indices, list) or len(mo_indices) != ncas or
            any(type(i) is not int or i < 1 for i in mo_indices) or len(set(mo_indices)) != ncas):
        raise ValueError('initial_active_mo_1based must list distinct initial RHF orbital indices')
    minimum = _real(data.get('min_state_overlap', .7), 'min_state_overlap', positive=True)
    etfs = data.get('use_etfs', True)
    if minimum > 1 or type(etfs) is not bool:
        raise ValueError('min_state_overlap must be <= 1 and use_etfs must be boolean')
    config = {**data, 'symbols': [symbol for symbol, _ in atoms], 'states': states, 'use_etfs': etfs,
              'multiplicity': multiplicity,
              'min_state_overlap': minimum, 'threads': _integer(data.get('threads', 1), 'threads'),
              'max_cycles': _integer(data.get('max_cycles', 100), 'max_cycles'),
              'max_memory_mb': _integer(data.get('max_memory_mb', 2000), 'max_memory_mb'),
              'orbital_gradient_tolerance': _real(data.get('orbital_gradient_tolerance', 1e-6),
                                                  'orbital_gradient_tolerance', positive=True),
              'min_gap_hartree': _real(data.get('min_gap_hartree', 1e-6), 'min_gap_hartree', positive=True)}
    q = tuple(x/BOHR_ANGSTROM for _, coords in atoms for x in coords)
    m = tuple(mass * U_TO_ELECTRON_MASS for mass in masses for _ in range(3))
    return config, q, m, {str(path): hashlib.sha256(raw).hexdigest(), str(xyz_path): xyz_hash}
