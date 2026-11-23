"""State-pair four-point algebra, independent of quantum engines and NumPy.

Inputs are total energies of the same two state surfaces at their two
reference geometries, in one common unit and a consistently defined environment.
This module does not generate geometries, identify states or verify minima.
"""
from dataclasses import asdict, dataclass
import math
from numbers import Real


def _energy(value):
    if isinstance(value, bool) or not isinstance(value, Real) or not math.isfinite(value):
        raise ValueError('four-point energies must be finite real numbers')
    return float(value)


def four_point_relaxations(e_a_ra, e_b_ra, e_b_rb, e_a_rb):
    """Signed deformations on surfaces A and B, plus their sum.

    All values use caller units. The sum is not automatically the Marcus
    reorganization energy of a single transition; that convention belongs
    to the physical process consuming these surface deformations.
    """
    aa, ba, bb, ab = map(_energy, (e_a_ra, e_b_ra, e_b_rb, e_a_rb))
    on_a, on_b = ab-aa, ba-bb
    result = {'on_a': on_a, 'on_b': on_b, 'sum': on_a+on_b}
    if not all(math.isfinite(value) for value in result.values()):
        raise ValueError('four-point energy differences overflowed')
    return result


@dataclass(frozen=True)
class ElectronicState:
    """An explicit state identity, not a promise of solver support.

    Root indices refer to the specified solver's state namespace. The same
    charge/multiplicity/root can describe distinct localized states, so the
    name must retain their identity. Electron-count parity needs the QM atoms
    and must be checked by an electronic-structure adapter.
    """
    name: str
    charge: int
    multiplicity: int
    root: int = 0
    solver: str = 'external'

    def __post_init__(self):
        for field in ('name', 'solver'):
            if not isinstance(getattr(self, field), str) or not getattr(self, field).strip():
                raise ValueError(f'state {field} must be a nonempty string')
        for field, minimum in (('charge', None), ('multiplicity', 1), ('root', 0)):
            value = getattr(self, field)
            if type(value) is not int or (minimum is not None and value < minimum):
                raise ValueError(f'invalid state {field}')


@dataclass(frozen=True)
class FourPointEnergies:
    e_a_ra: float
    e_b_ra: float
    e_b_rb: float
    e_a_rb: float
    unit: str

    def __post_init__(self):
        if self.unit not in ('hartree', 'eV'):
            raise ValueError('energy unit must be explicitly hartree or eV')
        for field in ('e_a_ra', 'e_b_ra', 'e_b_rb', 'e_a_rb'):
            object.__setattr__(self, field, _energy(getattr(self, field)))


def analyze_four_point_pair(state_a, state_b, energies):
    """Analyze externally supplied energies without running a QC calculation.

    No state-dependent energy alignment, unit conversion, absolute value,
    averaging of relaxation terms or implicit factor of two is applied.
    """
    if not isinstance(state_a, ElectronicState) or not isinstance(state_b, ElectronicState):
        raise TypeError('both states must be ElectronicState objects')
    if state_a.name == state_b.name:
        raise ValueError('state names must be distinct')
    if not isinstance(energies, FourPointEnergies):
        raise TypeError('energies must be FourPointEnergies')
    aa, ba, bb, ab = (energies.e_a_ra, energies.e_b_ra, energies.e_b_rb, energies.e_a_rb)
    terms = four_point_relaxations(aa, ba, bb, ab)
    gaps = {'vertical_at_ra': ba-aa, 'vertical_at_rb': bb-ab, 'adiabatic': bb-aa}
    if not all(math.isfinite(value) for value in gaps.values()):
        raise ValueError('four-point energy gaps overflowed')
    return {'schema': 'cube-state-pair-four-point-1', 'status': 'analyzed',
            'states': {'a': asdict(state_a), 'b': asdict(state_b)},
            'energy_unit': energies.unit, 'energy_kind': 'total',
            'energies': {'E_A_RA': aa, 'E_B_RA': ba, 'E_B_RB': bb, 'E_A_RB': ab},
            'relaxation': terms, 'gaps': gaps,
            'diagnostics': {'negative_relaxation_a': terms['on_a'] < 0,
                            'negative_relaxation_b': terms['on_b'] < 0,
                            'minima_verified': False,
                            'state_identity_verified': False,
                            'environment_consistency_verified': False}}
