"""Dependency-free electronic-structure settings shared by the CLI and API."""
from dataclasses import dataclass
import math

METHODS = ('dipro', 'frozen_density', 'esid', 'orbital_gmh')


@dataclass(frozen=True)
class SCFConfig:
    backend: str = 'cpu'
    basis: str = 'def2-svp'
    xc: str = 'b3lyp'
    auxbasis: str = 'def2-universal-jkfit'
    unit: str = 'Angstrom'
    conv_tol: float = 1e-10
    max_cycle: int = 100
    grid_level: int = 3
    threads: int = 4
    max_memory_mb: int = 4000
    degeneracy_tolerance_ev: float = 0.01
    site_tolerance_ev: float = 0.001
    max_pair_leakage: float = 0.15

    def validate(self):
        if self.backend not in ('cpu', 'gpu') or self.unit not in ('Angstrom', 'Bohr'):
            raise ValueError('backend must be cpu/gpu and unit Angstrom/Bohr')
        for name in ('basis', 'xc', 'auxbasis'):
            if not isinstance(getattr(self, name), str) or not getattr(self, name).strip():
                raise ValueError(f'{name} must be a nonempty string')
        for name in ('max_cycle', 'threads', 'max_memory_mb'):
            if type(getattr(self, name)) is not int or getattr(self, name) <= 0:
                raise ValueError(f'{name} must be a positive integer')
        if type(self.grid_level) is not int or not 0 <= self.grid_level <= 9:
            raise ValueError('grid_level must be an integer in [0,9]')
        for name in ('conv_tol', 'degeneracy_tolerance_ev', 'site_tolerance_ev'):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
                raise ValueError(f'{name} must be positive and finite')
        value = self.max_pair_leakage
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or not 0 <= value < 1:
            raise ValueError('max_pair_leakage must be in [0,1)')
