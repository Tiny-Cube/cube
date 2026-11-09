"""Dependency-free, single-frame XYZ input; coordinates are never optimized."""
import math
from pathlib import Path


def normalize_atoms(atoms):
    try:
        result = [(str(symbol), tuple(float(x) for x in xyz)) for symbol, xyz in atoms]
    except (TypeError, ValueError) as exc:
        raise ValueError('atoms must contain element and three coordinates') from exc
    if not result or any(not symbol.strip() or len(xyz) != 3 or not all(map(math.isfinite, xyz)) for symbol, xyz in result):
        raise ValueError('each fragment needs atoms with finite three-dimensional coordinates')
    return result


def read_xyz(path):
    """Read exactly one XYZ frame, retaining its atom ordering."""
    lines = Path(path).read_text().splitlines()
    if len(lines) < 3:
        raise ValueError('XYZ needs an atom count, comment line, and atoms')
    n = int(lines[0])
    if n <= 0 or len(lines) < n+2 or any(s.strip() for s in lines[n+2:]):
        raise ValueError('XYZ must contain exactly one complete, nonempty frame')
    atoms = []
    for line in lines[2:n+2]:
        fields = line.split()
        if len(fields) != 4:
            raise ValueError('XYZ atom rows must contain element x y z')
        atoms.append((fields[0], tuple(map(float, fields[1:]))))
    return normalize_atoms(atoms)


def normalize_fragment(value):
    return read_xyz(value) if isinstance(value, (str, Path)) else normalize_atoms(value)
