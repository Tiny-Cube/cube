"""Strict TOML controls for postprocessing three existing PySCF checkpoints.

Parsing requires only the standard library. File paths are relative to the
control file. Explicit atom and MO numbers are one-based in TOML.
"""
from dataclasses import dataclass
import hashlib
from pathlib import Path
import re
import tomllib

from .config import SCFConfig

CHECKPOINT_METHODS = ('dipro', 'esid', 'orbital_gmh')


def _table(value, name, allowed):
    if not isinstance(value, dict):
        raise ValueError(f'[{name}] must be a TOML table')
    unknown = set(value)-set(allowed)
    if unknown:
        raise ValueError(f'unknown keys in [{name}]: {", ".join(sorted(unknown))}')
    return value


def _path(value, base, label):
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f'{label} must be a nonempty path string')
    path = Path(value).expanduser()
    return (path if path.is_absolute() else base/path).resolve()


def _indices(value, label):
    if isinstance(value, str):
        result = []
        for token in value.split(','):
            match = re.fullmatch(r'\s*([1-9][0-9]*)(?:\s*-\s*([1-9][0-9]*))?\s*', token)
            if not match:
                raise ValueError(f'{label}: use positive one-based numbers/ranges, e.g. "1-6, 9"')
            first, last = int(match[1]), int(match[2] or match[1])
            if last < first or last-first > 100000:
                raise ValueError(f'{label}: invalid or oversized range')
            result.extend(range(first, last+1))
    elif isinstance(value, list):
        result = value
    else:
        raise ValueError(f'{label} must be an integer array or range string')
    if not result or any(type(i) is not int or i < 1 for i in result):
        raise ValueError(f'{label} must contain positive one-based integers')
    if len(set(result)) != len(result):
        raise ValueError(f'{label} contains duplicate numbers')
    return [i-1 for i in result]


def validate_output(path, protected_paths):
    path = path.resolve()
    for source in protected_paths:
        if path == source.resolve() or (path.exists() and source.exists() and path.samefile(source)):
            raise ValueError('output must not overwrite a control file or input file')
    if path.exists() and not path.is_file():
        raise ValueError('output path must refer to a file')
    ancestor = path.parent
    while not ancestor.exists():
        ancestor = ancestor.parent
    if not ancestor.is_dir():
        raise ValueError('output parent is a file, not a directory')
    return path


def _selection(value, label, allowed):
    if not isinstance(value, list) or not value or any(not isinstance(x, str) for x in value):
        raise ValueError(f'{label} must be a nonempty array of strings')
    if len(set(value)) != len(value) or any(x not in allowed for x in value):
        raise ValueError(f'{label} must select unique values from {", ".join(allowed)}')
    return tuple(value)


@dataclass(frozen=True)
class TransferControl:
    checkpoints: dict
    config: SCFConfig
    methods: tuple
    carriers: tuple
    assume_equivalent: bool
    assume_converged: bool
    orbital_indices: dict | None
    atom_indices: tuple | None
    output: Path
    provenance: dict

    def summary(self):
        return {'validation': 'paths_and_settings_only_no_checkpoint_matrices_read',
                'control': self.provenance,
                'checkpoints': {k: str(p) for k, p in self.checkpoints.items()},
                'backend': self.config.backend, 'methods': self.methods, 'carriers': self.carriers,
                'assume_equivalent': self.assume_equivalent, 'assume_converged': self.assume_converged,
                'orbital_indices_zero_based': self.orbital_indices,
                'atom_indices_zero_based': self.atom_indices, 'output': str(self.output)}


def load_control(filename):
    path = Path(filename).expanduser().resolve()
    raw = path.read_bytes()
    data = tomllib.loads(raw.decode('utf-8'))
    _table(data, 'root', ('schema_version', 'input', 'calculation', 'output', 'orbitals', 'fragments'))
    if type(data.get('schema_version')) is not int or data['schema_version'] != 1:
        raise ValueError('schema_version = 1 is required')
    inputs = _table(data.get('input'), 'input', ('fragment_a', 'fragment_b', 'dimer'))
    if set(inputs) != {'fragment_a', 'fragment_b', 'dimer'}:
        raise ValueError('[input] requires fragment_a, fragment_b and dimer checkpoint paths')
    checkpoints = {k: _path(v, path.parent, f'input.{k}') for k, v in inputs.items()}
    if len(set(checkpoints.values())) != 3:
        raise ValueError('A, B and dimer must refer to three distinct checkpoint files')
    for name, source in checkpoints.items():
        if not source.is_file():
            raise ValueError(f'input.{name} checkpoint not found: {source}')
        if source.suffix.lower() in ('.out', '.log', '.txt'):
            raise ValueError('text SCF logs lack complete orbital matrices; specify a PySCF .chk file')
    names = {'backend', 'threads', 'degeneracy_tolerance_ev', 'site_tolerance_ev', 'max_pair_leakage'}
    calculation = _table(data.get('calculation', {}), 'calculation', names | {'methods', 'carriers', 'assume_equivalent', 'assume_converged'})
    config = SCFConfig(**{k: v for k, v in calculation.items() if k in names})
    config.validate()
    methods = _selection(calculation.get('methods', ['dipro']), 'calculation.methods', CHECKPOINT_METHODS)
    carriers = _selection(calculation.get('carriers', ['hole', 'electron']), 'calculation.carriers', ('hole', 'electron'))
    equivalent, converged = (calculation.get(k, False) for k in ('assume_equivalent', 'assume_converged'))
    if type(equivalent) is not bool or type(converged) is not bool:
        raise ValueError('assume_equivalent and assume_converged must be boolean')
    orbital_indices = None
    if 'orbitals' in data:
        orbital_indices = {}
        for carrier, selection in _table(data['orbitals'], 'orbitals', carriers).items():
            selection = _table(selection, f'orbitals.{carrier}', ('a', 'b'))
            if set(selection) != {'a', 'b'}:
                raise ValueError(f'orbitals.{carrier} requires a and b')
            orbital_indices[carrier] = tuple(_indices(selection[k], f'orbitals.{carrier}.{k}') for k in ('a', 'b'))
    atom_indices = None
    if 'fragments' in data:
        fragments = _table(data['fragments'], 'fragments', ('a_atoms', 'b_atoms'))
        if set(fragments) != {'a_atoms', 'b_atoms'}:
            raise ValueError('[fragments] requires a_atoms and b_atoms in monomer atom order')
        atom_indices = tuple(_indices(fragments[k], f'fragments.{k}') for k in ('a_atoms', 'b_atoms'))
        if set(atom_indices[0]) & set(atom_indices[1]):
            raise ValueError('fragment atom selections overlap')
    output = _table(data.get('output'), 'output', ('file',))
    output_path = validate_output(_path(output.get('file'), path.parent, 'output.file'), [path, *checkpoints.values()])
    provenance = {'schema_version': 1, 'path': str(path), 'sha256': hashlib.sha256(raw).hexdigest(),
                  'atom_index_base': 1, 'orbital_index_base': 1}
    return TransferControl(checkpoints, config, methods, carriers, equivalent, converged,
                           orbital_indices, atom_indices, output_path, provenance)
