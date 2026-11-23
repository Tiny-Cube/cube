"""Dependency-free DOS control normalization for TOML and Wren sections."""
import math
from pathlib import Path

from wren.quantum.config import SCFConfig

FIELDS = {'do', 'mode', 'input', 'xyz', 'output', 'molecules', 'sigma', 'points', 'emin', 'emax',
          'kmesh', 'spin-degeneracy', 'dos-backend', 'qc-engine', 'method', 'basis-set', 'auxbasis',
          'threads', 'grid-level', 'conv-tol', 'max-cycle', 'max-memory-mb', 'pseudo', 'precision', 'max-atoms'}


def normalize(raw, base):
    if set(raw)-FIELDS:
        raise ValueError('unknown DOS fields: '+', '.join(sorted(set(raw)-FIELDS)))
    def number(key, default, *, integer=False, minimum=0, strict=False):
        value = raw.get(key, default)
        if isinstance(value, bool):
            raise ValueError(f'invalid DOS {key}')
        result = float(value)
        if not math.isfinite(result) or (result <= minimum if strict else result < minimum):
            raise ValueError(f'invalid DOS {key}')
        if integer and result != int(result):
            raise ValueError(f'DOS {key} must be integer')
        return int(result) if integer else result
    mode = str(raw.get('mode', 'monomer-orbitals')).lower()
    if mode not in ('monomer-orbitals', 'periodic', 'tight-binding'):
        raise ValueError('DOS mode must be monomer-orbitals, periodic or tight-binding')
    result = {'mode': mode}
    for name in ('input', 'xyz', 'output'):
        if name in raw or name == 'output':
            p = Path(raw.get(name, 'dos-results')).expanduser()
            p = (base/p).resolve() if not p.is_absolute() else p.resolve()
            if name != 'output' and not p.is_file():
                raise ValueError(f'DOS {name} file does not exist: {p}')
            result[name] = str(p)
    if sum(name in result for name in ('input', 'xyz')) != 1:
        raise ValueError('DOS requires exactly one input manifest/network or xyz')
    if mode != 'monomer-orbitals' and 'xyz' in result:
        raise ValueError('periodic/tight-binding DOS requires an input dataset with a defined cell/model')
    for name in ('input', 'xyz'):
        if result.get(name) == result['output']:
            raise ValueError('DOS output must not overwrite an input')
    if 'molecules' in raw:
        value = raw['molecules']
        selected = list(map(str, value)) if isinstance(value, list) else [s.strip() for s in str(value).split(',')]
        if not selected or any(not s for s in selected) or len(set(selected)) != len(selected):
            raise ValueError('DOS molecules need a unique nonempty list')
        if mode != 'monomer-orbitals':
            raise ValueError('molecule selection applies only to monomer DOS')
        result['molecules'] = selected
    result['sigma'] = number('sigma', .05, strict=True)
    if 'points' in raw:
        result['points'] = number('points', None, integer=True, minimum=2)
    for name in ('emin', 'emax'):
        if name in raw:
            result[name] = float(raw[name])
            if not math.isfinite(result[name]):
                raise ValueError('DOS energy bounds must be finite')
    if 'emin' in result and 'emax' in result and result['emin'] >= result['emax']:
        raise ValueError('DOS emin must be below emax')
    mesh = raw.get('kmesh', [1, 1, 1])
    mesh = mesh if isinstance(mesh, list) else str(mesh).split(',')
    try:
        result['kmesh'] = [int(str(n).strip()) for n in mesh]
    except ValueError as exc:
        raise ValueError('DOS kmesh needs three integers') from exc
    if len(mesh) != 3 or any(n < 1 for n in result['kmesh']):
        raise ValueError('DOS kmesh needs three positive integers')
    result['spin-degeneracy'] = number('spin-degeneracy', 1, integer=True, minimum=1)
    if result['spin-degeneracy'] not in (1, 2):
        raise ValueError('spin-degeneracy must be 1 or 2')
    result['dos-backend'] = str(raw.get('dos-backend', 'cpu')).lower()
    result['qc-engine'] = str(raw.get('qc-engine', 'pyscf')).lower()
    if result['dos-backend'] not in ('cpu', 'gpu') or result['qc-engine'] not in ('pyscf', 'gpu-pyscf'):
        raise ValueError('invalid DOS backend or QC engine')
    if mode == 'periodic' and result['qc-engine'] != 'pyscf':
        raise ValueError('periodic GPU DOS is not validated; explicitly request pyscf')
    for key, default in (('method', 'b3lyp'), ('basis-set', 'gth-szv' if mode == 'periodic' else '6-31G(d,p)'),
                         ('auxbasis', 'auto' if mode == 'periodic' else 'def2-universal-jkfit'), ('pseudo', 'gth-pbe')):
        result[key] = str(raw.get(key, default))
    for key, default, minimum in (('threads', 4, 1), ('grid-level', 3, 0), ('max-cycle', 100, 1),
                                   ('max-memory-mb', 4000, 1), ('max-atoms', 500, 1)):
        result[key] = number(key, default, integer=True, minimum=minimum)
    result['conv-tol'] = number('conv-tol', 1e-9, strict=True)
    result['precision'] = number('precision', 1e-8, strict=True)
    SCFConfig(basis=result['basis-set'], xc=result['method'], auxbasis=result['auxbasis'],
              grid_level=result['grid-level'], threads=result['threads'], max_cycle=result['max-cycle'],
              max_memory_mb=result['max-memory-mb'], conv_tol=result['conv-tol']).validate()
    return result
