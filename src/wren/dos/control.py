"""Dependency-free DOS control normalization for TOML and Wren sections."""
import math
from pathlib import Path

from wren.quantum.config import SCFConfig

FIELDS = {'do', 'mode', 'input', 'xyz', 'output', 'molecules', 'sigma', 'points', 'emin', 'emax',
          'kmesh', 'kmesh-study', 'spin-degeneracy', 'dos-backend', 'qc-engine', 'method', 'basis-set', 'auxbasis',
          'threads', 'grid-level', 'conv-tol', 'max-cycle', 'max-memory-mb', 'pseudo', 'precision', 'max-atoms',
          'site-projections', 'degeneracy-tol-ev', 'element-projections', 'solver',
          'kpm-moments', 'kpm-trace', 'kpm-vectors', 'kpm-seed', 'kpm-bounds-ev',
          'periodic-reference', 'reference-file', 'reference-label', 'block-size-frames'}


def _mesh(parts):
    try:
        mesh = [int(str(n).strip()) for n in parts]
    except (TypeError, ValueError) as exc:
        raise ValueError('DOS kmesh needs three integers') from exc
    if len(mesh) != 3 or min(mesh) < 1:
        raise ValueError('DOS kmesh needs three positive integers')
    return mesh


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
    molecular = mode in ('monomer-orbitals', 'vertical-dscf')
    if not molecular and mode not in ('periodic', 'tight-binding'):
        raise ValueError('DOS mode must be monomer-orbitals, vertical-dscf, periodic or tight-binding')
    result = {'mode': mode}
    solver = str(raw.get('solver', 'diagonalization')).lower()
    if 'solver' in raw or any(key.startswith('kpm-') for key in raw):
        if mode != 'tight-binding' or solver not in ('diagonalization', 'kpm'):
            raise ValueError('DOS solver requires tight-binding and diagonalization or kpm')
        result['solver'] = solver
    if solver == 'kpm':
        result['kpm-moments'] = number('kpm-moments', 128, integer=True, minimum=2)
        trace = str(raw.get('kpm-trace', 'stochastic')).lower()
        if trace not in ('stochastic', 'exact'):
            raise ValueError('kpm-trace must be stochastic or exact')
        result['kpm-trace'] = trace
        if trace == 'stochastic':
            result['kpm-vectors'] = number('kpm-vectors', 16, integer=True, minimum=2)
            try:
                seed = int(str(raw.get('kpm-seed', 0)))
            except ValueError as exc:
                raise ValueError('kpm-seed must be a nonnegative integer') from exc
            if seed < 0:
                raise ValueError('kpm-seed must be a nonnegative integer')
            result['kpm-seed'] = seed
        elif {'kpm-vectors', 'kpm-seed'} & raw.keys():
            raise ValueError('kpm-vectors and kpm-seed require stochastic trace')
        if 'kpm-bounds-ev' in raw:
            bounds = raw['kpm-bounds-ev']
            bounds = bounds if isinstance(bounds, list) else str(bounds).split(',')
            values = [float(v) for v in bounds]
            if len(values) != 2 or not all(map(math.isfinite, values)) or values[0] >= values[1]:
                raise ValueError('kpm-bounds-ev needs two increasing finite energies')
            result['kpm-bounds-ev'] = values
    elif any(key.startswith('kpm-') for key in raw):
        raise ValueError('kpm parameters require solver=kpm')
    def flag(name):
        value = str(raw.get(name, False)).strip().lower()
        if value not in ('true', 'false', '1', '0'):
            raise ValueError(f'{name} must be true/false or 1/0')
        return value in ('true', '1')
    if 'element-projections' in raw:
        if mode not in ('periodic', 'monomer-orbitals'):
            raise ValueError('element-projections requires periodic or monomer-orbitals mode')
        result['element-projections'] = flag('element-projections')
    if 'site-projections' in raw or 'degeneracy-tol-ev' in raw:
        if mode != 'tight-binding':
            raise ValueError('site-projections and degeneracy-tol-ev require tight-binding mode')
        result['site-projections'] = flag('site-projections')
        if result['site-projections'] and solver == 'kpm':
            raise ValueError('site-projections and participation require diagonalization')
        if result['site-projections']:
            result['degeneracy-tol-ev'] = number('degeneracy-tol-ev', 1e-8)
        elif 'degeneracy-tol-ev' in raw:
            raise ValueError('degeneracy-tol-ev requires enabled site-projections')
    if {'periodic-reference', 'reference-file', 'reference-label'} & raw.keys():
        reference = str(raw.get('periodic-reference', 'external')).strip().lower()
        if mode != 'periodic' or reference not in ('external', 'frame-vbm', 'cell-mean-local'):
            raise ValueError('periodic-reference requires external, frame-vbm or cell-mean-local in periodic mode')
        if (reference == 'external') != ('reference-file' in raw):
            raise ValueError('reference-file must be supplied for an external reference only')
        if 'reference-label' in raw:
            label = str(raw['reference-label']).strip()
            if reference != 'external' or not label:
                raise ValueError('reference-label describes a nonempty external energy reference')
            result['reference-label'] = label
        result['periodic-reference'] = reference
    if 'block-size-frames' in raw:
        if mode not in ('periodic', 'monomer-orbitals'):
            raise ValueError('block-size-frames requires periodic or monomer-orbitals DOS')
        result['block-size-frames'] = number('block-size-frames', None, integer=True, minimum=1)
    for name in ('input', 'xyz', 'output', 'reference-file'):
        if name in raw or name == 'output':
            p = Path(raw.get(name, 'dos-results')).expanduser()
            p = (base/p).resolve() if not p.is_absolute() else p.resolve()
            if name != 'output' and not p.is_file():
                raise ValueError(f'DOS {name} file does not exist: {p}')
            result[name] = str(p)
    if sum(name in result for name in ('input', 'xyz')) != 1:
        raise ValueError('DOS requires exactly one input manifest/network or xyz')
    if not molecular and 'xyz' in result:
        raise ValueError('periodic/tight-binding DOS requires an input dataset with a defined cell/model')
    for name in ('input', 'xyz', 'reference-file'):
        if result.get(name) == result['output']:
            raise ValueError('DOS output must not overwrite an input')
    if 'molecules' in raw:
        value = raw['molecules']
        selected = list(map(str, value)) if isinstance(value, list) else [s.strip() for s in str(value).split(',')]
        if not selected or any(not s for s in selected) or len(set(selected)) != len(selected):
            raise ValueError('DOS molecules need a unique nonempty list')
        if not molecular:
            raise ValueError('molecule selection applies only to monomer DOS')
        result['molecules'] = selected
    if solver == 'kpm':
        if 'sigma' in raw:
            raise ValueError('KPM uses Jackson damping; set kpm-moments for energy resolution')
    else:
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
    result['kmesh'] = _mesh(mesh)
    if 'kmesh-study' in raw:
        entries = raw['kmesh-study']
        entries = entries if isinstance(entries, list) else str(entries).split(',')
        if mode != 'periodic' or len(entries) < 2:
            raise ValueError('kmesh-study needs periodic mode and at least two distinct meshes')
        meshes = [_mesh(item if isinstance(item, list) else str(item).lower().split('x')) for item in entries]
        for old, new in zip(meshes, meshes[1:]):
            if old == new or any(a < b for a, b in zip(new, old)):
                raise ValueError('each study mesh must refine at least one axis and never coarsen another')
        result['kmesh-study'] = meshes
    result['spin-degeneracy'] = number('spin-degeneracy', 1, integer=True, minimum=1)
    if result['spin-degeneracy'] not in (1, 2):
        raise ValueError('spin-degeneracy must be 1 or 2')
    if mode == 'vertical-dscf' and (result['spin-degeneracy'] != 1 or result['kmesh'] != [1, 1, 1]):
        raise ValueError('vertical-dscf counts one transition per molecule and channel, with no k mesh or spin multiplier')
    result['dos-backend'] = str(raw.get('dos-backend', 'cpu')).lower()
    result['qc-engine'] = str(raw.get('qc-engine', 'pyscf')).lower()
    if result['dos-backend'] not in ('cpu', 'gpu') or result['qc-engine'] not in ('pyscf', 'gpu-pyscf'):
        raise ValueError('invalid DOS backend or QC engine')
    default_basis = {'periodic': 'gth-szv', 'vertical-dscf': 'aug-cc-pvdz'}.get(mode, '6-31G(d,p)')
    for key, default in (('method', 'b3lyp'), ('basis-set', default_basis),
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
