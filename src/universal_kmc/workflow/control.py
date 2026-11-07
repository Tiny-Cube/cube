"""Parse $section ... $end input controls without scientific dependencies."""
from dataclasses import dataclass
import hashlib
import math
from pathlib import Path
import re

from universal_kmc.quantum.config import SCFConfig

SECTIONS = ('get_dimer', 'two-mol', 'one-mol', 'reorganization-energy')
COMMON = {'qc-engine', 'method', 'basis-set', 'auxbasis', 'threads', 'grid-level',
          'conv-tol', 'max-cycle', 'max-memory-mb', 'input', 'output'}
ALLOWED = {
    'get_dimer': {'tpr', 'gro', 'xtc', 'frame', 'cutoff', 'fragment-file', 'output', 'format'},
    'one-mol': COMMON | {'molecules', 'xyz'},
    'two-mol': COMMON | {'dimers', 'transfer-integrals-method', 'carriers', 'assume-equivalent'},
    'reorganization-energy': COMMON | {'molecules', 'xyz', 're-method', 'environment', 'carriers', 'opt-maxsteps',
                                      'mm-source', 'mm-charges', 'mm-cutoff', 'mm-frozen'},
}


def _uncomment(line):
    quote = None
    escaped = False
    for index, char in enumerate(line):
        if escaped:
            escaped = False
        elif char == '\\' and quote:
            escaped = True
        elif quote:
            if char == quote:
                quote = None
        elif char in ('"', "'"):
            quote = char
        elif char == '#':
            return line[:index].strip()
    if quote:
        raise ValueError('unterminated quoted value')
    return line.strip()


def _value(text):
    if text.startswith(('"', "'")):
        if text[-1:] != text[0]:
            raise ValueError('quoted values must occupy the whole right-hand side')
        return text[1:-1]
    return text


def _csv(value):
    result = [s.strip() for s in value.split(',')]
    if not result or any(not s for s in result) or len(set(result)) != len(result):
        raise ValueError('expected a nonempty comma-separated list without duplicates')
    return result


def _integer(value, name, minimum=1):
    if not re.fullmatch(r'-?\d+', value):
        raise ValueError(f'{name} must be an integer')
    number = int(value)
    if number < minimum:
        raise ValueError(f'{name} must be >= {minimum}')
    return number


def _positive(value, name):
    number = float(value)
    if not math.isfinite(number) or number <= 0:
        raise ValueError(f'{name} must be positive and finite')
    return number


def _boolean(value, name):
    if value.lower() not in ('true', 'false'):
        raise ValueError(f'{name} must be true or false')
    return value.lower() == 'true'


def resolve_path(value, base, name, *, required_file=False):
    if not value.strip():
        raise ValueError(f'{name} needs a file/directory path')
    path = Path(value).expanduser()
    path = (path if path.is_absolute() else base/path).resolve()
    if required_file and not path.is_file():
        raise ValueError(f'{name} file not found: {path}')
    return path


@dataclass(frozen=True)
class WorkflowInput:
    path: Path
    sha256: str
    sections: dict

    def plan(self, only=None):
        selected = set(only or self.sections)
        if not selected <= self.sections.keys():
            raise ValueError('--only must name a section present in input')
        order = [s for s in ('get_dimer', 'one-mol', 'two-mol', 'reorganization-energy') if s in selected]
        return {'schema': 'cube-workflow-1', 'control': str(self.path), 'control_sha256': self.sha256,
                'tasks': order, 'sections': self.sections,
                'note': 'plan only; no trajectory extraction, SCF, optimization or GPU probe performed'}


def load_input(filename='input'):
    path = Path(filename).expanduser().resolve()
    raw = path.read_bytes()
    sections, active = {}, None
    for line_number, original in enumerate(raw.decode('utf-8-sig').splitlines(), 1):
        try:
            line = _uncomment(original)
            if not line:
                continue
            if line.startswith('$'):
                name = line[1:].strip().lower()
                if name == 'end':
                    if active is None:
                        raise ValueError('$end without an open section')
                    active = None
                else:
                    if name not in SECTIONS or active is not None or name in sections:
                        raise ValueError(f'unknown, nested or duplicate section ${name}')
                    sections[name], active = {}, name
                continue
            if active is None or '=' not in line:
                raise ValueError('expected key = value inside a section')
            key, value = (s.strip() for s in line.split('=', 1))
            key = key.lower()
            if key not in ALLOWED[active] or key in sections[active]:
                raise ValueError(f'unknown or duplicate key {key!r} in ${active}')
            sections[active][key] = _value(value)
        except ValueError as exc:
            raise ValueError(f'{path}:{line_number}: {exc}') from exc
    if active is not None:
        raise ValueError(f'missing $end for ${active}')
    if not sections:
        raise ValueError('input contains no tasks')
    normalized = {}
    if 'get_dimer' in sections:
        s = sections['get_dimer']
        n = {key: str(resolve_path(s.get(key, ''), path.parent, key, required_file=True))
             for key in ('tpr', 'gro', 'xtc', 'fragment-file')}
        n['frame'] = _integer(s.get('frame') or '-1', 'frame', -1)
        n['cutoff'] = _positive(s.get('cutoff', '6'), 'cutoff')
        n['output'] = str(resolve_path(s.get('output', 'mol'), path.parent, 'output'))
        n['format'] = [v.lower() for v in _csv(s.get('format', 'xyz'))]
        if len(set(n['format'])) != len(n['format']) or not set(n['format']) <= {'xyz', 'gjf'}:
            raise ValueError('format must select xyz and/or gjf')
        if any(Path(n['output']) == Path(n[k]) for k in ('tpr', 'gro', 'xtc', 'fragment-file')):
            raise ValueError('output directory must not overwrite an input file')
        normalized['get_dimer'] = n
    inherited = sections.get('one-mol', {}).get('qc-engine', 'pyscf')
    for name in SECTIONS[1:]:
        if name not in sections:
            continue
        s = sections[name]
        engine = s.get('qc-engine', inherited).lower()
        if engine not in ('pyscf', 'gpu-pyscf'):
            raise ValueError(f'${name}: qc-engine must be pyscf or gpu-pyscf')
        n = {'qc-engine': engine, 'method': s.get('method', 'b3lyp'),
             'basis-set': s.get('basis-set', '6-31G(d,p)'),
             'threads': _integer(s.get('threads', '4'), 'threads'),
             'grid-level': _integer(s.get('grid-level', '3'), 'grid-level', 0),
             'conv-tol': _positive(s.get('conv-tol', '1e-10'), 'conv-tol'),
             'max-cycle': _integer(s.get('max-cycle', '100'), 'max-cycle'),
             'max-memory-mb': _integer(s.get('max-memory-mb', '4000'), 'max-memory-mb'),
             'auxbasis': s.get('auxbasis', 'def2-universal-jkfit')}
        scf_config(n).validate()
        for key in ('input', 'output'):
            if key in s:
                n[key] = str(resolve_path(s[key], path.parent, key, required_file=(key == 'input')))
        if 'xyz' in s:
            n['xyz'] = [str(resolve_path(v, path.parent, 'xyz', required_file=True)) for v in _csv(s['xyz'])]
        if not ('get_dimer' in normalized or 'input' in n or 'xyz' in n):
            raise ValueError(f'${name} needs $get_dimer, input = manifest.json, or xyz = ...')
        selection = 'dimers' if name == 'two-mol' else 'molecules'
        if selection in s:
            n[selection] = _csv(s[selection])
        if name == 'two-mol':
            methods = [v.lower() for v in _csv(s.get('transfer-integrals-method', 'DIPRO'))]
            if not set(methods) <= {'dipro', 'esid', 'orbital_gmh'}:
                raise ValueError('transfer-integrals-method must select DIPRO, ESID or orbital_gmh')
            n['transfer-integrals-method'] = methods
            n['assume-equivalent'] = _boolean(s.get('assume-equivalent', 'false'), 'assume-equivalent')
        if name in ('two-mol', 'reorganization-energy'):
            n['carriers'] = [v.lower() for v in _csv(s.get('carriers', 'hole,electron'))]
            if not set(n['carriers']) <= {'hole', 'electron'}:
                raise ValueError('carriers must select hole and/or electron')
        if name == 'reorganization-energy':
            n['environment'] = s.get('environment', 'vacuum').lower()
            n['re-method'] = s.get('re-method', '4-points').lower()
            n['opt-maxsteps'] = _integer(s.get('opt-maxsteps', '100'), 'opt-maxsteps')
            if n['environment'] not in ('vacuum', 'qmmm') or n['re-method'] != '4-points':
                raise ValueError('use environment=vacuum or qmmm and re-method=4-points')
            mm_keys = {'mm-source', 'mm-charges', 'mm-cutoff', 'mm-frozen'}
            if n['environment'] == 'vacuum' and mm_keys.intersection(s):
                raise ValueError('MM settings require environment=qmmm')
            if n['environment'] == 'qmmm':
                n['mm-frozen'] = _boolean(s.get('mm-frozen', 'true'), 'mm-frozen')
                if not n['mm-frozen']:
                    raise ValueError('QM/MM four-point currently requires mm-frozen=true')
                source = s.get('mm-source', 'file' if 'mm-charges' in s else 'manifest').lower()
                if source not in ('manifest', 'file'):
                    raise ValueError('mm-source must be manifest or file')
                n['mm-source'] = source
                if source == 'file':
                    if 'mm-cutoff' in s:
                        raise ValueError('mm-cutoff applies only to mm-source=manifest; file charges are all included')
                    n['mm-charges'] = str(resolve_path(s.get('mm-charges', ''), path.parent,
                                                       'mm-charges', required_file=True))
                    if 'xyz' not in n or len(n['xyz']) != 1:
                        raise ValueError('explicit mm-charges requires exactly one xyz in the same coordinate frame')
                else:
                    if 'mm-charges' in s or 'xyz' in n:
                        raise ValueError('mm-source=manifest requires a GROMACS extraction manifest, without xyz/mm-charges')
                    n['mm-cutoff'] = _positive(s.get('mm-cutoff', '12'), 'mm-cutoff')
        normalized[name] = n
    if 'two-mol' in normalized and 'one-mol' in normalized:
        for key in ('method', 'basis-set', 'auxbasis', 'grid-level', 'conv-tol'):
            if str(normalized['two-mol'][key]).lower() != str(normalized['one-mol'][key]).lower():
                raise ValueError(f'DIPRO requires matching one-mol/two-mol {key}')
    return WorkflowInput(path, hashlib.sha256(raw).hexdigest(), normalized)


def scf_config(settings):
    return SCFConfig(backend='gpu' if settings['qc-engine'] == 'gpu-pyscf' else 'cpu',
                     basis=settings['basis-set'], xc=settings['method'], auxbasis=settings['auxbasis'],
                     threads=settings['threads'], grid_level=settings['grid-level'],
                     conv_tol=settings['conv-tol'], max_cycle=settings['max-cycle'],
                     max_memory_mb=settings['max-memory-mb'])
