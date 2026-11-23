"""Parse $section ... $end input controls without scientific dependencies."""
from dataclasses import dataclass, field
import hashlib
import math
from pathlib import Path
import re

from wren.quantum.config import SCFConfig

from .reorganization_methods import ANALYSIS_METHODS
from .calculations import (CALCULATIONS, EXTRACTION, KMC, DOS, SURFACE_DYNAMICS, QUANTUM_SECTIONS, route_sections, section_names)
from wren.transport.config import FIELDS as KMC_FIELDS
from wren.dos.control import FIELDS as DOS_FIELDS

COMMON = {'do', 'qc-engine', 'method', 'basis-set', 'auxbasis', 'threads', 'grid-level',
          'conv-tol', 'max-cycle', 'max-memory-mb', 'input', 'output'}
ALLOWED = {
    DOS: DOS_FIELDS,
    KMC: KMC_FIELDS,
    SURFACE_DYNAMICS: {'do', 'input', 'output'},
    EXTRACTION: {'do', 'tpr', 'gro', 'xtc', 'frame', 'frame_step', 'dimers_cutoff', 'fragment-file',
                 'monomers-qmmm-cutoff', 'output', 'format', 'lammps-data', 'lammps-dump',
                 'lammps-units', 'lammps-timestep', 'element-map'},
    'one-mol': COMMON | {'molecules', 'xyz'},
    'two-mol': COMMON | {'dimers', 'transfer-integrals-method', 'carriers', 'assume-equivalent'},
    'reorganization-energy': COMMON | {'molecules', 'xyz', 're-method', 'environment', 'carriers', 'opt-maxsteps',
                                      'mm-source', 'mm-charges', 'mm-cutoff', 'mm-frozen',
                                      'mm-vdw', 'mm-vdw-file', 'gmx-executable', 'energy-source',
                                      'state-a', 'state-b', 'root-gap-ev', 'root-min-overlap',
                                      'coordinates', 'topology', 'forcefield', 'electrostatics',
                                      'ensemble-a', 'ensemble-b', 'temperature', 'block-size'},
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


def needs_extraction(settings):
    return (not {'input', 'xyz', 'coordinates'}.intersection(settings)
            and settings.get('re-method', '').lower() != 'linear-response')


def _normalize_analysis(s, base):
    """Return a data-analysis request, or None for a quantum calculation."""
    method = s.get('re-method', '4-points').lower()
    source = s.get('energy-source', 'calculate' if method == '4-points' else 'external').lower()
    if method not in ANALYSIS_METHODS or source not in ('calculate', 'external'):
        raise ValueError(f'choose re-method from {", ".join(ANALYSIS_METHODS)} and energy-source=calculate or external')
    if source == 'calculate' and method != '4-points':
        raise ValueError(f'{method} needs external input data')
    if source == 'calculate':
        if {'ensemble-a', 'ensemble-b', 'temperature', 'block-size'}.intersection(s):
            raise ValueError('ensemble settings require re-method=linear-response')
        return None
    allowed = {'do', 're-method', 'energy-source', 'output'}
    allowed.update({'input'} if not ANALYSIS_METHODS[method].ensemble_input else
                   {'state-a', 'state-b', 'ensemble-a', 'ensemble-b', 'temperature', 'block-size', 'environment'})
    if extra := s.keys() - allowed:
        raise ValueError(f'{method} external analysis cannot use: {", ".join(sorted(extra))}')
    n = {'re-method': method, 'energy-source': source,
         'output': str(resolve_path(s.get('output', 'cube-results/reorganization'), base, 'output'))}
    if not ANALYSIS_METHODS[method].ensemble_input:
        n['input'] = str(resolve_path(s.get('input', ''), base, 'input', required_file=True))
    else:
        from .state_pairs import explicit_pair
        n.update(explicit_pair(s))
        n['temperature'] = _positive(s.get('temperature', ''), 'temperature (K)')
        n['environment'] = s.get('environment', '').strip()
        if not n['environment']:
            raise ValueError('linear-response requires an explicit environment description')
        for key in ('ensemble-a', 'ensemble-b'):
            n[key] = str(resolve_path(s.get(key, ''), base, key, required_file=True))
        if Path(n['ensemble-a']).samefile(n['ensemble-b']):
            raise ValueError('A and B must refer to separate equilibrium ensemble files')
        if 'block-size' in s:
            n['block-size'] = _integer(s['block-size'], 'block-size')
    return n


@dataclass(frozen=True)
class WorkflowInput:
    path: Path
    sha256: str
    sections: dict
    inactive_sections: dict = field(default_factory=dict)
    module_switches: dict = field(default_factory=dict)

    def plan(self, only=None):
        entries = {name: self.sections[name] for name in CALCULATIONS if name in self.sections}
        available = self.sections.keys() | self.inactive_sections.keys()
        route = route_sections(available, entries, only, module_switches=self.module_switches)
        # An input-less quantum module needs the selected entry's shared extraction.
        if EXTRACTION in self.sections and any(
                task in QUANTUM_SECTIONS and needs_extraction(self.sections[task]) for task in route['tasks']):
            if EXTRACTION not in route['tasks']:
                owners = [name for name in route['enabled_calculations']
                          if EXTRACTION in CALCULATIONS[name].modules
                          and entries[name].get(EXTRACTION, 1) == 1]
                if any('standalone' in route['module_owners'].get(task, []) for task in route['tasks']
                       if task in QUANTUM_SECTIONS and needs_extraction(self.sections[task])):
                    owners.append('standalone')
                if self.module_switches.get(EXTRACTION) == 0:
                    raise ValueError(f'${EXTRACTION} do=0 but a selected module needs extracted data; '
                                     'enable extraction or provide input/xyz')
                if entries and not owners:
                    raise ValueError('selected calculation needs an explicit data source or an owned extraction module')
                route['tasks'].insert(0, EXTRACTION)
                route['skipped_modules'].pop(EXTRACTION, None)
                route['module_owners'][EXTRACTION] = owners or ['legacy']
        return {'schema': 'cube-workflow-1', 'control': str(self.path), 'control_sha256': self.sha256,
                **route, 'sections': self.sections,
                'inactive_sections': self.inactive_sections,
                'note': 'plan only; no trajectory extraction, SCF, optimization or GPU probe performed'}

    def execution_plan(self, only=None):
        plan = self.plan(only)
        for name in plan['calculation_tasks']:
            for module, selected in self.sections[name].items():
                if module == 'do' or selected == 0:
                    continue
                if self.module_switches.get(module) == 0:
                    raise ValueError(f'${name} requests {module}=1 but ${module} do=0; resolve the conflicting switches')
                if module not in self.sections:
                    raise ValueError(f'${name} requests {module}=1 but ${module} settings are missing; no MD/QC tasks were started')
            missing = set(CALCULATIONS[name].required_modules)-set(plan['tasks'])
            if missing:
                raise ValueError(f'${name} requires {", ".join("$"+m for m in sorted(missing))}; '
                                 'no MD/QC tasks were started')
        if plan['unsupported_tasks']:
            details = '; '.join(f'${s}: {CALCULATIONS[s].unavailable_reason}'
                                for s in plan['unsupported_tasks'])
            raise NotImplementedError(
                f'{details}; these calculation entries are not yet connected to cube input; '
                'no MD/QC tasks were started. Use --only for a supported shared module or disable the entry.')
        if KMC in plan['tasks']:
            from wren.transport.config import preflight
            kmc = self.sections[KMC]
            preflight(kmc, plan['tasks'])
            if not ('input' in kmc or 'reorganization' in kmc) and 'reorganization-energy' in plan['tasks']:
                reorg = self.sections['reorganization-energy']
                if 'carriers' not in reorg:
                    raise ValueError('KMC automatic assembly needs carrier four-point results; '
                                     'supply a carrier result explicitly or use carriers in $reorganization-energy')
        return plan


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
                if name == 'get_dimer':  # Historical controls remain readable.
                    name = EXTRACTION
                if name == 'end':
                    if active is None:
                        raise ValueError('$end without an open section')
                    active = None
                else:
                    if name not in section_names() or active is not None or name in sections:
                        raise ValueError(f'unknown, nested or duplicate section ${name}')
                    sections[name], active = {}, name
                continue
            if active is None or '=' not in line:
                raise ValueError('expected key = value inside a section')
            key, value = (s.strip() for s in line.split('=', 1))
            key = key.lower()
            if active == EXTRACTION and key == 'cutoff':
                key = 'dimers_cutoff'
            allowed = {'do', *CALCULATIONS[active].modules} if active in CALCULATIONS else ALLOWED[active]
            if key not in allowed or key in sections[active]:
                raise ValueError(f'unknown or duplicate key {key!r} in ${active}')
            sections[active][key] = _value(value)
        except ValueError as exc:
            raise ValueError(f'{path}:{line_number}: {exc}') from exc
    if active is not None:
        raise ValueError(f'missing $end for ${active}')
    if not sections:
        raise ValueError('input contains no tasks')
    normalized = {}
    for name in CALCULATIONS:
        if name not in sections:
            continue
        value = sections[name].get('do', '0')
        if value not in ('0', '1'):
            raise ValueError(f'${name}: do must be 0 or 1')
        normalized[name] = {'do': int(value)}
        for module in CALCULATIONS[name].modules:
            if module in sections[name]:
                selection = sections[name][module]
                if selection not in ('0', '1'):
                    raise ValueError(f'${name}: {module} must be 0 or 1')
                normalized[name][module] = int(selection)
    switches = {}
    for name, settings in sections.items():
        if name in CALCULATIONS:
            continue
        value = settings.get('do')
        if value is not None and value not in ('0', '1'):
            raise ValueError(f'${name}: do must be 0 or 1')
        switches[name] = int(value) if value is not None else None
    route = route_sections(sections, normalized, module_switches=switches)
    active_modules = set(route['tasks'])
    # Independent requests share the same data prerequisites as main workflows.
    extraction_consumers = [name for name in QUANTUM_SECTIONS if name in active_modules
                        and needs_extraction(sections[name])]
    if extraction_consumers and EXTRACTION in sections:
        if switches.get(EXTRACTION) == 0:
            raise ValueError(f'${EXTRACTION} do=0 but selected modules need extracted data; '
                             'enable extraction or provide input/xyz')
        if any('standalone' in route['module_owners'][name] for name in extraction_consumers):
            active_modules.add(EXTRACTION)
    inactive = {name: settings for name, settings in sections.items()
                if name not in CALCULATIONS and name not in active_modules}
    if DOS in active_modules:
        from wren.dos.control import normalize as normalize_dos
        normalized[DOS] = normalize_dos(sections[DOS], path.parent)
    if KMC in active_modules:
        from wren.transport.config import normalize
        normalized[KMC]=normalize(sections[KMC],path.parent)
    if SURFACE_DYNAMICS in active_modules:
        from wren.surface_hopping.control import load_control
        s = sections[SURFACE_DYNAMICS]
        source = resolve_path(s.get('input', ''), path.parent, 'surface-hopping input', required_file=True)
        output = resolve_path(s['output'], path.parent, 'surface-hopping output') if 'output' in s else None
        _, directory, _, _ = load_control(source, output)
        normalized[SURFACE_DYNAMICS] = {'input': str(source), 'output': str(directory)}
    if EXTRACTION in active_modules:
        s = sections[EXTRACTION]
        lammps = {'lammps-data', 'lammps-dump', 'lammps-units', 'lammps-timestep', 'element-map'} & s.keys()
        if lammps and {'tpr', 'gro', 'xtc'} & s.keys():
            raise ValueError('choose either LAMMPS data/dump or GROMACS tpr/gro/xtc for extraction')
        source_files = ('lammps-data', 'lammps-dump', 'element-map') if lammps else ('tpr', 'gro', 'xtc')
        n = {key: str(resolve_path(s.get(key, ''), path.parent, key, required_file=True))
             for key in (*source_files, 'fragment-file')}
        if lammps:
            n['lammps-units'] = s.get('lammps-units', '').lower()
            if n['lammps-units'] not in ('real', 'metal'):
                raise ValueError('specify lammps-units = real or metal')
            if 'lammps-timestep' in s:
                n['lammps-timestep'] = _positive(s['lammps-timestep'], 'lammps-timestep')
        n['frame'] = _integer(s.get('frame') or '-1', 'frame', -1)
        n['frame_step'] = _integer(s.get('frame_step') or '1', 'frame_step')
        n['dimers_cutoff'] = _positive(s.get('dimers_cutoff', '6'), 'dimers_cutoff')
        n['monomers-qmmm-cutoff'] = _positive(s.get('monomers-qmmm-cutoff') or '12', 'monomers-qmmm-cutoff')
        n['output'] = str(resolve_path(s.get('output', 'mol'), path.parent, 'output'))
        n['format'] = [v.lower() for v in _csv(s.get('format', 'xyz'))]
        if len(set(n['format'])) != len(n['format']) or not set(n['format']) <= {'xyz', 'gjf'}:
            raise ValueError('format must select xyz and/or gjf')
        if any(Path(n['output']) == Path(n[k]) for k in (*source_files, 'fragment-file')):
            raise ValueError('output directory must not overwrite an input file')
        normalized[EXTRACTION] = n
    inherited = sections.get('one-mol', {}).get('qc-engine', 'pyscf') if 'one-mol' in active_modules else 'pyscf'
    for name in QUANTUM_SECTIONS:
        if name not in active_modules:
            continue
        s = sections[name]
        pair = None
        adapter = None
        if name == 'reorganization-energy':
            analysis = _normalize_analysis(s, path.parent)
            if analysis is not None:
                normalized[name] = analysis
                continue
            if {'state-a', 'state-b'}.intersection(s):
                from .state_pairs import explicit_pair, quantum_pair_adapter
                if 'carriers' in s:
                    raise ValueError('choose either state-a/state-b or carriers')
                pair = explicit_pair(s)
                adapter = quantum_pair_adapter(pair)
            if adapter != 'tda' and {'root-gap-ev', 'root-min-overlap'}.intersection(s):
                raise ValueError('root diagnostics apply only to the TDA state-pair adapter')
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
             'auxbasis': s.get('auxbasis', 'none' if adapter == 'tda' else 'def2-universal-jkfit')}
        scf_config(n).validate()
        for key in ('input', 'output'):
            if key in s:
                n[key] = str(resolve_path(s[key], path.parent, key, required_file=(key == 'input')))
        if 'xyz' in s:
            n['xyz'] = [str(resolve_path(v, path.parent, 'xyz', required_file=True)) for v in _csv(s['xyz'])]
        separated = bool({'coordinates', 'topology', 'forcefield', 'electrostatics'}.intersection(s))
        if separated:
            if not {'coordinates', 'topology'} <= s.keys() or {'input', 'xyz', 'molecules'}.intersection(s):
                raise ValueError('a separated QM/MM cluster needs coordinates and topology as its sole geometry source')
            if s.get('environment', '').lower() != 'qmmm':
                raise ValueError('a separated QM/MM cluster requires environment=qmmm')
            for key in ('coordinates', 'topology', 'forcefield'):
                if key in s:
                    n[key] = str(resolve_path(s[key], path.parent, key, required_file=True))
            if s.get('electrostatics', 'point-charge').lower() != 'point-charge':
                raise ValueError('separated QM/MM supports electrostatics=point-charge')
        if not (EXTRACTION in normalized or 'input' in n or 'xyz' in n or separated):
            raise ValueError(f'${name} needs ${EXTRACTION}, input = manifest.json, or xyz = ...')
        selection = 'dimers' if name == 'two-mol' else 'molecules'
        if selection in s:
            n[selection] = _csv(s[selection])
        if name == 'two-mol':
            methods = [v.lower() for v in _csv(s.get('transfer-integrals-method', 'DIPRO'))]
            if not set(methods) <= {'dipro', 'esid', 'orbital_gmh'}:
                raise ValueError('transfer-integrals-method must select DIPRO, ESID or orbital_gmh')
            n['transfer-integrals-method'] = methods
            n['assume-equivalent'] = _boolean(s.get('assume-equivalent', 'false'), 'assume-equivalent')
        if name in ('two-mol', 'reorganization-energy') and pair is None:
            n['carriers'] = [v.lower() for v in _csv(s.get('carriers', 'hole,electron'))]
            if not set(n['carriers']) <= {'hole', 'electron'}:
                raise ValueError('carriers must select hole and/or electron')
        if name == 'reorganization-energy':
            n['environment'] = s.get('environment', 'vacuum').lower()
            n['re-method'] = s.get('re-method', '4-points').lower()
            n['opt-maxsteps'] = _integer(s.get('opt-maxsteps', '100'), 'opt-maxsteps')
            if n['environment'] not in ('vacuum', 'qmmm') or n['re-method'] != '4-points':
                raise ValueError('use environment=vacuum or qmmm and re-method=4-points')
            if pair is not None:
                n.update(pair, adapter=adapter)
                if adapter == 'tda':
                    if engine != 'pyscf' or n['environment'] != 'vacuum' or n['auxbasis'].lower() != 'none':
                        raise ValueError('singlet TDA optimization requires qc-engine=pyscf, environment=vacuum and auxbasis=none')
                    n['root-gap-ev'] = _positive(s.get('root-gap-ev', '0.05'), 'root-gap-ev')
                    n['root-min-overlap'] = _positive(s.get('root-min-overlap', '0.7'), 'root-min-overlap')
                    if n['root-min-overlap'] > 1:
                        raise ValueError('root-min-overlap must be in (0,1]')
            mm_keys = {'mm-source', 'mm-charges', 'mm-cutoff', 'mm-frozen', 'mm-vdw', 'mm-vdw-file', 'gmx-executable'}
            if n['environment'] == 'vacuum' and mm_keys.intersection(s):
                raise ValueError('MM settings require environment=qmmm')
            if n['environment'] == 'qmmm':
                n['mm-frozen'] = _boolean(s.get('mm-frozen', 'true'), 'mm-frozen')
                if not n['mm-frozen']:
                    raise ValueError('QM/MM four-point currently requires mm-frozen=true')
                if separated:
                    if {'mm-source', 'mm-cutoff', 'mm-charges', 'mm-vdw-file', 'gmx-executable'}.intersection(s):
                        raise ValueError('the supplied cluster already defines MM membership and files')
                    n['mm-source'] = 'separated'
                    n['mm-vdw'] = s.get('mm-vdw', 'forcefield' if 'forcefield' in n else 'none').lower()
                    if n['mm-vdw'] != ('forcefield' if 'forcefield' in n else 'none'):
                        raise ValueError('mm-vdw=forcefield requires forcefield; omit the file for mm-vdw=none')
                    normalized[name] = n
                    continue
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
                        raise ValueError('mm-source=manifest requires an extraction manifest, without xyz/mm-charges')
                    if 'mm-cutoff' in s:
                        n['mm-cutoff'] = _positive(s['mm-cutoff'], 'mm-cutoff')
                    elif 'input' not in n and EXTRACTION in normalized:
                        n['mm-cutoff'] = normalized[EXTRACTION]['monomers-qmmm-cutoff']
                    # An explicit manifest carries its own extraction cutoff.
                n['mm-vdw'] = s.get('mm-vdw', 'none').lower()
                if n['mm-vdw'] not in ('none', 'tpr', 'file'):
                    raise ValueError('mm-vdw must be none, tpr or file')
                if n['mm-vdw'] == 'tpr':
                    if source != 'manifest' or 'mm-vdw-file' in s:
                        raise ValueError('mm-vdw=tpr requires mm-source=manifest without mm-vdw-file')
                    if 'input' not in n and 'lammps-data' in normalized.get(EXTRACTION, {}):
                        raise ValueError('mm-vdw=tpr requires a GROMACS TPR; LAMMPS force-field coefficients are not imported')
                    executable = s.get('gmx-executable', 'gmx')
                    n['gmx-executable'] = (str(resolve_path(executable, path.parent, 'gmx-executable', required_file=True))
                                           if '/' in executable else executable)
                elif 'gmx-executable' in s:
                    raise ValueError('gmx-executable requires mm-vdw=tpr')
                if n['mm-vdw'] == 'file':
                    if source != 'file':
                        raise ValueError('mm-vdw=file requires mm-charges and one xyz')
                    n['mm-vdw-file'] = str(resolve_path(s.get('mm-vdw-file', ''), path.parent, 'mm-vdw-file', required_file=True))
                elif 'mm-vdw-file' in s:
                    raise ValueError('mm-vdw-file requires mm-vdw=file')
        normalized[name] = n
    if 'two-mol' in normalized and 'one-mol' in normalized:
        for key in ('method', 'basis-set', 'auxbasis', 'grid-level', 'conv-tol'):
            if str(normalized['two-mol'][key]).lower() != str(normalized['one-mol'][key]).lower():
                raise ValueError(f'DIPRO requires matching one-mol/two-mol {key}')
    # Routing flags are reported separately and must not change scientific
    # fingerprints/checkpoints when switching between main and standalone use.
    return WorkflowInput(path, hashlib.sha256(raw).hexdigest(), normalized, inactive, switches)


def scf_config(settings):
    return SCFConfig(backend='gpu' if settings['qc-engine'] == 'gpu-pyscf' else 'cpu',
                     basis=settings['basis-set'], xc=settings['method'], auxbasis=settings['auxbasis'],
                     threads=settings['threads'], grid_level=settings['grid-level'],
                     conv_tol=settings['conv-tol'], max_cycle=settings['max-cycle'],
                     max_memory_mb=settings['max-memory-mb'])
