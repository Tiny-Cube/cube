"""Dependency-free, strict TOML controls for surface-hopping model runs."""
from dataclasses import dataclass, asdict
import hashlib
import json
import math
from pathlib import Path
import tomllib

ALGORITHM_VERSION = 'cube-real-fssh-1'


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            digest.update(block)
    return digest.hexdigest()


def _real(value, name, *, positive=False):
    if type(value) not in (int, float) or not math.isfinite(value) or (positive and value <= 0):
        raise ValueError(f'{name} must be finite'+(' and positive' if positive else ''))
    return float(value)


def _integer(value, name, minimum=1):
    if type(value) is not int or value < minimum:
        raise ValueError(f'{name} must be an integer >= {minimum}')
    return value


def _vector(value, name, *, positive=False):
    if not isinstance(value, list) or not value:
        raise ValueError(f'{name} must be a nonempty array')
    return tuple(_real(v, name, positive=positive) for v in value)


@dataclass(frozen=True)
class DynamicsConfig:
    model: str
    model_file: str | None
    coordinates: tuple[float, ...]
    momenta: tuple[float, ...]
    masses: tuple[float, ...]
    active: int
    dt: float
    substeps: int
    trajectories: int
    steps: int
    seed: int
    decoherence: str
    edc_constant: float
    stride: int
    boundaries: tuple[float, float] | None


def load_control(filename, output=None):
    path = Path(filename).expanduser().resolve()
    raw = tomllib.loads(path.read_text())
    fields = {'model': {'type', 'file'}, 'initial': {'coordinates', 'momenta', 'masses', 'active-state'},
              'dynamics': {'dt-au', 'electronic-substeps', 'trajectories', 'max-steps', 'seed', 'decoherence', 'edc-constant-hartree', 'output-stride'},
              'scattering': {'left', 'right'}, 'output': {'directory'}}
    if set(raw)-fields.keys():
        raise ValueError(f'unknown surface-hopping tables: {sorted(set(raw)-fields.keys())}')
    for table, allowed in fields.items():
        s = raw.get(table, {})
        if not isinstance(s, dict) or set(s)-allowed:
            raise ValueError(f'unknown/invalid fields in [{table}]')
    try:
        model, initial, dynamics = (raw[k] for k in ('model', 'initial', 'dynamics'))
        kind = model['type']
        if kind not in ('tully-1', 'tully-2', 'tully-3', 'lvc'):
            raise ValueError('model.type must be tully-1, tully-2, tully-3 or lvc; atomistic QC is not connected')
        model_path = None
        if kind == 'lvc':
            model_path = (path.parent/model['file']).resolve()
            if not model_path.is_file():
                raise ValueError('LVC model file missing')
        elif 'file' in model:
            raise ValueError('model.file is only valid for lvc')
        q = _vector(initial['coordinates'], 'initial.coordinates')
        p = _vector(initial['momenta'], 'initial.momenta')
        m = _vector(initial['masses'], 'initial.masses', positive=True)
        if len(q) != len(p) or len(q) != len(m) or (kind.startswith('tully') and len(q) != 1):
            raise ValueError('initial vector lengths differ, or Tully model is not one dimensional')
        active = _integer(initial.get('active-state', 0), 'active-state', 0)
        if kind.startswith('tully') and active > 1:
            raise ValueError('Tully active-state must be 0 or 1')
        dt = _real(dynamics.get('dt-au', .5), 'dt-au', positive=True)
        substeps = _integer(dynamics.get('electronic-substeps', 4), 'electronic-substeps')
        trajectories = _integer(dynamics.get('trajectories', 32), 'trajectories')
        steps = _integer(dynamics.get('max-steps', 4000), 'max-steps')
        seed = _integer(dynamics.get('seed', 2026), 'seed', 0)
        decoherence = dynamics.get('decoherence', 'none')
        if decoherence not in ('none', 'edc'):
            raise ValueError('decoherence must be none or edc')
        constant = _real(dynamics.get('edc-constant-hartree', .1), 'edc-constant-hartree', positive=True)
        stride = _integer(dynamics.get('output-stride', 100), 'output-stride')
        boundaries = None
        if 'scattering' in raw:
            boundaries = tuple(_real(raw['scattering'][k], f'scattering.{k}') for k in ('left', 'right'))
            if len(q) != 1 or not boundaries[0] < q[0] < boundaries[1]:
                raise ValueError('scattering needs one coordinate with left < initial q < right')
        config = DynamicsConfig(kind, str(model_path) if model_path else None, q, p, m, active, dt,
                                substeps, trajectories, steps, seed, decoherence, constant, stride, boundaries)
        directory = Path(output).expanduser().resolve() if output is not None else (path.parent/raw.get('output', {}).get('directory', 'results')).resolve()
    except KeyError as exc:
        raise ValueError(f'missing required surface-hopping control field {exc}') from exc
    sources = {str(path): sha256(path)}
    if model_path:
        sources[str(model_path)] = sha256(model_path)
    if any(directory == Path(source) for source in sources):
        raise ValueError('output directory must not overwrite control or model file')
    provenance = {'algorithm': ALGORITHM_VERSION, 'sources': sources, 'config': asdict(config), 'units': 'atomic'}
    key = hashlib.sha256(json.dumps(provenance, sort_keys=True).encode()).hexdigest()
    return config, directory, provenance, key
