"""Optimize any two distinct SCF charge/spin sectors in the same environment."""
from dataclasses import asdict
import json
from pathlib import Path

from wren import reorganization as core
from wren.quantum.cli import _write_result
from wren.quantum.transfer import _versions
from . import embedding as embedding_module, reorganization, scf
from .control import scf_config
from .dimers import file_hash, fingerprint


def calculate_scf_pair(atoms, settings, directory, *, state_a, state_b, resume=False, embedding=None):
    from pyscf.data import elements
    from pyscf.data.nist import HARTREE2EV
    states = {'a': state_a, 'b': state_b}
    nuclear_charge = sum(elements.charge(symbol) for symbol, _ in atoms)
    config = scf_config(settings)
    for state in states.values():
        electrons = nuclear_charge - state.charge
        unpaired = state.multiplicity - 1
        if electrons < 1 or electrons < unpaired or (electrons - unpaired) % 2:
            raise ValueError(f'{state.name}: QM electron count and multiplicity are incompatible')
    if (settings['environment'] == 'qmmm') != (embedding is not None):
        raise ValueError('SCF state-pair embedding must match the requested environment')
    directory = Path(directory)
    output = directory/'reorganization.json'
    request = {'states': {k: asdict(v) for k, v in states.items()}, 'atoms': atoms,
               'scf': asdict(config), 'opt_maxsteps': settings['opt-maxsteps'], 'versions': _versions(),
               'embedding': embedding.fingerprint if embedding else None,
               'implementation_sha256': {Path(p).name: file_hash(p) for p in (
                   __file__, scf.__file__, reorganization.__file__, embedding_module.__file__)},
               'algebra_sha256': file_hash(core.__file__)}
    key = fingerprint(request)
    if output.exists():
        old = json.loads(output.read_text())
        if not resume or old.get('fingerprint') != key:
            raise ValueError(f'existing state-pair task {directory}; use --resume with identical inputs')
        if old.get('status') != 'completed':
            raise ValueError('existing incomplete SCF state-pair job; use a fresh output directory')
        for name, digest in old['output_sha256'].items():
            if not (directory/name).is_file() or file_hash(directory/name) != digest:
                raise ValueError(f'cannot resume: state-pair output changed: {name}')
        return {**old, 'reused': True}
    elif directory.exists() and any(directory.iterdir()):
        raise ValueError(f'untracked state-pair outputs in {directory}')
    record = {'schema': 'wren-scf-pair-1', 'method': '4-points', 'status': 'running', 'fingerprint': key,
              'input': request, 'environment': settings['environment'], 'reused': False}
    directory.mkdir(parents=True, exist_ok=True)
    with output.open('x') as initial:
        json.dump(record, initial)
    try:
        if embedding is not None:
            embedding_module.snapshot(embedding, directory/'mm-environment.json')
        geometries, optimization = {}, {}
        for k, state in states.items():
            geometries[k], optimization[k] = reorganization._optimize(
                geometries.get('a', atoms), config, state.charge, directory, settings['opt-maxsteps'],
                spin=state.multiplicity - 1, label='state_'+k, resume=resume, embedding=embedding)
        points = {}
        for sector, geometry in (('a','a'), ('b','a'), ('b','b'), ('a','b')):
            state = states[sector]
            name = f'E_{sector.upper()}_R{geometry.upper()}'
            points[name] = scf.run_single_point(geometries[geometry], settings, directory/name,
                charge=state.charge, spin=state.multiplicity - 1, embedding=embedding, resume=resume)
        energies = core.FourPointEnergies(*(p['energy_hartree'] for p in points.values()), unit='hartree')
        analysis = core.analyze_four_point_pair(state_a, state_b, energies)
        relaxation_ev = {k: value * HARTREE2EV for k, value in analysis['relaxation'].items()}
        if min(relaxation_ev['on_a'], relaxation_ev['on_b']) < -1e-5:
            raise RuntimeError('negative SCF surface deformation; inspect geometries and SCF solutions')
        record.update(status='completed', analysis=analysis, relaxation_ev=relaxation_ev, points=points,
                      optimizations=optimization, output_sha256={p.name: file_hash(p) for p in directory.iterdir()
                                                                if p.is_file() and p != output})
    except Exception as exc:
        _write_result(output, {**record, 'status': 'failed', 'error': str(exc)})
        raise
    _write_result(output, record)
    return record
