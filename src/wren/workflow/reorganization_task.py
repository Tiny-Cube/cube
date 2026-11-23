"""Route reusable reorganization algorithms and project their results into workflow reports."""
import hashlib
import json
import math
import os
from pathlib import Path
import tempfile

from wren.quantum.cli import _write_result
from wren.reorganization import ElectronicState


def analyze_existing(settings, *, resume=False):
    from wren import reorganization
    from .reorganization_methods import ANALYSIS_METHODS
    method, record = ANALYSIS_METHODS[settings['re-method']].analyze(settings)
    count = len(record['pairs']) if 'pairs' in record else 1
    record['implementation_sha256'] = {
        p.name: hashlib.sha256(p.read_bytes()).hexdigest()
        for p in map(Path, (method.__file__, reorganization.__file__))
    }
    directory = Path(settings['output'])
    output = directory/'reorganization.json'
    directory.mkdir(parents=True, exist_ok=True)
    # The entire analysis is deterministic and cheap enough to regenerate. An
    # equality check covers the source hashes, numerical values and code version.
    reused = output.exists()
    if reused:
        if not resume or json.loads(output.read_bytes()) != record:
            raise ValueError(f'existing analysis differs or --resume missing: {output}')
    else:
        with tempfile.TemporaryDirectory(dir=directory) as scratch:
            temporary = Path(scratch)/output.name
            _write_result(temporary, record)
            os.link(temporary, output)
    return {'status': 'completed', 'method': settings['re-method'], 'energy_source': 'external',
            'environment': record['environment'], 'output': str(output), 'count': count, 'reused': reused}


def _environment_cutoff(settings, data, record):
    if 'mm-cutoff' in settings:
        return settings['mm-cutoff']
    extracted = data.get('provenance', {}).get('settings', {}).get('monomers-qmmm-cutoff', 12.)
    value = record.get('cluster', {}).get('cutoff_angstrom', extracted)
    if type(value) not in (int, float) or not math.isfinite(value) or value <= 0:
        raise ValueError('invalid stored monomer environment cutoff; specify mm-cutoff in Angstrom')
    return value


def run_quantum_reorganization(settings, data, records, atoms_of, *, resume=False):
    from .reorganization_csv import write_reorganization_csv, UNITS
    directory = Path(settings.get('output', data['root']/'reorganization'))
    frame = fixed_environment = pair_parameters = None
    if settings['environment'] == 'qmmm':
        from .embedding import read_charges, load_manifest_environment
        if settings['mm-source'] == 'separated':
            fixed_environment = data['frozen_environment']
        elif settings['mm-source'] == 'file':
            fixed_environment = read_charges(settings['mm-charges'])
        else:
            if settings.get('mm-vdw') == 'tpr' and 'lammps-data' in data.get('provenance', {}).get('settings', {}):
                raise ValueError('mm-vdw=tpr requires a GROMACS TPR; LAMMPS force-field coefficients are not imported')
            frame = load_manifest_environment(data)
        if settings.get('mm-vdw') == 'tpr':
            from .lj import read_gromacs_lj
            pair_parameters = read_gromacs_lj(data['provenance']['settings']['tpr'], settings['gmx-executable'])

    adapter = settings.get('adapter')
    if adapter:
        states = {name: ElectronicState(**settings[key]) for name, key in (
            ('state_a','state-a'), ('state_b','state-b'))}
        if adapter == 'tda':
            from .excited_reorganization import calculate_excited_reorganization as calculate
        else:
            from .scf_pair import calculate_scf_pair as calculate
    else:
        from .reorganization import calculate_reorganization as calculate
        states = {}

    outputs, csv_results = {}, {}
    report = {'status': 'completed', 'environment': settings['environment'], 'count': len(records),
              'outputs': outputs}
    for record in records:
        atoms = atoms_of(record)
        environment = fixed_environment
        if frame is not None:
            from .embedding import cluster_from_frame
            environment = cluster_from_frame(frame, record, _environment_cutoff(settings, data, record))
        if settings.get('mm-vdw') in ('tpr', 'file'):
            from dataclasses import replace
            from .lj import bind_gromacs_lj, read_lj_file
            vdw = (bind_gromacs_lj(pair_parameters, record, environment) if pair_parameters is not None else
                   read_lj_file(settings['mm-vdw-file'], atoms, environment))
            environment = replace(environment, vdw=vdw)
        identity = record['id']
        print(f'wren: {settings["environment"]} four-point reorganization {identity}', flush=True)
        kwargs = {} if adapter == 'tda' else {'embedding': environment}
        result = calculate(atoms, settings, directory/identity, resume=resume, **states, **kwargs)
        item = {'output': str(directory/identity/'reorganization.json'), 'reused': result['reused']}
        if adapter:
            item['relaxation_ev'] = result['relaxation_ev']
        else:
            item['reorganization_ev'] = {c: v['reorganization_energy_ev'] for c, v in result['results'].items()}
            csv_results[identity] = result
            report['csv'] = write_reorganization_csv(csv_results, directory, settings['carriers'])
            report['csv_units'] = UNITS
        outputs[identity] = item
    if not records:
        raise ValueError('no monomers selected for reorganization energy')
    return report
