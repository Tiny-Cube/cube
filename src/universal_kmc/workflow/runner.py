"""Dispatch input sections to extraction, single-point, coupling and RE tasks."""
import csv
from dataclasses import asdict
import json
from pathlib import Path
import re

from universal_kmc.quantum.cli import _write_result
from universal_kmc.quantum.geometry import normalize_atoms, read_xyz
from .control import EXTRACTION, scf_config
from .dimers import extract_dimers_monomers, file_hash, fingerprint


def _atoms(record):
    if len(record['elements']) != len(record['coordinates_angstrom']):
        raise ValueError('manifest atom/coordinate counts differ')
    return normalize_atoms(zip(record['elements'], record['coordinates_angstrom']))


def _select(records, selection, label):
    by_id = {str(r['id']): r for r in records}
    if len(by_id) != len(records) or any(not re.fullmatch(r'[A-Za-z0-9_-]+', key) for key in by_id):
        raise ValueError('manifest IDs must be unique safe names')
    if selection is None:
        return list(records)
    missing = set(selection)-by_id.keys()
    if missing:
        raise ValueError(f'unknown {label}: {", ".join(sorted(missing))}')
    return [by_id[key] for key in selection]


def _load_dataset(settings, extracted, control_dir):
    if 'xyz' in settings:
        monomers = []
        for name in settings['xyz']:
            atoms = read_xyz(name)
            identity = re.sub('[^A-Za-z0-9_-]', '_', Path(name).stem)
            monomers.append({'id': identity, 'elements': [a[0] for a in atoms],
                             'coordinates_angstrom': [a[1] for a in atoms]})
        _select(monomers, None, 'molecules')
        return {'monomers': monomers, 'dimers': [], 'root': control_dir/'cube-results'}
    source = Path(settings['input']) if 'input' in settings else extracted
    if source is None:
        raise ValueError('a quantum task needs an extraction manifest or xyz inputs')
    data = json.loads(source.read_text())
    if data.get('schema') != 'cube-dimers-1' or data.get('coordinate_unit') != 'Angstrom':
        raise ValueError('input must refer to a cube-dimers-1 Angstrom manifest')
    _select(data['monomers'], None, 'molecules')
    _select(data['dimers'], None, 'dimers')
    data['root'] = source.parent
    return data


def _dataset_sources(settings, extracted):
    if 'xyz' in settings:
        return [None]
    source = Path(settings['input']) if 'input' in settings else extracted
    if source is None:
        raise ValueError('a quantum task needs an extraction manifest or xyz inputs')
    index = json.loads(source.read_text())
    if index.get('schema') != 'cube-frames-1':
        return [source]
    entries = index.get('frames', [])
    if not entries or index.get('frame_count') != len(entries):
        raise ValueError('empty or inconsistent frame collection')
    sources, seen = [], set()
    for entry in entries:
        frame = entry.get('frame')
        if type(frame) is not int or frame < 0 or frame in seen:
            raise ValueError('frame collection needs unique nonnegative frame indices')
        seen.add(frame)
        relative = Path(entry['manifest'])
        path = (source.parent/relative).resolve()
        if relative.is_absolute() or not path.is_relative_to(source.parent.resolve()):
            raise ValueError('frame manifest must be within its collection directory')
        if not path.is_file() or file_hash(path) != entry['sha256']:
            raise ValueError('frame collection manifest missing/changed')
        if json.loads(path.read_text()).get('frame') != frame:
            raise ValueError('frame collection index and manifest disagree')
        sources.append(path)
    return sources


def run_workflow(control, *, only=None, resume=False):
    tasks = control.plan(only)['tasks']
    sections = control.sections
    extracted = None
    report = {'schema': 'cube-workflow-1', 'status': 'running', 'control': str(control.path),
              'control_sha256': control.sha256, 'tasks': {}, 'task_order': tasks}
    summary_path = None
    monomer_cache = {}

    def single(record, settings, prefix):
        from .scf import run_single_point
        atoms = _atoms(record)
        key = fingerprint({'prefix': str(prefix), 'atoms': atoms, 'config': asdict(scf_config(settings))})
        if key not in monomer_cache:
            print(f'cube: SCF {prefix}', flush=True)
            monomer_cache[key] = run_single_point(atoms, settings, prefix, resume=resume)
        return monomer_cache[key]

    try:
        if EXTRACTION in sections and (EXTRACTION in tasks or any('input' not in sections[t] and 'xyz' not in sections[t] for t in tasks if t != EXTRACTION)):
            print('cube: extracting complete molecules/dimers', flush=True)
            extracted = extract_dimers_monomers(sections[EXTRACTION], one_settings=sections.get('one-mol'),
                                       two_settings=sections.get('two-mol'), resume=resume)
            summary_path = extracted.parent/'workflow.json'
            manifest = json.loads(extracted.read_text())
            if manifest['schema'] == 'cube-frames-1':
                report['tasks'][EXTRACTION] = {'status': 'completed', 'manifest': str(extracted),
                    'frame_count': manifest['frame_count'], 'frames': manifest['frames']}
            else:
                report['tasks'][EXTRACTION] = {'status': 'completed', 'manifest': str(extracted),
                    'frame': manifest['frame'], 'molecule_count': manifest['molecule_count'], 'dimer_count': manifest['dimer_count']}
            _write_result(summary_path, report)
        for task in tasks:
            if task == EXTRACTION:
                continue
            original_settings = sections[task]
            sources = _dataset_sources(original_settings, extracted)
            multi = len(sources) > 1
            results = {}
            if summary_path is None:
                summary_path = control.path.parent/'cube-results/workflow.json'
            for source in sources:
                settings = dict(original_settings)
                if source is not None:
                    settings['input'] = str(source)
                data = _load_dataset(settings, source, control.path.parent)
                frame_label = f"frame{data['frame']:04d}" if 'frame' in data else 'dataset'
                if multi and 'output' in settings:
                    settings['output'] = str(Path(settings['output'])/frame_label)
                if task == 'one-mol':
                    directory = Path(settings.get('output', data['root']/'qc/one-mol'))
                    records = _select(data['monomers'], settings.get('molecules'), 'molecules')
                    outputs = {r['id']: single(r, settings, directory/r['id']) for r in records}
                    task_result = {'status': 'completed', 'count': len(outputs), 'outputs': outputs}
                elif task == 'two-mol':
                    from universal_kmc.quantum.checkpoints import calculate_from_checkpoints
                    from universal_kmc.quantum.control import TransferControl
                    from .scf import translated_checkpoint
                    directory = Path(settings.get('output', data['root']/'qc/two-mol'))
                    mono_settings = dict(sections.get('one-mol', settings))
                    if multi and 'one-mol' in sections and 'output' in mono_settings:
                        mono_settings['output'] = str(Path(mono_settings['output'])/frame_label)
                    mono_directory = Path(mono_settings.get('output', data['root']/'qc/one-mol'))
                    monomers = {r['id']: r for r in data['monomers']}
                    pairs = _select(data['dimers'], settings.get('dimers'), 'dimers')
                    outputs, csv_rows = {}, []
                    for pair in pairs:
                        a_id, b_id = str(pair['mol_i']), str(pair['mol_j'])
                        if a_id not in monomers or b_id not in monomers:
                            raise ValueError('dimer refers to an absent monomer')
                        a = single(monomers[a_id], mono_settings, mono_directory/a_id)
                        b = single(monomers[b_id], mono_settings, mono_directory/b_id)
                        pair_dir = directory/pair['id']
                        dimer = single(pair, settings, pair_dir/'ab')
                        atoms = _atoms(pair)
                        n_a = pair['atom_count_a']
                        if type(n_a) is not int or n_a != len(_atoms(monomers[a_id])) or len(atoms)-n_a != len(_atoms(monomers[b_id])):
                            raise ValueError('dimer/monomer atom counts do not match')
                        output = pair_dir/'transfer.json'
                        task_key = fingerprint({'a': a['checkpoint_sha256'], 'b': b['checkpoint_sha256'],
                                                'dimer': dimer['checkpoint_sha256'], 'settings': settings})
                        result = None
                        if output.exists():
                            previous = json.loads(output.read_text())
                            if not resume or previous.get('workflow_fingerprint') != task_key:
                                raise ValueError(f'existing transfer task {output}; use --resume with identical inputs')
                            if previous['status'] == 'completed':
                                for source in previous['input']['checkpoints'].values():
                                    if not Path(source['path']).is_file() or file_hash(source['path']) != source['sha256']:
                                        raise ValueError('cannot resume: coupling checkpoint missing/changed')
                                result = previous
                        was_reused = result is not None
                        if result is None:
                            a_path = translated_checkpoint(a['checkpoint'], atoms[:n_a], pair_dir/'a.chk')
                            b_path = translated_checkpoint(b['checkpoint'], atoms[n_a:], pair_dir/'b.chk')
                            provenance = {'schema_version': 1, 'format': 'cube-sections-1', 'path': str(control.path),
                                          'sha256': control.sha256, 'dimer': pair['id']}
                            job = TransferControl({'fragment_a': a_path, 'fragment_b': b_path, 'dimer': Path(dimer['checkpoint'])},
                                scf_config(settings), tuple(settings['transfer-integrals-method']), tuple(settings['carriers']),
                                settings['assume-equivalent'], False, None, None, output, provenance)
                            result = calculate_from_checkpoints(job)
                            result['workflow_fingerprint'] = task_key
                            _write_result(output, result)
                        outputs[pair['id']] = {'status': result['status'], 'output': str(output.resolve()), 'reused': was_reused}
                        row = {'dimer': pair['id'], 'mol_i': a_id, 'mol_j': b_id,
                               'core_min_angstrom': pair.get('core_min_angstrom'), 'result_json': str(output.resolve())}
                        method = settings['transfer-integrals-method'][0]
                        for carrier, column in (('hole', 'VH_meV'), ('electron', 'VL_meV')):
                            value = result['results'].get(carrier, {}).get('methods', {}).get(method, {})
                            block = value.get('coupling_mev')
                            row[column] = (block[0][0] if block and len(block) == len(block[0]) == 1 else
                                           value.get('scalar_magnitude_mev'))
                        csv_rows.append(row)
                    directory.mkdir(parents=True, exist_ok=True)
                    with (directory/'transfer_integrals.csv').open('w', newline='') as file:
                        writer = csv.DictWriter(file, fieldnames=['dimer','mol_i','mol_j','core_min_angstrom','VH_meV','VL_meV','result_json'])
                        writer.writeheader()
                        writer.writerows(csv_rows)
                    complete = all(v['status'] == 'completed' for v in outputs.values())
                    task_result = {'status': 'completed' if complete else 'partial', 'count': len(outputs),
                                             'csv': str(directory/'transfer_integrals.csv'), 'outputs': outputs}
                else:
                    from .reorganization import calculate_reorganization
                    directory = Path(settings.get('output', data['root']/'reorganization'))
                    records = _select(data['monomers'], settings.get('molecules'), 'molecules')
                    outputs = {}
                    frame = fixed_environment = pair_parameters = None
                    if settings['environment'] == 'qmmm':
                        from .embedding import read_charges, load_gromacs_environment, cluster_from_frame
                        if settings['mm-source'] == 'file':
                            fixed_environment = read_charges(settings['mm-charges'])
                        else:
                            frame = load_gromacs_environment(data)
                        if settings.get('mm-vdw') == 'tpr':
                            from .lj import read_gromacs_lj
                            pair_parameters = read_gromacs_lj(data['provenance']['settings']['tpr'], settings['gmx-executable'])
                    for record in records:
                        environment = (cluster_from_frame(frame, record, settings['mm-cutoff'])
                                       if frame is not None else fixed_environment)
                        if settings.get('mm-vdw') in ('tpr', 'file'):
                            from dataclasses import replace
                            from .lj import bind_gromacs_lj, read_lj_file
                            vdw = (bind_gromacs_lj(pair_parameters, record, environment) if pair_parameters is not None else
                                   read_lj_file(settings['mm-vdw-file'], _atoms(record), environment))
                            environment = replace(environment, vdw=vdw)
                        print(f'cube: {settings["environment"]} four-point reorganization {record["id"]}', flush=True)
                        result = calculate_reorganization(_atoms(record), settings, directory/record['id'], resume=resume,
                                                          embedding=environment)
                        outputs[record['id']] = {'output': str(directory/record['id']/'reorganization.json'),
                            'reorganization_ev': {c:v['reorganization_energy_ev'] for c,v in result['results'].items()},
                            'reused': result['reused']}
                    task_result = {'status': 'completed', 'environment': settings['environment'], 'count': len(outputs), 'outputs': outputs}
                results[frame_label] = task_result
                report['tasks'][task] = (task_result if not multi else {
                    'status': ('running' if len(results) < len(sources) else
                               'completed' if all(r['status'] == 'completed' for r in results.values()) else 'partial'),
                    'frame_count': len(results), 'count': sum(r['count'] for r in results.values()), 'frames': results})
                _write_result(summary_path, report)
        complete = all(t['status'] == 'completed' for t in report['tasks'].values())
        report['status'] = 'completed' if complete else 'partial'
        _write_result(summary_path, report)
        print(f'cube: {report["status"]}; report: {summary_path}', flush=True)
        return 0 if complete else 2
    except Exception as exc:
        if summary_path is not None:
            report.update(status='failed', error={'type': type(exc).__name__, 'message': str(exc)})
            _write_result(summary_path, report)
        raise
