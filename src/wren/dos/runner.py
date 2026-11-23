"""MD manifests/XYZ -> PySCF spectra -> normalized electronic DOS.

Monomer spectra omit intermolecular mixing. Periodic spectra use full-cell
restricted PySCF k-point SCF. Tight binding requires explicit signed couplings.
"""
from dataclasses import asdict
import csv
import json
from pathlib import Path
import time

import numpy as np

from wren.quantum.cli import _write_result
from wren.quantum.config import SCFConfig
from wren.quantum.geometry import read_xyz
from wren.workflow.dimers import file_hash, fingerprint
from .core import broaden, bloch_eigensystems, disorder_statistics, element_probabilities, kmesh, participation_groups


def _write_total_dos(directory, settings, grid, spectra, weights, frames):
    curves = [broaden(e, grid, settings['sigma'], weights=w, backend=settings['dos-backend'])
              for e, w in zip(spectra, weights, strict=True)]
    mean = np.mean(curves, axis=0)
    with (directory/'dos.csv').open('w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['energy_ev', 'dos', *[f'frame_{i}' for i in range(len(curves))]])
        writer.writerows(zip(grid, mean, *curves))
    expected = float(np.mean([w.sum() for w in weights]))
    integral = float(np.trapezoid(mean, grid))
    summary = {'expected_integral': expected, 'finite_grid_integral': integral,
               'integral_relative_error': abs(integral-expected)/max(expected, 1e-300)}
    if 'block-size-frames' in settings:
        from .sampling import block_curve_statistics
        indices, times = _sampling_axis(frames)
        blocks = block_curve_statistics(curves, indices, settings['block-size-frames'], times_ps=times)
        with (directory/'dos-blocks.csv').open('w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(['energy_ev', 'complete_blocks_mean', 'standard_error',
                             *[f'block_{i}' for i in range(len(blocks['block_means']))]])
            writer.writerows(zip(grid, blocks['retained_mean'], blocks['standard_error'], *blocks['block_means']))
        summary['sampling_uncertainty_estimated'] = True
        summary['sampling'] = {**blocks['metadata'], 'file': 'dos-blocks.csv'}
    return summary, mean


def _sampling_axis(frames):
    times = [frame.get('time_ps') for frame in frames]
    return [frame.get('frame') for frame in frames], times if any(t is not None for t in times) else None


def _write_pdos(directory, settings, grid, labels, spectra, projected_weights, total):
    """Apply the same frame average and finite energy window as the total DOS."""
    curves = np.zeros((len(grid), len(labels)))
    counts = np.zeros(len(labels))
    for energies, weights in zip(spectra, projected_weights, strict=True):
        curves += broaden(energies, grid, settings['sigma'], weights=weights, backend=settings['dos-backend'])
        counts += weights.sum(axis=0)
    curves /= len(spectra)
    counts /= len(spectra)
    residual = float(np.max(np.abs(curves.sum(axis=1)-total)))
    if not np.allclose(curves.sum(axis=1), total, rtol=1e-8, atol=1e-12):
        raise ValueError('sum of projected curves disagrees with total DOS')
    with (directory/'pdos.csv').open('w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['energy_ev', *labels])
        writer.writerows((energy, *row) for energy, row in zip(grid, curves))
    return counts.tolist(), np.trapezoid(curves, grid, axis=0).tolist(), residual


def _write_site_dos(directory, settings, spectra, grid, site_weights, total, nonorthogonal):
    labels = [f'site_{i}' for i in range(site_weights[0].shape[1])]
    _, integrals, residual = _write_pdos(directory, settings, grid, labels, spectra,
                                         [np.concatenate(site_weights)], total)
    return {'pdos': 'pdos.csv', 'participation': 'participation.csv',
            'basis': 'symmetric_lowdin_sites_within_model_cell' if nonorthogonal else 'orthonormal_sites_within_model_cell',
            'site_index_base': 0,
            'degeneracy_tol_ev': settings['degeneracy-tol-ev'],
            'degeneracy_grouping': 'consecutive_gaps_at_or_below_tolerance_within_each_k',
            'ipr_convention': 'sum_squared_normalized_site_probabilities',
            'degenerate_group_statistic': 'sum_i (P_ii/rank)^2',
            'expected_integral_per_site': settings['spin-degeneracy'],
            'finite_grid_integral_per_site': integrals,
            'max_abs_site_sum_residual_per_ev': residual}


def _write_element_dos(directory, settings, spectra, grid, frames, weights, total):
    labels = sorted({symbol for frame in frames for symbols, _ in frame for symbol in symbols})
    def frame_weights():
        for frame, line_weights in zip(frames, weights, strict=True):
            blocks = []
            for symbols, populations in frame:
                block = np.zeros((len(populations), len(labels)))
                block[:, [labels.index(symbol) for symbol in symbols]] = populations
                blocks.append(block)
            yield np.vstack(blocks) * line_weights[:, None]
    counts, integrals, residual = _write_pdos(directory, settings, grid, labels, spectra, frame_weights(), total)
    return {'pdos': 'pdos.csv', 'basis': 'symmetric_lowdin_atomic_orbitals', 'elements': labels,
            'expected_integral_per_element': dict(zip(labels, counts)),
            'finite_grid_integral_per_element': dict(zip(labels, integrals)),
            'max_abs_element_sum_residual_per_ev': residual}


def _network_terms(network):
    schema = network.get('schema')
    if schema == 'wren-dos-tb-1':
        if network.get('coupling_convention') != 'signed_real_consistent_orbital_phases':
            raise ValueError('orthonormal model requires an explicit signed-coupling convention')
        if {'overlap_pairs', 'hamiltonian_pairs_ev', 'matrix_element_convention'} & network.keys():
            raise ValueError('raw Hamiltonian and overlap pairs require schema wren-dos-tb-2')
        return network['pairs'], None
    if schema == 'wren-dos-tb-2':
        if network.get('matrix_element_convention') != 'raw_signed_H_and_S_in_same_orbital_gauge':
            raise ValueError('nonorthogonal model requires raw signed H and S in the same orbital gauge')
        if 'pairs' in network or 'coupling_convention' in network:
            raise ValueError('use hamiltonian_pairs_ev for the raw nonorthogonal Hamiltonian')
        if not all(isinstance(network.get(k), list) for k in ('hamiltonian_pairs_ev', 'overlap_pairs')):
            raise ValueError('nonorthogonal model requires hamiltonian_pairs_ev and overlap_pairs lists')
        return network['hamiltonian_pairs_ev'], network['overlap_pairs']
    raise ValueError('unknown tight-binding schema; expected wren-dos-tb-1 or wren-dos-tb-2')


def _dataset(settings):
    if 'xyz' in settings:
        atoms = read_xyz(settings['xyz'])
        return [{'monomers': [{'id': Path(settings['xyz']).stem,
                              'elements': [e for e, r in atoms],
                              'coordinates_angstrom': [r for e, r in atoms]}]}]
    from wren.workflow.runner import _dataset_sources
    sources = _dataset_sources({'input': settings['input']}, None)
    result = []
    for source in sources:
        data = json.loads(source.read_text())
        if data.get('schema') != 'cube-dimers-1' or data.get('coordinate_unit') != 'Angstrom':
            raise ValueError('DOS requires an Angstrom extraction manifest')
        data['_source_sha256'] = file_hash(source)
        result.append(data)
    return result


def _molecules(data, selected):
    records = data['monomers']
    ids = [str(r['id']) for r in records]
    if any(not name or name in ('.', '..') or Path(name).name != name for name in ids):
        raise ValueError('DOS molecule IDs must be nonempty single path components')
    if len(set(ids)) != len(ids):
        raise ValueError('duplicate DOS molecule IDs')
    if selected:
        if not set(selected) <= set(ids):
            raise ValueError('unknown DOS molecule selection')
        records = [records[ids.index(name)] for name in selected]
    return records


def _scf_settings(settings):
    cfg = SCFConfig(backend='gpu' if settings['qc-engine'] == 'gpu-pyscf' else 'cpu',
                    basis=settings['basis-set'], xc=settings['method'],
                    threads=settings['threads'], grid_level=settings['grid-level'],
                    conv_tol=settings['conv-tol'], max_cycle=settings['max-cycle'],
                    max_memory_mb=settings['max-memory-mb'], auxbasis=settings['auxbasis'])
    cfg.validate()
    return cfg


def _complete(directory, result, started):
    output = directory/'result.json'
    result.setdefault('sampling_uncertainty_estimated', False)
    result.update(status='completed',
                  elapsed_seconds=time.perf_counter()-started,
                  output_sha256={str(p.relative_to(directory)): file_hash(p)
                                 for p in directory.rglob('*') if p.is_file() and p != output})
    _write_result(output, result)
    return {**result, 'output': str(output), 'reused': False}


def run(settings, *, resume=False):
    from . import core
    if 'kmesh-study' in settings:
        from .mesh_study import run_study
        single = settings.copy()
        meshes = single.pop('kmesh-study')
        return run_study(single, meshes, resume=resume)
    started = time.perf_counter()
    mode = settings['mode']
    use_kpm = settings.get('solver') == 'kpm'
    datasets = _dataset(settings) if mode != 'tight-binding' else []
    network = json.loads(Path(settings['input']).read_text()) if mode == 'tight-binding' else None
    sources = {key: file_hash(settings[key]) for key in ('input', 'xyz', 'reference-file') if key in settings}
    request = {'settings': settings, 'datasets': datasets, 'network': network, 'source_sha256': sources,
               'implementation_sha256': {p: file_hash(p) for p in (__file__, core.__file__)},
               'numpy_version': np.__version__}
    if 'block-size-frames' in settings:
        from . import sampling
        indices, times = _sampling_axis(datasets)
        sampling.validate_sampling(indices, settings['block-size-frames'], times_ps=times)
        request['implementation_sha256'][sampling.__file__] = file_hash(sampling.__file__)
    if use_kpm:
        from . import kpm
        import scipy
        del request['network']  # The input hash identifies large sparse networks without duplicating them.
        request['implementation_sha256'][kpm.__file__] = file_hash(kpm.__file__)
        request['scipy_version'] = scipy.__version__
        if settings['dos-backend'] == 'gpu':
            import cupy
            request['cupy_version'] = cupy.__version__
    if mode != 'tight-binding':
        import pyscf
        request['pyscf_version'] = pyscf.__version__
        from wren.workflow import scf as scf_module
        from wren.quantum import checkpoints
        request['implementation_sha256'].update({p: file_hash(p) for p in (scf_module.__file__, checkpoints.__file__)})
    if mode == 'periodic':
        from . import alignment, periodic
        references = alignment.reference_values(settings, datasets)
        from wren.quantum.transfer import _versions
        for module in (periodic, alignment):
            request['implementation_sha256'][module.__file__] = file_hash(module.__file__)
        request['versions'] = _versions()
    if mode == 'vertical-dscf':
        from . import charging
        from wren.quantum import geometry
        from wren.workflow import control as workflow_control
        from wren.quantum.transfer import _versions
        request['implementation_sha256'].update({p: file_hash(p) for p in
                                                (charging.__file__, geometry.__file__, workflow_control.__file__)})
        request['versions'] = _versions()
    key = fingerprint(request)
    directory = Path(settings['output'])
    output = directory/'result.json'
    if output.exists():
        old = json.loads(output.read_text())
        if not resume or old.get('fingerprint') != key or old.get('status') != 'completed':
            raise ValueError('existing DOS task; resume requires identical completed inputs')
        for relative, expected in old['output_sha256'].items():
            if not (directory/relative).is_file() or file_hash(directory/relative) != expected:
                raise ValueError('DOS output missing or changed')
        return {**old, 'reused': True, 'output': str(output)}
    if directory.exists() and any(directory.iterdir()):
        raise ValueError('nonempty DOS output directory; choose a new directory')
    directory.mkdir(parents=True, exist_ok=True)
    result = {'schema': 'wren-dos-1', 'status': 'running', 'fingerprint': key, 'request': request,
              'mode': mode, 'energy_unit': 'eV', 'energy_reference': 'common_unshifted_input_zero',
              'backend': settings['dos-backend'], 'frames': []}
    if mode == 'periodic' and 'periodic-reference' in settings:
        reference = settings['periodic-reference']
        result['energy_reference'] = {'frame-vbm': 'per_frame_vbm_shape_only',
                                      'cell-mean-local': 'per_cell_mean_hartree_plus_local_gth',
                                      'external': 'external_frame_reference'}[reference]
        result['alignment'] = {'mode': reference, 'energy_axis': 'raw_orbital_energy_ev - reference_ev',
                               'label': settings.get('reference-label', reference)}
    if not use_kpm:
        result['broadening'] = {'kernel': 'Gaussian', 'sigma_ev': settings['sigma']}
    _write_result(output, result)
    spectra, weights, frame_records, homo_table, lumo_table = [], [], [], [], []
    site_weights, element_frames = [], []
    try:
        if mode == 'tight-binding':
            pairs, overlaps = _network_terms(network)
            if use_kpm:
                if overlaps:
                    raise ValueError('sparse KPM requires an orthonormal Hamiltonian; use diagonalization for overlap pairs')
                result.update(kpm.calculate(directory, network, pairs, settings))
                return _complete(directory, result, started)
            points, kw = kmesh(settings['kmesh'])
            systems = bloch_eigensystems(network['site_energies_ev'], pairs, points, overlap_pairs=overlaps,
                                         backend=settings['dos-backend'],
                                         project_sites=settings.get('site-projections', False))
            bands, groups = [], []
            for index, (energies, probabilities) in enumerate(systems):
                bands.append(energies)
                if probabilities is not None:
                    site_weights.append(probabilities * (kw[index] * settings['spin-degeneracy']))
                    groups.extend({'k_index': index, **group} for group in participation_groups(
                        energies, probabilities, settings['degeneracy-tol-ev']))
            if groups:
                with (directory/'participation.csv').open('w', newline='') as file:
                    writer = csv.DictWriter(file, fieldnames=list(groups[0]))
                    writer.writeheader()
                    writer.writerows(groups)
            values = np.asarray(bands)
            spectra.append(values.ravel())
            weights.append(np.repeat(kw*settings['spin-degeneracy'], values.shape[1]))
            frame_records.append({'eigenvalues_ev': values.tolist(), 'kpoints_fractional': points.tolist(),
                                  'kpoint_weights': kw.tolist()})
            result.update(scope='orthonormal_effective_Hamiltonian' if overlaps is None else 'nonorthogonal_effective_Hamiltonian',
                          normalization='states_per_model_cell_per_eV')
        else:
            from pyscf.data.nist import HARTREE2EV
            from pyscf.scf import chkfile
            from wren.workflow.scf import run_single_point
            cfg = _scf_settings(settings)
            ids_first = None
            for index, data in enumerate(datasets):
                records = _molecules(data, settings.get('molecules'))
                ids = [str(r['id']) for r in records]
                if not records or (ids_first is not None and ids != ids_first):
                    raise ValueError('DOS frames require the same ordered nonempty molecule selection')
                ids_first = ids
                label = f"frame{data.get('frame', index):04d}"
                here = directory/label
                here.mkdir()
                frame_info = {'frame': data.get('frame'), 'time_ps': data.get('time_ps'), 'molecules': []}
                if settings.get('element-projections'):
                    element_frames.append([])
                if mode == 'vertical-dscf':
                    molecules = charging.calculate_frame(records, settings, here)
                    frame_info['molecules'] = molecules
                    spectra.append(np.array([[m['removal_ev'], m['addition_ev']] for m in molecules]))
                elif mode == 'monomer-orbitals':
                    levels, spectral_weights, homos, lumos = [], [], [], []
                    for r in records:
                        atoms = list(zip(r['elements'], r['coordinates_angstrom'], strict=True))
                        status = run_single_point(atoms, settings, here/str(r['id']))
                        mol, checkpoint = chkfile.load_scf(status['checkpoint'])
                        e = np.asarray(checkpoint['mo_energy'])*HARTREE2EV
                        occ = np.asarray(checkpoint['mo_occ'])
                        if e.ndim != 1 or not np.any(occ > 0) or not np.any(occ == 0):
                            raise ValueError('DOS currently requires restricted closed-shell orbitals with virtual states')
                        if settings.get('element-projections'):
                            element_frames[-1].append(element_probabilities(mol, checkpoint['mo_coeff'],
                                mol.intor_symmetric('int1e_ovlp'), backend=settings['dos-backend']))
                        homo, lumo = float(e[occ > 0].max()), float(e[occ == 0].min())
                        homos.append(homo); lumos.append(lumo)
                        levels.extend(e.tolist()); spectral_weights.extend([2/len(records)]*len(e))
                        frame_info['molecules'].append({'id': str(r['id']), 'homo_ev': homo, 'lumo_ev': lumo,
                            'orbital_energies_ev': e.tolist(), 'occupations': occ.tolist(), 'scf': status})
                    spectra.append(np.array(levels)); weights.append(np.array(spectral_weights))
                    homo_table.append(homos); lumo_table.append(lumos)
                    result.update(scope='isolated_molecular_orbital_ensemble_no_inter_molecular_coupling',
                                  normalization='spin_states_per_molecule_per_eV')
                else:
                    periodic_frame, levels, band_weights, groups = periodic.calculate(data, records, settings, cfg, here/'pbc.log')
                    frame_info.update(periodic_frame)
                    spectra.append(alignment.aligned_levels(settings, frame_info, levels, references))
                    weights.append(band_weights)
                    if groups:
                        element_frames[-1].extend(groups)
                    result.update(scope='periodic_restricted_HF_or_Kohn_Sham_spectrum_not_quasiparticle',
                                  normalization='spin_states_per_cell_per_eV')
                frame_records.append(frame_info)
                print(f'wren: DOS {label} spectrum completed', flush=True)
        minimum = min(float(e.min()) for e in spectra)-8*settings['sigma']
        maximum = max(float(e.max()) for e in spectra)+8*settings['sigma']
        lower, upper = settings.get('emin', minimum), settings.get('emax', maximum)
        count = settings.get('points', max(2, int(np.ceil((upper-lower)/(settings['sigma']/5)))+1))
        if lower >= upper or count > 1000000 or (upper-lower)/(count-1) > settings['sigma']/2:
            raise ValueError('DOS grid must resolve sigma with step <= sigma/2 and <= 1000000 points; adjust bounds or points')
        grid = np.linspace(lower, upper, count)
        if mode == 'vertical-dscf':
            result.update(charging.write_distributions(directory, frame_records, grid, settings))
        else:
            summary, total = _write_total_dos(directory, settings, grid, spectra, weights, frame_records)
            result.update(summary)
            if site_weights:
                result['site_projections'] = _write_site_dos(directory, settings, spectra, grid, site_weights,
                                                            total, overlaps is not None)
            if element_frames:
                result['element_projections'] = _write_element_dos(directory, settings, spectra, grid,
                                                                  element_frames, weights, total)
        _write_result(directory/'spectra.json', {'schema': 'wren-dos-spectra-1', 'frames': frame_records})
        if homo_table:
            result['frontier_disorder'] = {'homo': disorder_statistics(homo_table), 'lumo': disorder_statistics(lumo_table)}
        result.update(frame_count=len(frame_records), frames='spectra.json', grid_spacing_ev=float(grid[1]-grid[0]))
        return _complete(directory, result, started)
    except Exception as exc:
        result.update(status='failed', error={'type': type(exc).__name__, 'message': str(exc)})
        _write_result(output, result)
        raise
