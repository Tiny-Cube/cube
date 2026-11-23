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
from .core import broaden, bloch_spectra, disorder_statistics, kmesh


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


def run(settings, *, resume=False):
    from . import core
    started = time.perf_counter()
    mode = settings['mode']
    datasets = _dataset(settings) if mode != 'tight-binding' else []
    if mode == 'periodic' and len(datasets) > 1:
        raise ValueError('multi-frame periodic DOS needs explicit potential alignment, not yet implemented; run each frame separately')
    network = json.loads(Path(settings['input']).read_text()) if mode == 'tight-binding' else None
    sources = {key: file_hash(settings[key]) for key in ('input', 'xyz') if key in settings}
    request = {'settings': settings, 'datasets': datasets, 'network': network, 'source_sha256': sources,
               'implementation_sha256': {p: file_hash(p) for p in (__file__, core.__file__)},
               'numpy_version': np.__version__}
    if mode != 'tight-binding':
        import pyscf
        request['pyscf_version'] = pyscf.__version__
        from wren.workflow import scf as scf_module
        from wren.quantum import checkpoints
        request['implementation_sha256'].update({p: file_hash(p) for p in (scf_module.__file__, checkpoints.__file__)})
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
              'broadening': {'kernel': 'Gaussian', 'sigma_ev': settings['sigma']},
              'backend': settings['dos-backend'], 'frames': []}
    _write_result(output, result)
    spectra, weights, frame_records, homo_table, lumo_table = [], [], [], [], []
    try:
        if mode == 'tight-binding':
            if network.get('schema') != 'wren-dos-tb-1' or network.get('coupling_convention') != 'signed_real_consistent_orbital_phases':
                raise ValueError('tight-binding input needs wren-dos-tb-1 and explicit signed-coupling convention')
            points, kw = kmesh(settings['kmesh'])
            values = bloch_spectra(network['site_energies_ev'], network['pairs'], points,
                                   backend=settings['dos-backend'])
            spectra.append(values.ravel())
            weights.append(np.repeat(kw*settings['spin-degeneracy'], values.shape[1]))
            frame_records.append({'eigenvalues_ev': values.tolist(), 'kpoints_fractional': points.tolist(),
                                  'kpoint_weights': kw.tolist()})
            result.update(scope='orthonormal_effective_Hamiltonian', normalization='states_per_model_cell_per_eV')
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
                if mode == 'monomer-orbitals':
                    levels, spectral_weights, homos, lumos = [], [], [], []
                    for r in records:
                        atoms = list(zip(r['elements'], r['coordinates_angstrom']))
                        status = run_single_point(atoms, settings, here/str(r['id']))
                        _, checkpoint = chkfile.load_scf(status['checkpoint'])
                        e = np.asarray(checkpoint['mo_energy'])*HARTREE2EV
                        occ = np.asarray(checkpoint['mo_occ'])
                        if e.ndim != 1 or not np.any(occ > 0) or not np.any(occ == 0):
                            raise ValueError('DOS currently requires restricted closed-shell orbitals with virtual states')
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
                    if cfg.backend != 'cpu':
                        raise ValueError('periodic GPU DOS is not validated; request pyscf explicitly')
                    if settings.get('molecules'):
                        raise ValueError('periodic DOS must use the full cell, not a molecule subset')
                    if 'lattice_angstrom' not in data:
                        raise ValueError('periodic DOS needs an MD manifest with a full lattice')
                    atoms = []
                    for r in records:
                        center = np.asarray(r.get('transport_position_angstrom'), dtype=float)
                        if center.shape != (3,):
                            raise ValueError('periodic DOS needs original transport_position_angstrom; re-extract this frame')
                        atoms.extend(zip(r['elements'], (np.asarray(r['coordinates_angstrom'])+center).tolist()))
                    if len(atoms) > settings['max-atoms']:
                        raise ValueError('periodic cell exceeds max-atoms; choose an explicitly justified smaller cell or raise the limit')
                    from pyscf.pbc import gto, scf, dft
                    cell = gto.Cell(atom=atoms, a=data['lattice_angstrom'], unit='Angstrom',
                                    basis=cfg.basis, pseudo=settings['pseudo'], precision=settings['precision'],
                                    max_memory=cfg.max_memory_mb, verbose=4, output=str(here/'pbc.log'))
                    mf = None
                    try:
                        cell.build()
                        points, kw = kmesh(settings['kmesh'])
                        ks = cell.get_abs_kpts(points)
                        mf = (scf.KRHF(cell, ks) if cfg.xc.lower() == 'hf' else dft.KRKS(cell, ks, xc=cfg.xc)).density_fit(
                            auxbasis=None if cfg.auxbasis == 'auto' else cfg.auxbasis)
                        mf.conv_tol, mf.max_cycle = cfg.conv_tol, cfg.max_cycle
                        if hasattr(mf, 'grids'):
                            mf.grids.level = cfg.grid_level
                        from pyscf import lib
                        lib.num_threads(cfg.threads)
                        mf.kernel()
                        if not mf.converged:
                            raise RuntimeError('periodic SCF did not converge')
                        values = np.asarray(mf.mo_energy)*HARTREE2EV
                        spectra.append(values.ravel()); weights.append(np.repeat(2*kw, values.shape[1]))
                        frame_info.update(orbital_energies_ev=values.tolist(), occupations=np.asarray(mf.mo_occ).tolist(),
                                          kpoints_fractional=points.tolist(), kpoint_weights=kw.tolist(),
                                          total_energy_hartree=float(mf.e_tot), scf_converged=True)
                        result.update(scope='periodic_restricted_HF_or_Kohn_Sham_spectrum_not_quasiparticle',
                                      normalization='spin_states_per_cell_per_eV')
                    finally:
                        if cell.stdout is not None:
                            cell.stdout.close()
                frame_records.append(frame_info)
                print(f'wren: DOS {label} spectrum completed', flush=True)
        minimum = min(float(e.min()) for e in spectra)-8*settings['sigma']
        maximum = max(float(e.max()) for e in spectra)+8*settings['sigma']
        lower, upper = settings.get('emin', minimum), settings.get('emax', maximum)
        count = settings.get('points', max(2, int(np.ceil((upper-lower)/(settings['sigma']/5)))+1))
        if lower >= upper or count > 1000000 or (upper-lower)/(count-1) > settings['sigma']/2:
            raise ValueError('DOS grid must resolve sigma with step <= sigma/2 and <= 1000000 points; adjust bounds or points')
        grid = np.linspace(lower, upper, count)
        curves = [broaden(e, grid, settings['sigma'], weights=w, backend=settings['dos-backend'])
                  for e, w in zip(spectra, weights)]
        mean = np.mean(curves, axis=0)
        with (directory/'dos.csv').open('w', newline='') as file:
            writer = csv.writer(file); writer.writerow(['energy_ev', 'dos', *[f'frame_{i}' for i in range(len(curves))]])
            writer.writerows(zip(grid, mean, *curves))
        _write_result(directory/'spectra.json', {'schema': 'wren-dos-spectra-1', 'frames': frame_records})
        if homo_table:
            result['frontier_disorder'] = {'homo': disorder_statistics(homo_table), 'lumo': disorder_statistics(lumo_table)}
        expected = float(np.mean([w.sum() for w in weights]))
        result.update(status='completed', frame_count=len(curves), frames='spectra.json',
                      expected_integral=expected, finite_grid_integral=float(np.trapezoid(mean, grid)),
                      integral_relative_error=abs(float(np.trapezoid(mean, grid))-expected)/max(expected, 1e-300),
                      grid_spacing_ev=float(grid[1]-grid[0]),
                      sampling_uncertainty_estimated=False, elapsed_seconds=time.perf_counter()-started,
                      output_sha256={str(p.relative_to(directory)): file_hash(p) for p in directory.rglob('*') if p.is_file() and p != output})
        _write_result(output, result)
        return {**result, 'output': str(output), 'reused': False}
    except Exception as exc:
        result.update(status='failed', error={'type': type(exc).__name__, 'message': str(exc)})
        _write_result(output, result)
        raise
