"""Single-point jobs, frozen electrostatic embedding and vacuum translations."""
from dataclasses import asdict
import json
import math
from pathlib import Path
import time

import numpy as np

from universal_kmc.quantum.checkpoints import save_transfer_checkpoint
from universal_kmc.quantum.cli import _write_result
from universal_kmc.quantum.coupling import host
from universal_kmc.quantum.transfer import _versions
from .control import scf_config
from .dimers import file_hash, fingerprint


def make_mean_field(atoms, config, *, charge=0, spin=0, log=None, embedding=None):
    from pyscf import gto, lib
    config.validate()
    lib.num_threads(config.threads)
    mol = gto.M(atom=atoms, basis=config.basis, unit='Angstrom', charge=charge, spin=spin,
                cart=False, verbose=4 if log else 0, output=str(log) if log else None,
                max_memory=config.max_memory_mb)
    if config.xc.lower() == 'hf':
        mf = mol.UHF() if spin else mol.RHF()
    else:
        mf = mol.UKS(xc=config.xc) if spin else mol.RKS(xc=config.xc)
    mf = mf.density_fit(auxbasis=config.auxbasis)
    mf.conv_tol, mf.max_cycle, mf.direct_scf_tol = config.conv_tol, config.max_cycle, 1e-13
    if hasattr(mf, 'grids'):
        mf.grids.level = config.grid_level
    if config.backend == 'gpu':
        import gpu4pyscf
        mf = mf.to_gpu()
        if not type(mf).__module__.startswith('gpu4pyscf.'):
            raise RuntimeError('GPU request did not produce a GPU4PySCF mean-field object')
    if embedding is not None:
        embedding.validate_qm(atoms)
        if config.backend == 'gpu':
            from gpu4pyscf.qmmm.itrf import add_mm_charges
        else:
            from pyscf.qmmm import add_mm_charges
        mf = add_mm_charges(mf, embedding.coordinates_angstrom, embedding.charges_e, unit='Angstrom')
        embedding.assert_attached(mf)
    return mf


def scf_record(mf, backend):
    if not mf.converged or not math.isfinite(float(mf.e_tot)):
        raise RuntimeError('SCF did not converge; no successful checkpoint is emitted')
    device = None
    if backend == 'gpu':
        import cupy as cp
        if not isinstance(mf.mo_coeff, cp.ndarray):
            raise RuntimeError('GPU request produced host orbital coefficients')
        device = cp.cuda.runtime.getDeviceProperties(cp.cuda.runtime.getDevice())['name']
        device = device.decode() if isinstance(device, bytes) else str(device)
        cp.cuda.get_current_stream().synchronize()
    return {'converged': True, 'energy_hartree': float(mf.e_tot), 'backend': backend, 'device': device,
            'source_class': type(mf).__module__+'.'+type(mf).__name__,
            'charge': mf.mol.charge, 'spin_2s': mf.mol.spin,
            'spin_square': float(mf.spin_square()[0]) if mf.mol.spin else 0.0}


def run_single_point(atoms, settings, prefix, *, resume=False, charge=0, spin=0, embedding=None):
    from pyscf.scf import chkfile
    prefix = Path(prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    config = scf_config(settings)
    request = {'atoms': atoms, 'settings': asdict(config), 'charge': charge, 'spin_2s': spin,
               'versions': _versions(), 'implementation_sha256': file_hash(__file__)}
    if embedding is not None:
        from . import embedding as embedding_module
        request.update(environment='qmmm', embedding_fingerprint=embedding.fingerprint,
                       embedding_implementation_sha256=file_hash(embedding_module.__file__))
    key = fingerprint(request)
    checkpoint, status, log = (Path(str(prefix)+extension) for extension in ('.chk', '.status.json', '.log'))
    if status.exists():
        old = json.loads(status.read_text())
        if not resume or old.get('fingerprint') != key:
            raise ValueError(f'existing QC task {prefix}; use --resume with identical settings or a different output')
        if old.get('status') == 'completed':
            if not checkpoint.is_file() or file_hash(checkpoint) != old.get('checkpoint_sha256'):
                raise ValueError(f'cannot resume: checkpoint changed or missing: {checkpoint}')
            return {**old, 'reused': True}
    elif checkpoint.exists() or log.exists():
        raise ValueError(f'untracked output exists for {prefix}; choose another output')
    checkpoint.unlink(missing_ok=True)
    record = {'fingerprint': key, 'request': request, 'status': 'running',
              'checkpoint': str(checkpoint.resolve()), 'log': str(log.resolve())}
    _write_result(status, record)
    started = time.perf_counter()
    mf = None
    try:
        mf = make_mean_field(atoms, config, charge=charge, spin=spin, log=log, embedding=embedding)
        mf.kernel()
        details = scf_record(mf, config.backend)
        if embedding is not None:
            embedding.assert_attached(mf)
            details.update(environment='qmmm', mm_frozen=True, embedding_fingerprint=embedding.fingerprint)
        if charge == 0 and spin == 0 and embedding is None:
            save_transfer_checkpoint(mf, checkpoint)
        else:
            chkfile.dump_scf(mf.mol, str(checkpoint), float(mf.e_tot), host(mf.mo_energy), host(mf.mo_coeff), host(mf.mo_occ))
            marker = {'converged': True, 'charge': charge, 'spin_2s': spin,
                      'environment': 'vacuum' if embedding is None else 'qmmm'}
            if embedding is not None:
                marker.update(embedding_fingerprint=embedding.fingerprint, mm_frozen=True,
                              mm_coordinates_angstrom=np.asarray(embedding.coordinates_angstrom),
                              mm_charges_e=np.asarray(embedding.charges_e))
            chkfile.save(str(checkpoint), 'cube_reorganization', marker)
        record.update(details, status='completed', checkpoint_sha256=file_hash(checkpoint),
                      elapsed_seconds=time.perf_counter()-started, reused=False)
        _write_result(status, record)
        return record
    except Exception as exc:
        checkpoint.unlink(missing_ok=True)
        record.update(status='failed', error={'type': type(exc).__name__, 'message': str(exc)})
        _write_result(status, record)
        raise
    finally:
        if mf is not None and mf.mol.stdout is not None:
            mf.mol.stdout.close()


def translated_checkpoint(source, target_atoms, destination):
    """Rebase an isolated neutral molecule by a PURE TRANSLATION, never rotation.

    Vacuum atom-centered AO coefficients/Fock are invariant to rigid translation.
    Source files remain unchanged. Every source atom must match its target.
    This function is used only for cube's vacuum jobs, not embedded Hamiltonians.
    """
    from pyscf.scf import chkfile
    source, destination = Path(source), Path(destination)
    if source.resolve() == destination.resolve():
        raise ValueError('rebased checkpoint must not overwrite its source')
    mol, data = chkfile.load_scf(str(source))
    metadata = chkfile.load(str(source), 'cube_transfer')
    if mol.charge != 0 or mol.spin != 0 or not metadata or not metadata.get('converged'):
        raise ValueError('rebase requires a successful neutral closed-shell cube vacuum checkpoint')
    symbols = [a[0] for a in target_atoms]
    if symbols != [mol.atom_symbol(i) for i in range(mol.natm)]:
        raise ValueError('rebased checkpoint element/order mismatch')
    coordinates = np.array([a[1] for a in target_atoms], dtype=float)
    shift = coordinates-mol.atom_coords(unit='Angstrom')
    if not np.allclose(shift, shift[0], atol=1e-7, rtol=0):
        raise ValueError('monomer and dimer fragment differ by more than a pure translation')
    mol.verbose = 0
    mol.set_geom_(coordinates, unit='Angstrom')
    destination.parent.mkdir(parents=True, exist_ok=True)
    chkfile.dump_scf(mol, str(destination), data['e_tot'], data['mo_energy'], data['mo_coeff'], data['mo_occ'])
    metadata.update({'vacuum_translation_angstrom': shift[0], 'translation_source_sha256': file_hash(source),
                     'translation_source_path': str(source.resolve())})
    chkfile.save(str(destination), 'cube_transfer', metadata)
    return destination
