"""Couplings from existing PySCF checkpoints: no SCF kernel is called here.

Standard PySCF SCF checkpoints contain mol, mo_coeff, mo_energy, mo_occ,
e_tot (https://pyscf.org/_modules/pyscf/scf/chkfile.html). For a complete
S-orthonormal MO basis, F = S C diag(epsilon) C.T S. This reproduces the
saved eigenoperator without guessing a functional or rebuilding SCF.
"""
import hashlib
import math
from pathlib import Path
import time
from types import SimpleNamespace

import numpy as np

from .control import CHECKPOINT_METHODS, validate_output
from .coupling import host
from .transfer import _evaluate_transfer


def save_transfer_checkpoint(mf, filename):
    """Save a converged RHF/RKS result, including convergence and Fock data.

    Call after mf.kernel(); GPU4PySCF arrays are transferred explicitly.
    The checkpoint remains readable by pyscf.scf.chkfile.load_scf.
    """
    from pyscf.scf import chkfile
    if not bool(mf.converged) or not math.isfinite(float(mf.e_tot)):
        raise ValueError('refusing to export an unconverged SCF result')
    if mf.mol.charge != 0 or mf.mol.spin != 0 or host(mf.mo_coeff).ndim != 2:
        raise ValueError('checkpoint export supports neutral restricted closed-shell RHF/RKS only')
    if getattr(mf, 'mm_mol', None) is not None or getattr(mf, 'with_solvent', None) is not None:
        raise ValueError('this output template does not support QM/MM or implicit-solvent calculations')
    path = Path(filename).expanduser().resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    # Evaluate on the original backend before transferring small metadata/arrays.
    fock = host(mf.get_fock(dm=mf.make_rdm1()))
    chkfile.dump_scf(mf.mol, str(path), float(mf.e_tot), host(mf.mo_energy),
                     host(mf.mo_coeff), host(mf.mo_occ))
    chkfile.save(str(path), 'cube_transfer', {
        'schema_version': 1, 'converged': True,
        'model': str(getattr(mf, 'xc', 'hf')).lower(),
        'source_class': type(mf).__module__+'.'+type(mf).__name__,
        'fock_hartree': fock,
    })
    return path


def _real_array(value, label):
    array = np.asarray(value)
    if np.iscomplexobj(array) or not np.issubdtype(array.dtype, np.number):
        raise ValueError(f'{label} must contain real numerical data')
    array = np.asarray(array, dtype=float)
    if not np.all(np.isfinite(array)):
        raise ValueError(f'{label} contains nonfinite data')
    return array


def _text(value):
    return value.decode() if isinstance(value, bytes) else str(value)


def _load_checkpoint(path, assume_converged):
    from pyscf.scf import chkfile
    try:
        mol, data = chkfile.load_scf(str(path))
        metadata = chkfile.load(str(path), 'cube_transfer')
    except (OSError, KeyError, TypeError, ValueError) as exc:
        raise ValueError(f'cannot read PySCF SCF checkpoint {path}: {exc}') from exc
    if not isinstance(data, dict) or not {'e_tot', 'mo_coeff', 'mo_energy', 'mo_occ'} <= data.keys():
        raise ValueError(f'{path}: missing complete SCF orbital data')
    if mol.charge != 0 or mol.spin != 0 or mol.cart or mol._ecp or getattr(mol, '_pseudo', {}):
        raise ValueError('saved-result mode requires neutral closed-shell spherical all-electron molecular calculations')
    if any(mol.atom_charge(i) <= 0 for i in range(mol.natm)):
        raise ValueError('ghost atoms are unsupported in this fragment basis convention')
    n = mol.nao_nr()
    c, e, occ = (_real_array(data[key], key) for key in ('mo_coeff', 'mo_energy', 'mo_occ'))
    if c.shape != (n, n) or e.shape != (n,) or occ.shape != (n,):
        raise ValueError('complete square restricted MO coefficients and full energy/occupation arrays are required; truncated checkpoints are unsupported')
    if not np.all(np.isin(occ, [0., 2.])) or occ.sum() != mol.nelectron:
        raise ValueError('restricted integer closed-shell occupations are required')
    if np.any(np.diff(e) < -1e-10) or np.any(np.diff(occ) > 0):
        raise ValueError('MOs must be energy ordered with occupied orbitals before virtual orbitals')
    energy = _real_array(data['e_tot'], 'e_tot')
    if energy.shape != ():
        raise ValueError('e_tot must be a scalar')
    if metadata is not None:
        if not isinstance(metadata, dict) or metadata.get('schema_version') != 1:
            raise ValueError('unsupported cube_transfer checkpoint metadata')
        flag = metadata.get('converged')
        if not isinstance(flag, (bool, np.bool_)) or not flag:
            raise ValueError('checkpoint explicitly reports unconverged SCF; assumption cannot override it')
        convergence = 'verified_export_marker'
    elif not assume_converged:
        raise ValueError('standard PySCF .chk has no convergence flag; use save_transfer_checkpoint or explicitly set calculation.assume_converged=true after verifying the SCF log')
    else:
        convergence = 'user_assumed_from_external_log'
    s = mol.intor_symmetric('int1e_ovlp')
    metric_error = float(np.max(np.abs(c.T@s@c-np.eye(n))))
    if metric_error > 1e-7:
        raise ValueError(f'MO coefficients fail AO-metric orthonormality: {metric_error}')
    f = (s@c*e[None, :])@c.T@s
    fock_residual = None
    operator_source = 'complete_saved_dimer_MO_eigenoperator'
    if metadata is not None and 'fock_hartree' in metadata:
        saved_fock = _real_array(metadata['fock_hartree'], 'fock_hartree')
        if saved_fock.shape != (n, n) or not np.allclose(saved_fock, saved_fock.T, atol=1e-11, rtol=1e-10):
            raise ValueError('saved Fock matrix has invalid shape or symmetry')
        fock_residual = float(np.max(np.abs(saved_fock-f)))
        if fock_residual > 1e-6:
            raise ValueError('saved Fock and complete MO eigenoperator disagree by more than 1e-6 Hartree')
        f = saved_fock
        operator_source = 'saved_fock_matrix'
    record = {'path': str(path), 'sha256': _sha256(path), 'energy_hartree': float(energy),
              'convergence': convergence, 'converged': True, 'ao_count': n,
              'mo_metric_error': metric_error, 'saved_fock_residual_hartree': fock_residual,
              'operator_source': operator_source,
              'model': _text(metadata['model']) if metadata and 'model' in metadata else None,
              'class': _text(metadata['source_class']) if metadata and 'source_class' in metadata else 'unknown_from_standard_checkpoint'}
    return SimpleNamespace(mol=mol, mo_coeff=c, mo_energy=e, mo_occ=occ), s, f, record


def _sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as file:
        for chunk in iter(lambda: file.read(1024*1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def _ao_mapping(ma, mb, md, selections):
    if selections is None:
        selections = (list(range(ma.natm)), list(range(ma.natm, ma.natm+mb.natm)))
    ia, ib = selections
    flat = ia+ib
    if len(ia) != ma.natm or len(ib) != mb.natm or len(flat) != md.natm or set(flat) != set(range(md.natm)):
        raise ValueError('fragment atom mappings must cover the dimer exactly, in each monomer checkpoint atom order')
    order = []
    for monomer, mapping in ((ma, ia), (mb, ib)):
        for source, target in enumerate(mapping):
            if monomer.atom_charge(source) != md.atom_charge(target) or not np.allclose(monomer.atom_coord(source), md.atom_coord(target), atol=1e-7, rtol=0):
                raise ValueError('monomer atoms must have the same elements and frozen coordinates as their mapped dimer atoms (common coordinate frame)')
            shells_a = [i for i in range(monomer.nbas) if monomer.bas_atom(i) == source]
            shells_d = [i for i in range(md.nbas) if md.bas_atom(i) == target]
            if len(shells_a) != len(shells_d):
                raise ValueError('fragment/dimer per-atom basis sets differ')
            for sa, sd in zip(shells_a, shells_d):
                if (monomer.bas_angular(sa), monomer.bas_nprim(sa), monomer.bas_nctr(sa)) != (md.bas_angular(sd), md.bas_nprim(sd), md.bas_nctr(sd)):
                    raise ValueError('fragment/dimer basis shell layouts differ')
                if not np.allclose(monomer.bas_exp(sa), md.bas_exp(sd), atol=1e-12, rtol=1e-12) or not np.allclose(monomer.bas_ctr_coeff(sa), md.bas_ctr_coeff(sd), atol=1e-12, rtol=1e-12):
                    raise ValueError('fragment/dimer basis exponents or contractions differ')
            start, stop = md.aoslice_by_atom()[target, 2:4]
            order.extend(range(int(start), int(stop)))
    if len(order) != ma.nao_nr()+mb.nao_nr() or len(order) != md.nao_nr():
        raise ValueError('AO mapping is incomplete')
    if md.nelectron != ma.nelectron+mb.nelectron:
        raise ValueError('dimer electron count differs from the sum of the fragments')
    return order


def calculate_from_checkpoints(job):
    """Evaluate a parsed TransferControl without constructing or running SCF."""
    from pyscf import lib
    started = time.perf_counter()
    job.config.validate()
    if any(method not in CHECKPOINT_METHODS for method in job.methods):
        raise ValueError('saved checkpoints do not define a frozen-density Hamiltonian')
    validate_output(job.output, [Path(job.provenance['path']), *job.checkpoints.values()])
    lib.num_threads(job.config.threads)
    loaded = {name: _load_checkpoint(job.checkpoints[name], job.assume_converged)
              for name in ('fragment_a', 'fragment_b', 'dimer')}
    mfa, sa, _, _ = loaded['fragment_a']
    mfb, sb, _, _ = loaded['fragment_b']
    mfd, s, fock, _ = loaded['dimer']
    records = {name: value[3] for name, value in loaded.items()}
    models = [record['model'] for record in records.values()]
    known_models = {model for model in models if model is not None}
    if len(known_models) > 1:
        raise ValueError('fragment and dimer checkpoint electronic methods differ')
    order = _ao_mapping(mfa.mol, mfb.mol, mfd.mol, job.atom_indices)
    s, fock = s[np.ix_(order, order)], fock[np.ix_(order, order)]
    mfd.mo_coeff = mfd.mo_coeff[order, :]
    na = mfa.mol.nao_nr()
    if not np.allclose(s[:na, :na], sa, atol=1e-11, rtol=1e-10) or not np.allclose(s[na:, na:], sb, atol=1e-11, rtol=1e-10):
        raise ValueError('fragment and dimer AO overlap metrics do not match')
    xp, device, gpu_runtime = np, None, None
    if job.config.backend == 'gpu':
        import cupy as xp
        device = xp.cuda.runtime.getDeviceProperties(xp.cuda.runtime.getDevice())['name']
        device = _text(device)
        gpu_runtime = {'driver_version': xp.cuda.runtime.driverGetVersion(),
                       'runtime_version': xp.cuda.runtime.runtimeGetVersion(),
                       'projection_array_module': xp.__name__}
    s, fock = xp.asarray(s), xp.asarray(fock)
    a, b = ([(mol.atom_symbol(i), tuple(mol.atom_coord(i, unit='Angstrom')))
             for i in range(mol.natm)] for mol in (mfa.mol, mfb.mol))
    result = _evaluate_transfer(mfa, mfb, mfd, {'dipro': fock}, s, config=job.config,
                                methods=job.methods, carriers=job.carriers,
                                assume_equivalent=job.assume_equivalent, orbital_indices=job.orbital_indices,
                                xp=xp, scf_records=records, timings={}, started=started, device=device,
                                gpu_runtime=gpu_runtime, a=a, b=b, ao_order=order)
    # SCF settings come from the checkpoints, never from SCFConfig defaults.
    result['input'].pop('config')
    result['input']['postprocessing'] = {
        'backend': job.config.backend, 'threads': job.config.threads,
        'degeneracy_tolerance_ev': job.config.degeneracy_tolerance_ev,
        'site_tolerance_ev': job.config.site_tolerance_ev,
        'max_pair_leakage': job.config.max_pair_leakage,
        'atom_indices_zero_based': job.atom_indices,
    }
    result['input']['checkpoints'] = {name: {'path': record['path'], 'sha256': record['sha256']}
                                      for name, record in records.items()}
    import json
    result['input_sha256'] = hashlib.sha256(json.dumps(result['input'], sort_keys=True).encode()).hexdigest()
    result.update({'workflow': 'saved_pyscf_checkpoints', 'scf_executed': False,
                   'scope': 'neutral_closed_shell_frozen_geometry_saved_operator',
                   'environment_model': 'not_reconstructed_or_verified_from_saved_results',
                   'operator_source': records['dimer']['operator_source'],
                   'control': job.provenance,
                   'model_consistency': 'verified_export_metadata' if all(models) else 'not_available_in_standard_checkpoint',
                   'checkpoint_loading_and_ao_integrals_backend': 'cpu',
                   'projection_residual_note': ('independent saved Fock versus complete MO projection'
                       if records['dimer']['operator_source'] == 'saved_fock_matrix' else
                       'consistency identity for a reconstructed eigenoperator; not an independent SCF convergence test')})
    return result
