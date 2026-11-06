"""Frozen-geometry fragment couplings from CPU PySCF or GPU4PySCF.

DIPRO uses converged dimer orbitals/Hamiltonian (Baumeier et al. 2010,
DOI:10.1039/C002337J). ``frozen_density`` uses a different Hamiltonian:
one Fock/KS construction at the sum of isolated neutral fragment densities.
The latter is an explicitly defined neutral frozen-density approximation,
NOT a charged FO-DFT method or a claim to reproduce every FO-DFT flavour.
See Schober, Reuter & Oberhofer, JCP 144, 054103 (2016),
DOI:10.1063/1.4940920, Sec. II (preprint Sec. II.3.1),
https://arxiv.org/html/1512.00200, and FHI-aims manual, Sec. 3.48.
Original abstract excerpt: "the construction of the Hamiltonian".
"""
from dataclasses import asdict, dataclass
import hashlib
import importlib.metadata
import json
import math
from pathlib import Path
import time

import numpy as np

from .coupling import ApplicabilityError, energy_splitting, gmh, host, lowdin

METHODS = ('dipro', 'frozen_density', 'esid', 'orbital_gmh')


@dataclass(frozen=True)
class SCFConfig:
    backend: str = 'cpu'
    basis: str = 'def2-svp'
    xc: str = 'b3lyp'  # 'hf' selects restricted Hartree-Fock
    auxbasis: str = 'def2-universal-jkfit'
    unit: str = 'Angstrom'
    conv_tol: float = 1e-10
    max_cycle: int = 100
    grid_level: int = 3
    threads: int = 4
    max_memory_mb: int = 4000
    degeneracy_tolerance_ev: float = 0.01
    site_tolerance_ev: float = 0.001
    max_pair_leakage: float = 0.15

    def validate(self):
        if self.backend not in ('cpu', 'gpu') or self.unit not in ('Angstrom', 'Bohr'):
            raise ValueError('backend must be cpu/gpu and unit Angstrom/Bohr')
        for name in ('basis', 'xc', 'auxbasis'):
            if not isinstance(getattr(self, name), str) or not getattr(self, name).strip():
                raise ValueError(f'{name} must be a nonempty string')
        for name in ('max_cycle', 'threads', 'max_memory_mb'):
            if type(getattr(self, name)) is not int or getattr(self, name) <= 0:
                raise ValueError(f'{name} must be a positive integer')
        if type(self.grid_level) is not int or not 0 <= self.grid_level <= 9:
            raise ValueError('grid_level must be an integer in [0,9]')
        for name in ('conv_tol', 'degeneracy_tolerance_ev', 'site_tolerance_ev'):
            if not math.isfinite(getattr(self, name)) or getattr(self, name) <= 0:
                raise ValueError(f'{name} must be positive and finite')
        if not math.isfinite(self.max_pair_leakage) or not 0 <= self.max_pair_leakage < 1:
            raise ValueError('max_pair_leakage must be in [0,1)')


def read_xyz(path):
    """Read exactly one XYZ frame, retaining its atom ordering."""
    lines = Path(path).read_text().splitlines()
    if len(lines) < 3:
        raise ValueError('XYZ needs an atom count, comment line, and atoms')
    n = int(lines[0])
    if n <= 0 or len(lines) < n+2 or any(s.strip() for s in lines[n+2:]):
        raise ValueError('XYZ must contain exactly one complete, nonempty frame')
    atoms = []
    for line in lines[2:n+2]:
        fields = line.split()
        if len(fields) != 4:
            raise ValueError('XYZ atom rows must contain element x y z')
        atoms.append((fields[0], tuple(map(float, fields[1:]))))
    return atoms


def _atoms(value):
    value = read_xyz(value) if isinstance(value, (str, Path)) else value
    atoms = [(str(symbol), tuple(float(x) for x in xyz)) for symbol, xyz in value]
    if not atoms or any(len(xyz) != 3 or not all(map(math.isfinite, xyz)) for _, xyz in atoms):
        raise ValueError('each fragment needs atoms with finite three-dimensional coordinates')
    return atoms


def _versions():
    versions = {}
    for name in ('numpy', 'pyscf', 'gpu4pyscf-cuda13x', 'cupy-cuda13x', 'cutensor-cu13'):
        try:
            versions[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            pass
    return versions


def calculate_transfer(fragment_a, fragment_b, *, config=None, methods=('dipro',),
                       carriers=('hole', 'electron'), assume_equivalent=False,
                       orbital_indices=None):
    """Compute couplings; explicit orbital_indices map carrier -> (A indices, B indices).

    Orbital indices are zero-based, in each isolated fragment's own MO order.
    Default near-degenerate HOMO/LUMO manifolds are selected together. Coupling
    blocks are returned without compressing them into an ambiguous scalar rate.
    Fragments are neutral closed-shell molecules; unsupported open shells fail.
    No coordinates are optimized, and no external MM embedding is implied.
    """
    started = time.perf_counter()
    config = config or SCFConfig()
    config.validate()
    methods, carriers = tuple(methods), tuple(carriers)
    if not methods or len(set(methods)) != len(methods) or any(m not in METHODS for m in methods):
        raise ValueError(f'methods must be a nonempty unique selection from {METHODS}')
    if not carriers or len(set(carriers)) != len(carriers) or any(c not in ('hole', 'electron') for c in carriers):
        raise ValueError('carriers must select hole and/or electron without duplicates')
    if type(assume_equivalent) is not bool:
        raise ValueError('assume_equivalent must be boolean')
    if orbital_indices is not None:
        if not isinstance(orbital_indices, dict) or any(k not in carriers for k in orbital_indices):
            raise ValueError('orbital_indices must map requested carriers to (A indices, B indices)')
        if any(not isinstance(v, (tuple, list)) or len(v) != 2 for v in orbital_indices.values()):
            raise ValueError('each explicit selection needs both A and B orbital indices')
    a, b = _atoms(fragment_a), _atoms(fragment_b)
    from pyscf import gto, lib
    from pyscf.data.nist import HARTREE2EV
    lib.num_threads(config.threads)
    xp = np
    device = None
    gpu_runtime = None
    if config.backend == 'gpu':
        import cupy as xp
        import gpu4pyscf  # registers PySCF to_gpu conversions; no fallback on failure
        device = xp.cuda.runtime.getDeviceProperties(xp.cuda.runtime.getDevice())['name']
        device = device.decode() if isinstance(device, bytes) else str(device)
        from gpu4pyscf.lib import cutensor
        gpu_runtime = {'driver_version': xp.cuda.runtime.driverGetVersion(),
                       'runtime_version': xp.cuda.runtime.runtimeGetVersion(),
                       'tensor_contraction_engine': cutensor.contract_engine or 'cutensor',
                       'projection_array_module': xp.__name__}

    def sync():
        if config.backend == 'gpu':
            xp.cuda.get_current_stream().synchronize()

    def molecule(atoms):
        # Same spherical/cartesian convention, basis and atom order in every SCF.
        mol = gto.M(atom=atoms, basis=config.basis, unit=config.unit,
                    charge=0, spin=0, cart=False, verbose=0,
                    max_memory=config.max_memory_mb)
        if any(mol.atom_charge(i) <= 0 for i in range(mol.natm)):
            raise ValueError('ghost atoms are not supported in this isolated-fragment basis convention')
        return mol

    ma, mb, md = molecule(a), molecule(b), molecule(a+b)
    na, nb = ma.nao_nr(), mb.nao_nr()
    if md.nao_nr() != na+nb or md.aoslice_by_atom()[len(a)-1, 3] != na:
        raise ValueError('fragment/dimer AO mapping is inconsistent')

    def mean_field(mol):
        mf = mol.RHF() if config.xc.lower() == 'hf' else mol.RKS(xc=config.xc)
        mf = mf.density_fit(auxbasis=config.auxbasis)
        mf.conv_tol, mf.max_cycle = config.conv_tol, config.max_cycle
        mf.direct_scf_tol = 1e-13
        if hasattr(mf, 'grids'):
            mf.grids.level = config.grid_level
        return mf.to_gpu() if config.backend == 'gpu' else mf

    timings, scf_records = {}, {}

    def run(mol, label):
        mf = mean_field(mol)
        sync()
        t = time.perf_counter()
        mf.kernel()
        sync()
        timings[label] = time.perf_counter()-t
        if not mf.converged or not math.isfinite(float(mf.e_tot)):
            raise RuntimeError(f'{label} SCF did not converge; no coupling is emitted')
        if config.backend == 'gpu' and not isinstance(mf.mo_coeff, xp.ndarray):
            raise RuntimeError('GPU request produced host MO coefficients')
        scf_records[label] = {'converged': bool(mf.converged), 'energy_hartree': float(mf.e_tot),
                              'cycles': int(mf.cycles) if hasattr(mf, 'cycles') else None,
                              'class': type(mf).__module__+'.'+type(mf).__name__}
        return mf

    mfa, mfb = run(ma, 'fragment_a'), run(mb, 'fragment_b')
    need_dimer = any(m != 'frozen_density' for m in methods)
    mfd = run(md, 'dimer') if need_dimer else mean_field(md)
    s = xp.asarray(mfd.get_ovlp())
    if not bool(xp.allclose(s[:na, :na], xp.asarray(mfa.get_ovlp()), atol=1e-11, rtol=1e-10)) or not bool(xp.allclose(s[na:, na:], xp.asarray(mfb.get_ovlp()), atol=1e-11, rtol=1e-10)):
        raise ValueError('isolated and dimer AO metrics do not match')
    operators = {}
    if need_dimer:
        operators['dipro'] = xp.asarray(mfd.get_fock(dm=mfd.make_rdm1()))
    if 'frozen_density' in methods:
        sync()
        t = time.perf_counter()
        dm = xp.zeros((na+nb, na+nb))
        dm[:na, :na], dm[na:, na:] = mfa.make_rdm1(), mfb.make_rdm1()
        # Baumeier et al. (2010), noSCF approximation, and the neutral-density
        # Hamiltonian distinctions discussed by Schober et al. (2016).
        # Here P=P_A (direct sum) P_B is explicitly the AO-density sum, not a
        # determinant formed by orthogonalizing all occupied fragment orbitals.
        # No call to kernel() on this density, and no occupation reset/charge claim.
        operators['frozen_density'] = xp.asarray(mfd.get_fock(dm=dm))
        sync()
        timings['frozen_density_fock'] = time.perf_counter()-t

    def select(mf, carrier, explicit=None):
        energies, occupied = host(mf.mo_energy), host(mf.mo_occ) > 0
        sector = np.flatnonzero(occupied if carrier == 'hole' else ~occupied)
        if not len(sector):
            raise ApplicabilityError(f'no {carrier} frontier orbitals available')
        frontier = sector[-1] if carrier == 'hole' else sector[0]
        if explicit is None:
            indices = sector[abs(energies[sector]-energies[frontier])*HARTREE2EV <= config.degeneracy_tolerance_ev].tolist()
        else:
            indices = list(explicit)
            if not indices or len(set(indices)) != len(indices) or any(type(i) is not int or i not in sector for i in indices):
                raise ValueError('explicit orbital indices must be unique zero-based indices in the carrier occupation sector')
        c = xp.asarray(mf.mo_coeff[:, indices]).copy()
        for j in range(c.shape[1]):
            anchor = int(xp.argmax(xp.abs(c[:, j])))
            if float(c[anchor, j]) < 0:
                c[:, j] *= -1  # deterministic display gauge; |V| is gauge invariant
        return indices, c

    output = {}
    max_projection_residual = 0.0
    for carrier in carriers:
        explicit = None if orbital_indices is None else orbital_indices.get(carrier)
        ia, ca = select(mfa, carrier, None if explicit is None else explicit[0])
        ib, cb = select(mfb, carrier, None if explicit is None else explicit[1])
        n_a, n_b = len(ia), len(ib)
        c = xp.zeros((na+nb, n_a+n_b))
        c[:na, :n_a], c[na:, n_a:] = ca, cb
        overlap = c.T @ s @ c
        projected, effective, orthogonalizer = {}, {}, None
        for name, operator in operators.items():
            h = c.T @ operator @ c
            eff, orthogonalizer, eigenvalues = lowdin(h, overlap)
            projected[name], effective[name] = h, eff
        results = {}
        for method in methods:
            if method in ('dipro', 'frozen_density'):
                block = effective[method][:n_a, n_a:]*HARTREE2EV
                values = host(block)
                results[method] = {
                    'status': 'ok', 'coupling_ev': values.tolist(),
                    'coupling_mev': (values*1000).tolist(),
                    'scalar_magnitude_ev': float(abs(values[0, 0])) if values.shape == (1, 1) else None,
                    'frobenius_norm_ev': float(xp.linalg.norm(block)),
                    'singular_values_ev': host(xp.linalg.svd(block, compute_uv=False)).tolist(),
                    'raw_hamiltonian_hartree': host(projected[method]).tolist(),
                    'orthogonal_hamiltonian_hartree': host(effective[method]).tolist(),
                    'density': 'converged_dimer' if method == 'dipro' else 'sum_isolated_neutral_fragments',
                }
                continue
            try:
                if n_a != 1 or n_b != 1:
                    raise ApplicabilityError('two-state method rejected: selected fragment manifolds contain more than one orbital per fragment')
                occ = host(mfd.mo_occ) > 0
                sector = np.flatnonzero(occ if carrier == 'hole' else ~occ)
                pair = sector[-2:] if carrier == 'hole' else sector[:2]
                if len(pair) != 2:
                    raise ApplicabilityError('two dimer frontier orbitals are required')
                pair_coeff = xp.asarray(mfd.mo_coeff[:, pair])
                weights = (c @ orthogonalizer).T @ s @ pair_coeff
                singular = xp.linalg.svd(weights, compute_uv=False)
                leakage = max(0.0, 1-float(singular[-1])**2)
                if leakage > config.max_pair_leakage:
                    raise ApplicabilityError('dimer frontier pair does not span the selected fragment subspace')
                energies = host(mfd.mo_energy[pair])
                if method == 'esid':
                    value = energy_splitting(energies, host(xp.diag(effective['dipro'])),
                                             equivalent=assume_equivalent,
                                             tolerance=config.site_tolerance_ev/HARTREE2EV)
                    details = {'equivalent_fragment_assumption': assume_equivalent}
                else:
                    # AO position integrals use CPU libcint; SCF and projection are GPU
                    # when requested. Nuclear dipole is a common diagonal offset and
                    # cancels in the GMH expression for orthonormal orbital states.
                    ao_dipole = -xp.asarray(md.intor('int1e_r', comp=3, hermi=1))
                    dipoles = xp.einsum('pi,xpq,qj->xij', pair_coeff, ao_dipole, pair_coeff)
                    details = gmh(energies, host(dipoles))
                    value = details.pop('magnitude')
                    details.update({'approximation': 'single_particle_orbital_GMH_not_many_electron_GMH',
                                    'dipoles_au': host(dipoles).tolist(), 'ao_dipole_integrals_backend': 'cpu_libcint'})
                results[method] = {'status': 'ok', 'scalar_magnitude_ev': value*HARTREE2EV,
                                   'scalar_magnitude_mev': value*HARTREE2EV*1000,
                                   'dimer_orbital_indices': pair.tolist(), 'pair_subspace_leakage': leakage, **details}
            except ApplicabilityError as exc:
                results[method] = {'status': 'not_applicable', 'reason': str(exc)}
        if need_dimer:
            # Complete, not truncated, dimer MO expansion checks AO mapping:
            # F_sub = U eps U^T, U=C_frag^T S C_dimer. Same operator as DIPRO.
            u = c.T @ s @ xp.asarray(mfd.mo_coeff)
            reconstructed = (u*xp.asarray(mfd.mo_energy)[None, :]) @ u.T
            residual = float(xp.max(xp.abs(reconstructed-projected['dipro'])))
            max_projection_residual = max(max_projection_residual, residual)
            if residual > 1e-6:
                raise RuntimeError('dimer Fock/MO projection residual exceeds 1e-6 Hartree')
        output[carrier] = {'orbital_indices_a': ia, 'orbital_indices_b': ib,
                           'fragment_overlap': host(overlap).tolist(), 'methods': results}
    sync()
    request = {'fragment_a': a, 'fragment_b': b, 'config': asdict(config),
               'methods': methods, 'carriers': carriers, 'assume_equivalent': assume_equivalent,
               'orbital_indices': orbital_indices}
    complete = all(r['status'] == 'ok' for c in output.values() for r in c['methods'].values())
    return {'schema': 'cube-transfer-1', 'status': 'completed' if complete else 'partial', 'input': request,
            'input_sha256': hashlib.sha256(json.dumps(request, sort_keys=True).encode()).hexdigest(),
            'backend': config.backend, 'device': device, 'gpu_runtime': gpu_runtime, 'versions': _versions(),
            'scf': scf_records, 'timings_seconds': {**timings, 'total': time.perf_counter()-started},
            'maximum_complete_mo_projection_residual_hartree': max_projection_residual if need_dimer else None,
            'results': output,
            'scope': 'neutral_closed_shell_frozen_geometry_no_MM_embedding'}
