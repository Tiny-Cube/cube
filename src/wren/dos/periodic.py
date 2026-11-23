"""Restricted periodic HF/KS orbitals on a full MD or crystal cell.

Use Gaussian density fitting and the same atom-centred XC grid on both devices.
See Sun et al., J. Chem. Phys. 147, 164119 (2017), doi:10.1063/1.4998644.
"""
import time
import numpy as np

from wren.quantum.coupling import host
from .core import element_probabilities, kmesh


def _full_cell(data, records, max_atoms):
    lattice = np.asarray(data.get('lattice_angstrom'), dtype=float)
    if lattice.shape != (3, 3) or not np.isfinite(lattice).all() or np.linalg.det(lattice) <= 0:
        raise ValueError('periodic DOS needs a finite right-handed full lattice in Angstrom')
    atoms = []
    for r in records:
        center = np.asarray(r.get('transport_position_angstrom'), dtype=float)
        if center.shape != (3,) or not np.isfinite(center).all():
            raise ValueError('periodic DOS needs original transport_position_angstrom; re-extract this frame')
        local = np.asarray(r['coordinates_angstrom'], dtype=float)
        if local.shape != (len(r['elements']), 3) or not len(local) or not np.isfinite(local).all():
            raise ValueError('periodic molecular coordinates must match elements and be finite')
        atoms.extend(zip(r['elements'], (local+center).tolist(), strict=True))
    if len(atoms) > max_atoms:
        raise ValueError('periodic cell exceeds max-atoms; choose an explicitly justified smaller cell or raise the limit')
    return atoms, lattice


def band_edges(energies, occupations):
    """Frontier energies from sampled k points, with zero gap for partial filling."""
    occ_values, virtual_values, direct = [], [], []
    partial = False
    for i, (e, occ) in enumerate(zip(energies, occupations, strict=True)):
        if e.ndim != 1 or e.shape != occ.shape or not np.isfinite(e).all() or not np.isfinite(occ).all():
            raise ValueError('invalid periodic energies or occupations')
        if np.any(occ < -1e-8) or np.any(occ > 2+1e-8):
            raise ValueError('invalid restricted periodic orbital occupation')
        partial |= bool(np.any((occ > 1e-8) & (occ < 2-1e-8)))
        occupied, available = e[occ > 1e-8], e[occ < 2-1e-8]
        if occupied.size:
            occ_values.append((float(occupied.max()), i))
        if available.size:
            virtual_values.append((float(available.min()), i))
        if occupied.size and available.size:
            direct.append(max(0., float(available.min()-occupied.max())))
    vbm, vk = max(occ_values, key=lambda p: p[0]) if occ_values else (None, None)
    cbm, ck = min(virtual_values, key=lambda p: p[0]) if virtual_values else (None, None)
    gap = max(0., cbm-vbm) if cbm is not None and vbm is not None else None
    return {'vbm_ev': vbm, 'vbm_k_index': vk, 'cbm_ev': cbm, 'cbm_k_index': ck,
            'sampled_gap_ev': gap, 'minimum_sampled_direct_gap_ev': min(direct) if direct else None,
            'partial_occupation': partial}


def calculate(data, records, settings, cfg, log_path):
    from pyscf import lib
    from pyscf.data.nist import HARTREE2EV
    from pyscf.pbc import gto
    if settings.get('molecules'):
        raise ValueError('periodic DOS must use the full cell')
    atoms, lattice = _full_cell(data, records, settings['max-atoms'])
    gpu = cfg.backend == 'gpu'
    if gpu:
        import cupy
        from gpu4pyscf.pbc.scf.khf import KRHF
        from gpu4pyscf.pbc.dft.krks import KRKS
        from gpu4pyscf.pbc.dft.gen_grid import BeckeGrids
    else:
        from pyscf.pbc.scf.khf import KRHF
        from pyscf.pbc.dft.krks import KRKS
        from pyscf.pbc.dft.gen_grid import BeckeGrids
    lib.num_threads(cfg.threads)
    points, kw = kmesh(settings['kmesh'])
    with log_path.open('w') as log:
        cell = gto.Cell(atom=atoms, a=lattice, unit='Angstrom', basis=cfg.basis,
                        pseudo=settings['pseudo'], precision=settings['precision'],
                        max_memory=cfg.max_memory_mb, verbose=4, stdout=log).build()
        ks = cell.get_abs_kpts(points)
        hf = cfg.xc.lower() == 'hf'
        mf = KRHF(cell, ks) if hf else KRKS(cell, ks, xc=cfg.xc)
        mf = mf.density_fit(auxbasis=None if cfg.auxbasis == 'auto' else cfg.auxbasis)
        if not hf:
            mf.grids = BeckeGrids(cell)
            mf.grids.level = cfg.grid_level
        mf.conv_tol, mf.max_cycle = cfg.conv_tol, cfg.max_cycle
        started = time.perf_counter()
        mf.kernel()
        if not mf.converged or not np.isfinite(float(mf.e_tot)):
            raise RuntimeError('periodic SCF did not converge to a finite energy')
        device = None
        if gpu:
            if (not type(mf).__module__.startswith('gpu4pyscf.') or
                    not all(isinstance(c, cupy.ndarray) for c in mf.mo_coeff)):
                raise RuntimeError('periodic GPU SCF did not produce device orbitals')
            cupy.cuda.get_current_stream().synchronize()
            name = cupy.cuda.runtime.getDeviceProperties(cupy.cuda.runtime.getDevice())['name']
            device = name.decode() if isinstance(name, bytes) else str(name)
        seconds = time.perf_counter()-started
        energies = [host(e)*HARTREE2EV for e in mf.mo_energy]
        occupations = [host(occ) for occ in mf.mo_occ]
        if len(energies) != len(points) or not all(len(e) for e in energies):
            raise ValueError('missing periodic orbital spectrum at a sampled k point')
        edges = band_edges(energies, occupations)
        electron_count = float(kw @ [occ.sum() for occ in occupations])
        if not np.isclose(electron_count, cell.nelectron, rtol=0, atol=1e-6):
            raise ValueError(f'periodic orbital occupations give {electron_count:g} electrons; expected {cell.nelectron}')
        groups = []
        if settings.get('element-projections'):
            for c, s in zip(mf.mo_coeff, mf.get_ovlp(), strict=True):
                if settings['dos-backend'] == 'cpu' and gpu:
                    c, s = host(c), host(s)
                groups.append(element_probabilities(cell, c, s, backend=settings['dos-backend']))
        frame = {'orbital_energies_ev': [e.tolist() for e in energies],
                 'occupations': [occ.tolist() for occ in occupations], 'band_edges': edges,
                 'kpoints_fractional': points.tolist(), 'kpoint_weights': kw.tolist(),
                 'total_energy_hartree': float(mf.e_tot), 'scf_converged': True,
                 'scf': {'backend': cfg.backend, 'device': device, 'class': type(mf).__module__+'.'+type(mf).__name__,
                         'density_fit_class': type(mf.with_df).__module__+'.'+type(mf.with_df).__name__,
                         'xc_grid': None if hf else 'BeckeGrids', 'grid_level': None if hf else cfg.grid_level,
                         'elapsed_seconds': seconds, 'modeled_electrons_per_cell': cell.nelectron,
                         'occupation_integral_per_cell': electron_count, 'atoms_per_cell': cell.natm}}
    weights = np.concatenate([np.full(len(e), 2*w) for e, w in zip(energies, kw, strict=True)])
    return frame, np.concatenate(energies), weights, groups
