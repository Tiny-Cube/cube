"""On-the-fly adiabatic energies, state gradients and real SA-CASSCF NACs."""
from dataclasses import dataclass
from contextlib import nullcontext

import numpy as np

from .engine import AdiabaticFrame


@dataclass
class Wavefunctions:
    mol: object
    mo: np.ndarray
    ci: list
    ncore: int
    ncas: int
    nelecas: tuple


def wavefunction_overlaps(left, right):
    """Full CI-state overlaps, including both spin sectors of the inactive core."""
    from pyscf import fci, gto
    n = left.ncore
    occupied = n + left.ncas
    cross_ao = gto.intor_cross('int1e_ovlp', left.mol, right.mol)
    cross = left.mo[:, :occupied].T @ cross_ao @ right.mo[:, :occupied]
    det_core = 1.
    active = cross[n:, n:]
    if n:
        core = cross[:n, :n]
        if np.linalg.svd(core, compute_uv=False)[-1] < 1e-8:
            raise RuntimeError('inactive-core overlap nearly singular; inspect active space or reduce dt')
        # Each occupied determinant contains the same doubly occupied core.
        # Apply the Schur determinant identity separately for alpha and beta.
        det_core = np.linalg.det(core)
        active = active - cross[n:, :n] @ np.linalg.solve(core, cross[:n, n:])
    return det_core**2 * np.array([[fci.addons.overlap(a, b, left.ncas, left.nelecas, s=active)
                                     for b in right.ci] for a in left.ci])


def align_roots(previous, current, minimum):
    if previous is None:
        # Deterministic initial gauge; subsequent steps use many-electron overlaps.
        signs = [np.sign(c.flat[np.argmax(abs(c))]) for c in current.ci]
        overlap = None
    else:
        overlap = wavefunction_overlaps(previous, current)
        diagonal = np.diag(overlap)
        competitors = np.abs(overlap).copy()
        np.fill_diagonal(competitors, -1)
        if (not np.isfinite(overlap).all() or np.min(abs(diagonal)) < minimum or
                np.any(abs(diagonal) <= np.maximum(competitors.max(axis=0), competitors.max(axis=1)) + 1e-8)):
            raise RuntimeError('SA-CASSCF energy-ordered roots lose overlap/identity; reduce dt or review active space')
        signs = np.sign(diagonal)
    # Transform the CI coefficients themselves before differentiating the states.
    current.ci = [s*c for s, c in zip(signs, current.ci)]
    return None if overlap is None else (overlap * signs).tolist()


class SACASSCF:
    def __init__(self, config):
        self.config = config

    def evaluate_adiabatic(self, q, previous=None):
        from pyscf import fci, gto, lib, mcscf, scf
        from pyscf.nac.sacasscf import NonAdiabaticCouplings

        s = self.config
        coords = np.asarray(q, dtype=float).reshape(len(s['symbols']), 3)
        if not np.isfinite(coords).all():
            raise ValueError('nonfinite atomic coordinates')
        lib.num_threads(s['threads'])
        spin = s['multiplicity'] - 1
        mol = gto.M(atom=list(zip(s['symbols'], coords)), basis=s['basis'], charge=s['charge'],
                    spin=spin, unit='Bohr', verbose=0, max_memory=s['max_memory_mb'])
        if mol.nelectron < s['active_electrons'] or (mol.nelectron - s['active_electrons']) % 2:
            raise ValueError('electron counts must leave a doubly occupied inactive core')
        mf = scf.ROHF(mol) if spin else scf.RHF(mol)
        mf.conv_tol = 1e-11
        mf.max_cycle = s['max_cycles']
        mf.kernel()
        if not mf.converged:
            raise RuntimeError('reference SCF did not converge at this geometry')
        ne = s['active_electrons']
        mc = mcscf.CASSCF(mf, s['active_orbitals'], ((ne+spin)//2, (ne-spin)//2))
        if mc.ncore + mc.ncas > mf.mo_coeff.shape[1]:
            raise ValueError('basis has too few orbitals for the supplied core and active space')
        mc.fcisolver = fci.direct_spin1.FCI(mol) if spin else fci.direct_spin0.FCI(mol)
        mc.fcisolver.spin = spin
        mc.state_average_([1/s['states']] * s['states'])
        mc.fcisolver.conv_tol = 1e-12
        mc.conv_tol, mc.conv_tol_grad, mc.max_cycle_macro = 1e-11, s['orbital_gradient_tolerance'], s['max_cycles']
        prior = previous.wavefunctions if previous is not None else None
        if prior is not None:
            initial_mo = mcscf.addons.project_init_guess(mc, prior.mo, prev_mol=prior.mol,
                                                         priority='active', use_hf_core=False)
        elif s.get('initial_active_mo_1based'):
            if max(s['initial_active_mo_1based']) > mf.mo_coeff.shape[1]:
                raise ValueError('initial active MO index exceeds orbital count')
            initial_mo = mc.sort_mo(s['initial_active_mo_1based'], base=1)
        else:
            initial_mo = mf.mo_coeff
        mc.kernel(initial_mo)
        energies = np.asarray(mc.e_states)
        if (not mc.converged or not np.all(mc.fcisolver.converged) or energies.shape != (s['states'],) or
                not np.isfinite(energies).all()):
            raise RuntimeError('SA-CASSCF or its individual CI states did not converge')
        if np.min(np.diff(energies)) <= s['min_gap_hartree']:
            raise RuntimeError('SA-CASSCF unordered or near-degenerate states; singular single-surface propagation')
        spin_squared = [float(fci.spin_op.spin_square0(ci, mc.ncas, mc.nelecas)[0]) for ci in mc.ci]
        expected_spin_squared = .5*spin*(.5*spin+1)
        if not np.allclose(spin_squared, expected_spin_squared, atol=1e-6, rtol=0):
            raise RuntimeError('tracked roots left the requested total-spin sector')
        wave = Wavefunctions(mol, mc.mo_coeff, mc.ci, mc.ncore, mc.ncas, mc.nelecas)
        overlaps = align_roots(prior, wave, s['min_state_overlap'])
        mc.ci = wave.ci
        n = len(energies)
        gradients = np.empty((n, coords.size))
        nac = np.zeros((n, n, coords.size))
        # The restricted CAS NAC driver subtracts the derivative of a *closed*,
        # doubly occupied core. PySCF 2.15 dispatches that contraction through
        # the reference's gradient class, whose ROHF implementation expects
        # spin-separated densities. Supply the spatial RHF integral-gradient
        # driver just during the NAC evaluation. No RHF ground-state solve is
        # performed; the converged ROHF reference, doublet CI and CAS response
        # are retained. Signed overlap finite differences test this convention.
        core_grad = scf.hf.RHF(mol).nuc_grad_method if spin else None
        for i in range(n):
            grad = mc.nuc_grad_method()
            gradients[i] = grad.kernel(state=i).ravel()
            if not grad.converged:
                raise RuntimeError(f'SA-CASSCF state {i} gradient response did not converge')
            for j in range(i):
                response = NonAdiabaticCouplings(mc)
                # Request the energy-weighted interstate derivative. Impose our
                # own denominator for d_ji=<j|∂i>, verified against signed
                # central differences of full many-electron overlaps.
                with lib.temporary_env(mf, nuc_grad_method=core_grad) if spin else nullcontext():
                    weighted = response.kernel(state=(i,j), use_etfs=s['use_etfs'], mult_ediff=True)
                nac[j, i] = weighted.ravel() / (energies[i] - energies[j])
                if not response.converged:
                    raise RuntimeError(f'SA-CASSCF derivative coupling response failed for {j},{i}')
                nac[i, j] = -nac[j, i]
        if not np.isfinite(gradients).all() or not np.isfinite(nac).all():
            raise RuntimeError('SA-CASSCF returned a non-finite force or derivative coupling')
        diagnostics = {'adjacent_wavefunction_overlaps': overlaps, 'spin_squared': spin_squared,
                       'response_converged': True, 'minimum_gap_hartree': float(np.min(np.diff(energies)))}
        occupations = [np.linalg.eigvalsh(fci.direct_spin1.make_rdm1(c, mc.ncas, mc.nelecas))[::-1].tolist()
                       for c in wave.ci]
        diagnostics['active_natural_occupations'] = occupations
        if 'fragments' in s:
            from .fragment_charges import charge_operators
            charges, nuclei = charge_operators(wave, s['fragments'])
            diagnostics.update(fragment_charge_operators_e=charges, fragment_nuclear_charges_e=nuclei)
        if prior is None:
            active = slice(mc.ncore, mc.ncore+mc.ncas)
            projection = initial_mo[:, active].T @ mf.get_ovlp() @ mc.mo_coeff[:, active]
            reference = 'rohf' if spin else 'rhf'
            diagnostics[f'initial_active_{reference}_mo_1based'] = s.get('initial_active_mo_1based',
                                                                list(range(mc.ncore+1, mc.ncore+mc.ncas+1)))
            diagnostics[f'{reference}_to_casscf_active_subspace_squared_overlaps'] = (
                np.linalg.svd(projection, compute_uv=False)**2).tolist()
        return AdiabaticFrame(energies, None, gradients, nac, wave, diagnostics)
