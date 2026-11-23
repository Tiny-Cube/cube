"""Löwdin fragment charges projected into a real multistate CAS wavefunction space."""
import numpy as np


def charge_operators(wave, fragments):
    from pyscf import fci
    mol = wave.mol
    s = mol.intor_symmetric('int1e_ovlp')
    eigenvalues, u = np.linalg.eigh(s)
    if eigenvalues[0] <= 1e-10 * eigenvalues[-1]:
        raise ValueError('AO overlap nearly linearly dependent; fragment charge cannot be resolved reliably')
    lowdin = ((u * np.sqrt(eigenvalues)) @ u.T) @ wave.mo[:, :wave.ncore+wave.ncas]
    slices = mol.aoslice_by_atom()
    states = len(wave.ci)
    identity = np.eye(states)
    charges, nuclei = {}, {}
    for name, group in fragments.items():
        # P on this fragment is diagonal in the symmetric orthogonal AO basis.
        atoms = np.asarray(group) - 1
        rows = np.concatenate([np.arange(slices[i,2],slices[i,3]) for i in atoms])
        local = lowdin[rows]
        core_population = 2 * np.sum(local[:,:wave.ncore]**2)
        active_projector = local[:,wave.ncore:].T @ local[:,wave.ncore:]
        nuclei[name] = float(sum(mol.atom_charge(int(i)) for i in atoms))
        op = (nuclei[name] - core_population)*identity
        for i, bra in enumerate(wave.ci):
            for j in range(i+1):
                density = fci.direct_spin1.trans_rdm1(bra, wave.ci[j], wave.ncas, wave.nelecas)
                value = float(np.einsum('pq,qp->', active_projector, density))
                op[i,j] -= value
                op[j,i] = op[i,j]
        charges[name] = op.tolist()
    if not np.allclose(sum(np.asarray(v) for v in charges.values()), mol.charge*identity, atol=1e-8, rtol=0):
        raise RuntimeError('fragment charge operators do not sum to the molecular charge times identity')
    return charges, nuclei


def localize_charge(operators, name, carrier):
    """Extremize a fragment charge over the computed adiabatic-state space."""
    charges, vectors = np.linalg.eigh(operators[name])
    if carrier == 'hole':
        gap, selected = charges[-1]-charges[-2], -1
    else:
        gap, selected = charges[1]-charges[0], 0
    if not np.isfinite(charges).all() or gap <= 1e-8:
        raise ValueError('localized initial charge has no unique extremum in the supplied electronic subspace')
    c = vectors[:,selected]
    c *= np.sign(c[np.argmax(abs(c))])
    return c, float(charges[selected])


def fragment_observables(c, operators, nuclei):
    charges = {name: float(np.vdot(c, np.asarray(op) @ c).real) for name, op in operators.items()}
    return {'fragment_net_charges_e': charges,
            'fragment_electron_populations': {name: nuclei[name]-v for name, v in charges.items()}}
