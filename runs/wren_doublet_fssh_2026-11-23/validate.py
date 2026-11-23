"""Validate C2H4 dimer cation derivatives, state preparation and charge dynamics.

Run with Wren, NumPy and PySCF installed:
    python runs/wren_doublet_fssh_2026-11-23/validate.py build/dimer-hole-validation

Optional --trajectory-results COARSE FINE reads two already completed runs,
checking their trace hashes, model and XYZ hashes, configurations and end times.
"""
import argparse
from importlib.metadata import version
import json
from pathlib import Path

import numpy as np
from pyscf import fci, gto, scf

from wren.surface_hopping.casscf_input import read_casscf_model
from wren.surface_hopping.control import sha256
from wren.surface_hopping.pyscf_casscf import SACASSCF
from wren.surface_hopping.runner import run_control

EXAMPLES = Path(__file__).resolve().parents[2]/'examples/surface-hopping'
CONTROL = EXAMPLES/'ethene-hole.toml'
MODEL = EXAMPLES/'ethene-hole.json'


def direct_overlap(left, right):
    # 31 electrons, CAS(3,2), 14 closed-core orbitals. Both active alpha
    # orbitals are occupied, and the remaining beta electron occupies one.
    active = left.ncas
    nc = left.ncore
    if active != 2 or tuple(left.nelecas) != (2,1):
        raise ValueError('specific reference requires CAS(3,2)')
    s = left.mo.T @ gto.intor_cross('int1e_ovlp', left.mol, right.mol) @ right.mo
    alpha = np.linalg.det(s[:nc+active,:nc+active])
    beta = np.array([[np.linalg.det(s[np.ix_([*range(nc),nc+i],[*range(nc),nc+j])])
                      for j in range(active)] for i in range(active)])
    return alpha * left.ci[0].ravel() @ beta @ right.ci[1].ravel()


def quantum_check(settings, q):
    p = SACASSCF(dict(settings, use_etfs=False))
    center = p.evaluate_adiabatic(q)
    w = center.wavefunctions
    ia = int(np.argmax(abs(center.nac[0,1])))
    h = .002
    dq = np.zeros_like(q)
    dq[ia] = h
    plus = p.evaluate_adiabatic(q+dq, center)
    minus = p.evaluate_adiabatic(q-dq, center)
    fd = (plus.energies-minus.energies)/(2*h)
    nac = (direct_overlap(w, plus.wavefunctions) - direct_overlap(w, minus.wavefunctions))/(2*h)
    np.testing.assert_allclose(fd, center.gradients[:,ia], rtol=0, atol=5e-4)
    np.testing.assert_allclose(nac, center.nac[0,1,ia], rtol=0, atol=.02)
    # Independent full AO 1-RDM for a coherent (complex) state. Do not use
    # the previously constructed fragment transition-charge matrices here.
    c = np.array([.5, .5+np.sqrt(.5)*1j])
    combined = sum(z*ci for z,ci in zip(c,w.ci))
    d = sum(fci.direct_spin1.make_rdm1(z.copy(), w.ncas, w.nelecas) for z in (combined.real,combined.imag))
    core, active = w.mo[:,:w.ncore],w.mo[:,w.ncore:w.ncore+w.ncas]
    density = 2*core@core.T + active@d@active.T
    metric = w.mol.intor_symmetric('int1e_ovlp')
    eig, u = np.linalg.eigh(metric)
    sqrt_s = (u*np.sqrt(eig))@u.T
    atom_aos = w.mol.aoslice_by_atom()
    pops = np.diag(sqrt_s@density@sqrt_s)
    independent = {}
    for name, group in settings['fragments'].items():
        atoms = np.asarray(group)-1
        independent[name] = sum(w.mol.atom_charge(int(i))-pops[atom_aos[i,2]:atom_aos[i,3]].sum() for i in atoms)
    operators = center.diagnostics['fragment_charge_operators_e']
    charges = {k:np.vdot(c,np.asarray(op)@c).real for k,op in operators.items()}
    np.testing.assert_allclose(list(charges.values()),list(independent.values()),atol=1e-9)
    # Both ethylenes lie in x,y planes. Test actual active π character on the
    # four individual C 2p_z orbitals in the nonorthogonal AO metric.
    idx = [i for i,(_,el,shell,component) in enumerate(w.mol.ao_labels(fmt=False))
           if el=='C' and shell=='2p' and component=='z']
    mf = scf.ROHF(w.mol).run(conv_tol=1e-11)
    if not mf.converged:
        raise RuntimeError('orbital validation ROHF failed')
    def projected(mos):
        z = metric[idx]@mos
        return z.T@np.linalg.solve(metric[np.ix_(idx,idx)],z)
    return {'closed_core_orbitals':w.ncore,'active_electrons_alpha_beta':list(w.nelecas),
            'spin_squared':center.diagnostics['spin_squared'], 'starting_rohf_pi_weights':np.diag(projected(
                mf.mo_coeff[:,np.array(settings['initial_active_mo_1based'])-1])).tolist(),
            'optimized_active_pi_subspace_squared_overlaps':np.linalg.eigvalsh(projected(active)).tolist(),
            'coordinate':{'atom_1based':ia//3+1,'axis':'xyz'[ia%3]},'step_bohr':h,'etfs':False,
            'gradients_analytic_hartree_per_bohr':center.gradients[:,ia].tolist(),
            'gradients_fd_hartree_per_bohr':fd.tolist(),'nac_01_analytic_per_bohr':float(center.nac[0,1,ia]),
            'nac_01_fd_full_determinants_per_bohr':float(nac),
            'test_complex_amplitudes_real':c.real.tolist(),'test_complex_amplitudes_imag':c.imag.tolist(),
            'direct_AO_fragment_charges_e':independent,'transition_charge_expectations_e':charges}


def read_trajectory(directory, settings, sources, dt, steps):
    result = json.loads((directory/'result.json').read_text())
    trace = directory/'trajectories.jsonl'
    if result['status'] != 'completed' or sha256(trace)!=result['trajectory_sha256']:
        raise ValueError('only a completed trajectory with matching hash may be analyzed')
    c = result['config']
    if (result['quantum_settings']!=settings or any(result['provenance']['sources'][k]!=v for k,v in sources.items())
            or c['dt']!=dt or c['steps']!=steps or c['seed']!=2026 or c['charge_localization']!=['hole','A']
            or c['substeps']!=8 or c['trajectories']!=1 or c['decoherence']!='none' or any(c['momenta'])):
        raise ValueError('trajectory source or protocol mismatch')
    rows = [json.loads(line) for line in trace.read_text().splitlines()]
    if ([r['step'] for r in rows]!=list(range(steps+1))
            or [r['time_au'] for r in rows]!=[i*dt for i in range(steps+1)]):
        raise ValueError('sample times or trajectory end differ from controlled protocol')
    max_spin_error = max(abs(s2-.75) for r in rows for s2 in r['quantum']['spin_squared'])
    if max_spin_error>1e-6 or not all(r['quantum']['response_converged'] for r in rows):
        raise AssertionError('doublet spin or response convergence check failed')
    total_error = max(abs(sum(row['fragment_net_charges_e'].values())-1) for row in rows)
    total_electrons_error = max(abs(sum(row['fragment_electron_populations'].values())-31) for row in rows)
    if total_error>1e-8 or total_electrons_error>1e-8:
        raise AssertionError('total charge or electron count drifted')
    return {'dt_au':dt,'time_au':rows[-1]['time_au'], 'seed':c['seed'],'hops':result['hops'],
            'initial_localization':result['localization'],'initial_active_state':rows[0]['active_state'],
            'initial_energies_hartree':rows[0]['energies_hartree'],
            'initial_natural_occupations':rows[0]['quantum']['active_natural_occupations'],
            'initial_charge_operators_e':rows[0]['quantum']['fragment_charge_operators_e'],
            'initial_active_probabilities':result['initial_active_probabilities'],
            'max_charge_sum_error_e':total_error, 'max_electron_sum_error':total_electrons_error,
            'max_energy_drift_hartree':result['max_energy_drift_hartree'],
            'max_norm_error':result['max_norm_error'],'max_spin_squared_error':max_spin_error,
            'min_gap_hartree':result['minimum_state_gap_hartree'],'trajectory_sha256':result['trajectory_sha256'],
            'min_adjacent_overlap':result['minimum_adjacent_state_overlap'],
            'history':[{'t_au':r['time_au'],'charge_A_e':r['fragment_net_charges_e']['A'],
                        'charge_B_e':r['fragment_net_charges_e']['B'],'active_state':r['active_state']}
                       for r in rows]}


def run(out, previous):
    out.mkdir(parents=True,exist_ok=False)
    cfg,q,_,sources = read_casscf_model(MODEL)
    qc = quantum_check(cfg,np.array(q))
    (out/'derivatives.json').write_text(json.dumps(qc,indent=2,allow_nan=False)+'\n')
    if previous is None:
        fine = out/'half.toml'
        fine.write_text(CONTROL.read_text().replace('"ethene-hole.json"',json.dumps(str(MODEL)))
                        .replace('dt-au = 2.0','dt-au = 1.0').replace('max-steps = 6','max-steps = 12'))
        previous = [out/'coarse',out/'half']
        for control, target in zip([CONTROL,fine],previous):
            run_control(control,output=target)
    summarize(out, cfg, sources, qc, previous)


def summarize(out, cfg, sources, qc, previous):
    coarse = read_trajectory(previous[0],cfg,sources,2.,6)
    fine = read_trajectory(previous[1],cfg,sources,1.,12)
    report = {'system':'(C2H4)2 radical cation in vacuum, spin doublet, ROHF → two-state equal-weight SA-CASSCF(3,2)/STO-3G',
              'versions':{n:version(n) for n in ('pyscf','scipy','numpy')},
              'input_sha256':{Path(k).name:v for k,v in sources.items()}, 'quantum_check':qc,
              'coarse': coarse, 'fine':fine,'abs_final_charge_A_difference_e': abs(
                  coarse['history'][-1]['charge_A_e']-fine['history'][-1]['charge_A_e'])}
    (out/'validation.json').write_text(json.dumps(report, indent=2,allow_nan=False)+'\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output',type=Path)
    parser.add_argument('--trajectory-results',type=Path,nargs=2)
    args=parser.parse_args()
    run(args.output,args.trajectory_results)
