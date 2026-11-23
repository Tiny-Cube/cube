"""Reproduce LiH derivative checks and short LiH/C2H4 timestep comparisons.

Run with PySCF, NumPy and Wren installed:
    python runs/wren_casscf_fssh_2026-11-23/validate.py build/casscf-fssh-validation

Writes complete trajectories and a compact summary to a new output directory.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from pyscf import gto, lib, scf

from wren.surface_hopping.casscf_input import read_casscf_model
from wren.surface_hopping.pyscf_casscf import SACASSCF
from wren.surface_hopping.runner import run_control


EXAMPLES = Path(__file__).resolve().parents[2] / 'examples/surface-hopping'


def check_lih():
    from pyscf import gto
    cfg, q, _, _ = read_casscf_model(EXAMPLES/'lih-casscf.json')
    model = SACASSCF(dict(cfg, use_etfs=False))
    q = np.array(q)
    center = model.evaluate_adiabatic(q)

    def overlap(right):
        # Direct full Slater determinant contractions independent of the provider.
        left = center.wavefunctions
        s = left.mo.T @ gto.intor_cross('int1e_ovlp', left.mol, right.mol) @ right.mo
        m = np.array([[np.linalg.det(s[np.ix_([0,i+1],[0,j+1])]) for j in range(2)] for i in range(2)])
        return np.einsum('ab,ac,bd,cd', left.ci[0], m, m, right.ci[1])

    reports = {}
    for name, coordinate in [('Li_z', 2), ('H_z', 5)]:
        dq = np.zeros(6)
        dq[coordinate] = .002
        plus = model.evaluate_adiabatic(q + dq, center)
        minus = model.evaluate_adiabatic(q - dq, center)
        grad = (plus.energies - minus.energies)/.004
        signed_nac = (overlap(plus.wavefunctions) - overlap(minus.wavefunctions))/.004
        reports[name] = {'gradient_analytic': center.gradients[:,coordinate].tolist(), 'gradient_fd': grad.tolist(),
                         'nac_01_analytic': center.nac[0,1,coordinate], 'nac_01_fd': signed_nac}
        np.testing.assert_allclose(grad, center.gradients[:,coordinate], atol=1e-5, rtol=0)
        np.testing.assert_allclose(signed_nac, center.nac[0,1,coordinate], atol=3e-4, rtol=0)
    return {'step_bohr': .002, 'etfs': False, 'derivatives': reports}


def ethene_orbitals():
    cfg, q, _, _ = read_casscf_model(EXAMPLES/'ethene-casscf.json')
    model = SACASSCF(cfg)
    frame = model.evaluate_adiabatic(q)
    coords = np.reshape(q, (-1,3))
    mol = frame.wavefunctions.mol
    mf = scf.RHF(mol).run(conv_tol=1e-11)
    if not mf.converged:
        raise RuntimeError('reference RHF did not converge')
    # Minimal basis only: a local carbon p_pi AO follows the normal of the CH2 plane.
    # Orient both normals to the same side before inspecting relative MO phases.
    normals = []
    carbon = np.zeros((mol.nao_nr(),2))
    for c, h1, h2 in [(0,2,3),(1,4,5)]:
        normal = np.cross(coords[h1]-coords[c], coords[h2]-coords[c])
        normal /= np.linalg.norm(normal)
        if normals and normal @ normals[0] < 0:
            normal *= -1
        normals.append(normal)
        for ao, (atom, _, shell, component) in enumerate(mol.ao_labels(fmt=False)):
            if atom == c and shell == '2p':
                carbon[ao,c] = normal['xyz'.index(component)]
    metric = mf.get_ovlp()
    pi_metric = carbon.T @ metric @ carbon

    def projections(orbitals):
        cross = carbon.T @ metric @ orbitals
        coefficients = np.linalg.solve(pi_metric, cross)
        return cross.T @ coefficients, coefficients

    rhf_indices = np.array(cfg['initial_active_mo_1based']) - 1
    pi_projector, coeff = projections(mf.mo_coeff[:,rhf_indices])
    wfn = frame.wavefunctions
    active, _ = projections(wfn.mo[:,wfn.ncore:wfn.ncore+wfn.ncas])
    return {'reference': 'AO S-metric projection onto one local 2p normal per CH2, aligned normals; STO-3G',
            'angle_between_ch2_plane_normals_deg':np.degrees(np.arccos(np.clip(normals[0]@normals[1],-1.,1.))),
            'initial_rhf_mo_1based': list(cfg['initial_active_mo_1based']),
            'initial_rhf_mo_energy_hartree': mf.mo_energy[rhf_indices].tolist(),
            'initial_rhf_local_pi_weights': np.diag(pi_projector).tolist(),
            'relative_sign_of_local_carbon_p_components': np.sign(np.prod(coeff,axis=0)).tolist(),
            'optimized_active_pi_subspace_squared_overlaps': np.linalg.eigvalsh(active).tolist(),
            'natural_occupations_by_state': frame.diagnostics['active_natural_occupations']}


def trajectories(name, directory):
    config = EXAMPLES/f'{name}-casscf.toml'
    from wren.surface_hopping.control import load_control
    settings, _, _, _ = load_control(config)
    fine_control = directory/f'{name}-half-step.toml'
    # Keep the same starting structure, seed, end time, substeps per nuclear step.
    fine_control.write_text(config.read_text()
        .replace(f'"{name}-casscf.json"', json.dumps(str(config.with_suffix('.json'))))
        .replace(f'dt-au = {settings.dt}', f'dt-au = {settings.dt / 2}')
        .replace(f'max-steps = {settings.steps}', f'max-steps = {settings.steps * 2}'))
    records = []
    for dt, source, label in [(settings.dt, config, 'standard'), (settings.dt / 2, fine_control, 'half')]:
        result = run_control(source, output=directory/f'{name}-{label}')
        initial = json.loads(Path(result['trajectories']).read_text().splitlines()[0])
        end = result['final_states'][0]
        records.append({'dt_au': dt, 'steps': end['step'], 'time_au': end['time_au'],
                        'energy_initial_hartree': initial['energies_hartree'],
                        'initial_nac_01_norm_per_bohr': np.linalg.norm(initial['derivative_couplings_per_bohr'][0][1]),
                        'spin_squared_initial': initial['quantum']['spin_squared'],
                        'initial_rhf_to_casscf_active_squared_overlaps':
                            initial['quantum']['rhf_to_casscf_active_subspace_squared_overlaps'],
                        'end_coordinates_bohr':end['q_bohr'], 'end_populations': end['adiabatic_populations'],
                        **{k: result[k] for k in ('hops','max_norm_error','max_energy_drift_hartree','minimum_state_gap_hartree',
                              'minimum_adjacent_state_overlap','max_probability_sum')}})
    if records[0]['time_au'] != records[1]['time_au']:
        raise AssertionError('comparison requires identical end time')
    return {'traces': records, 'max_final_coordinate_difference_bohr': float(np.max(abs(np.subtract(
            records[0]['end_coordinates_bohr'], records[1]['end_coordinates_bohr'])))),
            'max_final_population_difference':float(np.max(abs(np.subtract(records[0]['end_populations'],
                                                                           records[1]['end_populations']))))}


if __name__ == '__main__':
    from importlib.metadata import version
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    output = parser.parse_args().output
    output.mkdir(parents=True, exist_ok=False)
    lib.num_threads(1)
    result = {'scope': 'isolated LiH and C2H4, equal-weight singlet CAS(2,2)/STO-3G; short deterministic-start tests',
              'versions': {p:version(p) for p in ('numpy','scipy','pyscf')},
              'lih_finite_differences':check_lih(), 'ethene_pi_character':ethene_orbitals()}
    for name in ('lih','ethene'):
        result[name] = trajectories(name, output)
        (output/'validation.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
