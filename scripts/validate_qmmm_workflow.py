#!/usr/bin/env python3
"""Real fixed-MM four-point, gradient, CPU/GPU and resume checks."""
import argparse
import json
from pathlib import Path
import sys
from dataclasses import replace

import numpy as np
from pyscf.data.nist import HARTREE2EV
from pyscf.scf import chkfile

from universal_kmc.workflow.control import load_input
from universal_kmc.workflow.dimers import file_hash
from universal_kmc.workflow.embedding import read_charges
from validate_cube_workflow import check_gradients, run


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--gpu', action='store_true')
    parser.add_argument('--lj', action='store_true', help='include a small test LJ12-6 environment')
    args = parser.parse_args()
    output = args.output.resolve()
    if output.exists():
        raise ValueError('choose a new output directory for this validation')
    output.mkdir(parents=True)
    (output/'F2.xyz').write_text('2\nF2 QM region; Angstrom\nF -.75 0 0\nF .75 0 0\n')
    charges = '# x y z q; Angstrom/e; fixed neutral quadrupolar test environment\n3 0 0 0.02\n-3 0 0 0.02\n0 3 0 -0.02\n0 -3 0 -0.02\n'
    pc = output/'mm.pc'
    pc.write_text(charges)
    lj_controls = ''
    if args.lj:
        from universal_kmc.quantum.geometry import read_xyz
        from universal_kmc.workflow.lj import LennardJones, geometry_key, sites_key
        environment = read_charges(pc)
        epsilon, sigma = .1, .24  # Synthetic numerical test parameters: kJ/mol, nm.
        model = LennardJones((0, 0), (0, 0, 0, 0), ((4*epsilon*sigma**6,),), ((4*epsilon*sigma**12,),),
                             {'source': 'synthetic_gradient_test_not_material_force_field',
                              'qm_initial_geometry_fingerprint': geometry_key(read_xyz(output/'F2.xyz')),
                              'mm_sites_fingerprint': sites_key(environment)})
        (output/'lj.json').write_text(json.dumps(model.payload(), indent=2)+'\n')
        lj_controls = 'mm-vdw = file\nmm-vdw-file = ../lj.json\n'
    results = {}
    for backend in (('cpu', 'gpu') if args.gpu else ('cpu',)):
        here = output/backend
        here.mkdir()
        control = here/'input'
        control.write_text(f'''$reorganization-energy
xyz = ../F2.xyz
output = reorganization
qc-engine = {'gpu-pyscf' if backend == 'gpu' else 'pyscf'}
method = b3lyp
basis-set = 6-31G(d,p)
re-method = 4-points
environment = qmmm
mm-charges = ../mm.pc
mm-frozen = true
{lj_controls}carriers = hole,electron
opt-maxsteps = 60
$end
''')
        elapsed = run(control, here/'console.log')
        root = here/'reorganization/F2'
        record = json.loads((root/'reorganization.json').read_text())
        assert record['status'] == 'completed' and record['environment'] == 'qmmm'
        env = read_charges(pc)
        if args.lj:
            from universal_kmc.quantum.geometry import read_xyz
            from universal_kmc.workflow.lj import read_lj_file
            env = replace(env, vdw=read_lj_file(output/'lj.json', read_xyz(output/'F2.xyz'), env))
        snapshot = json.loads((root/'mm-environment.json').read_text())
        assert snapshot['fingerprint'] == env.fingerprint
        np.testing.assert_array_equal(snapshot['coordinates_angstrom'], env.coordinates_angstrom)
        np.testing.assert_array_equal(snapshot['charges_e'], env.charges_e)
        for state in ('neutral', 'cation', 'anion'):
            trace = json.loads((root/f'{state}.trace.json').read_text())
            assert trace and all(len(step['coordinates_angstrom']) == 2 for step in trace)
            for step in trace:
                scf = step['scf']
                assert scf['mm_frozen'] and scf['embedding_fingerprint'] == env.fingerprint
                assert scf['backend'] == backend and scf['converged']
                if backend == 'gpu':
                    assert scf['device']
        for result in record['results'].values():
            assert result['charged_optimization']['converged']
            assert result['charged_optimization']['qm_mm_gradient_verified']
            assert result['reorganization_energy_ev'] > 0
            energies = result['energies_hartree']
            expected = (energies['Eq_R0']-energies['Eq_Rq']+energies['E0_Rq']-energies['E0_R0'])*HARTREE2EV
            np.testing.assert_allclose(expected, result['reorganization_energy_ev'], atol=1e-12)
            for point in result['single_points'].values():
                assert point['embedding_fingerprint'] == env.fingerprint
                marker = chkfile.load(point['checkpoint'], 'cube_reorganization')
                np.testing.assert_array_equal(marker['mm_coordinates_angstrom'], env.coordinates_angstrom)
                np.testing.assert_array_equal(marker['mm_charges_e'], env.charges_e)
                assert chkfile.load(point['checkpoint'], 'cube_transfer') is None
        hashes = {p: file_hash(p) for p in root.rglob('*') if p.is_file()}
        stamps = {p: p.stat().st_mtime_ns for p in hashes}
        resume_elapsed = run(control, here/'resume.log', resume=True)
        assert all(file_hash(p) == h and p.stat().st_mtime_ns == stamps[p] for p, h in hashes.items())
        gradient = check_gradients(load_input(control).sections['reorganization-energy'], here/'gradient-check', embedding=env)
        results[backend] = {'reorganization_ev': {c: v['reorganization_energy_ev'] for c, v in record['results'].items()},
                            'elapsed_seconds_observation': elapsed, 'resume_seconds_observation': resume_elapsed,
                            'gradients': gradient, 'mm_frozen_checked_every_evaluation': True,
                            'same_mm_at_all_four_points': True, 'resume_no_recomputation': True,
                            'neutral_optimization': record['neutral_optimization']}
    if args.gpu:
        for carrier in ('hole', 'electron'):
            assert abs(results['cpu']['reorganization_ev'][carrier]-results['gpu']['reorganization_ev'][carrier]) < 5e-3
    # Changing MM charges must never reuse an old Hamiltonian or overwrite results.
    for backend in results:
        control = output/backend/'input'
        record_path = output/backend/'reorganization/F2/reorganization.json'
        original_hash = file_hash(record_path)
        pc.write_text(charges.replace('0.02', '0.03'))
        try:
            run(control, output/backend/'changed-mm-resume.log', resume=True, expected=1)
            assert file_hash(record_path) == original_hash
        finally:
            pc.write_text(charges)
    failure = output/'failure'
    failure.mkdir()
    source = (output/'cpu/input').read_text().replace('opt-maxsteps = 60', 'opt-maxsteps = 1')
    (failure/'input').write_text(source)
    run(failure/'input', failure/'console.log', expected=1)
    failed = json.loads((failure/'reorganization/F2/reorganization.json').read_text())
    assert failed['status'] == 'failed' and not failed['results']
    summary = {'status': 'passed', 'scope': ('F2 frozen finite point-charge MM plus LJ12-6' if args.lj else
                                            'F2 frozen finite point-charge MM electrostatic embedding; no QM-MM vdw or MM relaxation'),
               'backends': results, 'mm_edit_invalidates_resume': True,
               'failed_optimization_no_reorganization': True, 'cpu_gpu_lambda_tolerance_ev': 5e-3}
    (output/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
