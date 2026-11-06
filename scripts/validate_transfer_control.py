#!/usr/bin/env python3
"""Validate the two-stage XYZ -> checkpoint -> control-file workflow."""
import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys
from unittest.mock import patch

import numpy as np
from pyscf import gto
from pyscf.scf import chkfile

from universal_kmc.quantum.checkpoints import calculate_from_checkpoints
from universal_kmc.quantum.control import load_control
from universal_kmc.quantum.transfer import SCFConfig, calculate_transfer, _versions

ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = ROOT/'examples/quantum'


def command(args, *, expected=0):
    process = subprocess.run([sys.executable, *map(str, args)], cwd='/tmp', capture_output=True, text=True)
    if process.returncode != expected:
        raise AssertionError(f'command failed ({process.returncode}, expected {expected}): {args}\n{process.stdout}\n{process.stderr}')
    return {'returncode': process.returncode, 'stdout': process.stdout, 'stderr': process.stderr}


def write_control(path, inputs, backend='cpu', *, equivalent=True, assume=False, atoms=None, orbitals=''):
    content = f'''schema_version = 1
[input]
fragment_a = {json.dumps(str(inputs[0]))}
fragment_b = {json.dumps(str(inputs[1]))}
dimer = {json.dumps(str(inputs[2]))}
[calculation]
backend = {json.dumps(backend)}
methods = ["dipro", "esid", "orbital_gmh"]
carriers = ["hole", "electron"]
assume_equivalent = {str(equivalent).lower()}
assume_converged = {str(assume).lower()}
'''
    if atoms:
        content += f'[fragments]\na_atoms = {json.dumps(atoms[0])}\nb_atoms = {json.dumps(atoms[1])}\n'
    content += orbitals
    content += f'[output]\nfile = {json.dumps(path.stem+".json")}\n'
    path.write_text(content)
    return path


def couplings(result):
    return np.array([result['results'][carrier]['methods'][method]['scalar_magnitude_ev']
                     for carrier in ('hole', 'electron') for method in ('dipro', 'esid', 'orbital_gmh')])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--gpu', action='store_true')
    parser.add_argument('--reuse', action='store_true', help='reuse this directory\'s existing successful stage-1 and direct-reference results')
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    backends = ('cpu', 'gpu') if args.gpu else ('cpu',)
    report = {'versions': _versions(), 'stage1': {}, 'postprocess': {}, 'rejections': {}}
    xyz = [EXAMPLE/f'ethylene_{name}.xyz' for name in ('a', 'b', 'ab')]
    for backend in backends:
        if args.reuse:
            for name in ('a', 'b', 'ab'):
                assert json.loads((output/backend/f'ethylene_{name}.status.json').read_text())['status'] == 'completed'
            report['stage1'][backend] = {'reused_existing_converged_outputs': True}
        else:
            report['stage1'][backend] = command([ROOT/'scripts/run_pyscf.py', *xyz,
                '--output-dir', output/backend, '--backend', backend])
    if args.reuse:
        reference = json.loads((output/'direct_reference.json').read_text())
    else:
        reference = calculate_transfer(xyz[0], xyz[1], config=SCFConfig(),
                                       methods=('dipro', 'esid', 'orbital_gmh'), assume_equivalent=True)
        (output/'direct_reference.json').write_text(json.dumps(reference, indent=2)+'\n')
    reference_values = couplings(reference)
    for source_backend in backends:
        # Relative checkpoint paths must resolve against the control, not /tmp.
        inputs = [Path(source_backend)/f'ethylene_{name}.chk' for name in ('a', 'b', 'ab')]
        for backend in backends:
            control = write_control(output/f'{source_backend}-source-{backend}-projection.toml', inputs, backend)
            log = command([ROOT/'scripts/calculate_transfer.py', control])
            result = json.loads(control.with_suffix('.json').read_text())
            if result['status'] != 'completed' or result['scf_executed'] is not False:
                raise AssertionError('saved-result workflow reported SCF or did not complete')
            values = couplings(result)
            np.testing.assert_allclose(values, reference_values, atol=2e-5, rtol=2e-4)
            report['postprocess'][control.stem] = {'maximum_coupling_difference_ev': float(np.max(abs(values-reference_values))),
                                                   'device': result['device'], **log}
    inputs = [output/'cpu'/f'ethylene_{name}.chk' for name in ('a', 'b', 'ab')]
    control = write_control(output/'no-scf.toml', inputs)
    with patch('pyscf.scf.hf.SCF.kernel', side_effect=AssertionError('postprocessing must not execute SCF')):
        result = calculate_from_checkpoints(load_control(control))
    np.testing.assert_allclose(couplings(result), reference_values, atol=2e-5, rtol=2e-4)
    report['scf_kernel_blocked_during_postprocessing'] = True

    # Permute dimer atoms and all corresponding AO rows; no fresh SCF.
    mol, data = chkfile.load_scf(str(inputs[2]))
    atom_order = [i for pair in zip(range(6), range(6, 12)) for i in pair]
    permuted_mol = gto.M(atom=[(mol.atom_symbol(i), mol.atom_coord(i)) for i in atom_order],
                         unit='Bohr', basis='def2-svp', verbose=0)
    order = [i for atom in atom_order for i in range(*map(int, mol.aoslice_by_atom()[atom, 2:4]))]
    permuted = output/'permuted.chk'
    chkfile.dump_scf(permuted_mol, str(permuted), data['e_tot'], data['mo_energy'], data['mo_coeff'][order], data['mo_occ'])
    metadata = chkfile.load(str(inputs[2]), 'cube_transfer')
    metadata['fock_hartree'] = metadata['fock_hartree'][np.ix_(order, order)]
    chkfile.save(str(permuted), 'cube_transfer', metadata)
    control = write_control(output/'permuted.toml', [*inputs[:2], permuted],
                            atoms=(list(range(1, 13, 2)), list(range(2, 13, 2))))
    command([ROOT/'scripts/calculate_transfer.py', control])
    permuted_result = json.loads(control.with_suffix('.json').read_text())
    np.testing.assert_allclose(couplings(permuted_result), couplings(result), atol=1e-11, rtol=0)
    report['permuted_atom_mapping'] = 'passed'

    # One-based explicit MOs must select exactly the automatically chosen states.
    orbital_table = '\n'.join(f'[orbitals.{carrier}]\na = [8]\nb = [8]' if carrier == 'hole'
                              else f'[orbitals.{carrier}]\na = [9]\nb = [9]' for carrier in ('hole', 'electron'))+'\n'
    control = write_control(output/'explicit-orbitals.toml', inputs, orbitals=orbital_table)
    command([ROOT/'scripts/calculate_transfer.py', control])
    explicit_result = json.loads(control.with_suffix('.json').read_text())
    np.testing.assert_allclose(couplings(explicit_result), couplings(result), atol=1e-12, rtol=0)
    report['one_based_explicit_orbitals'] = 'passed'

    import h5py
    standard = output/'native.chk'
    shutil.copyfile(inputs[0], standard)
    with h5py.File(standard, 'a') as file:
        del file['cube_transfer']
    control = write_control(output/'native-rejected.toml', [standard, *inputs[1:]])
    report['rejections']['native_without_convergence_confirmation'] = command([ROOT/'scripts/calculate_transfer.py', control], expected=1)
    control = write_control(output/'native-confirmed.toml', [standard, *inputs[1:]], assume=True)
    command([ROOT/'scripts/calculate_transfer.py', control])
    native = json.loads(control.with_suffix('.json').read_text())
    np.testing.assert_allclose(couplings(native), couplings(result), atol=1e-12, rtol=0)
    assert native['scf']['fragment_a']['convergence'] == 'user_assumed_from_external_log'
    assert native['model_consistency'] == 'not_available_in_standard_checkpoint'
    report['standard_native_checkpoint'] = 'passed_with_explicit_confirmation'
    native_dimer = output/'native-ab.chk'
    shutil.copyfile(inputs[2], native_dimer)
    with h5py.File(native_dimer, 'a') as file:
        del file['cube_transfer']
    control = write_control(output/'native-dimer-confirmed.toml', [*inputs[:2], native_dimer], assume=True)
    command([ROOT/'scripts/calculate_transfer.py', control])
    native_dimer_result = json.loads(control.with_suffix('.json').read_text())
    np.testing.assert_allclose(couplings(native_dimer_result), couplings(result), atol=2e-5, rtol=2e-4)
    assert native_dimer_result['operator_source'] == 'complete_saved_dimer_MO_eigenoperator'
    report['native_dimer_eigenoperator_maximum_difference_ev'] = float(np.max(abs(couplings(native_dimer_result)-couplings(result))))
    for kind in ('unconverged', 'bad_metric', 'truncated', 'wrong_frame', 'wrong_model', 'bad_fock'):
        bad = output/f'{kind}.chk'
        shutil.copyfile(inputs[0], bad)
        if kind == 'wrong_frame':
            m = chkfile.load_mol(str(bad))
            m.set_geom_(m.atom_coords()+[.1, 0., 0.], unit='Bohr')
            chkfile.save_mol(m, str(bad))
        else:
            with h5py.File(bad, 'a') as file:
                if kind == 'unconverged':
                    file['cube_transfer/converged'][...] = False
                elif kind == 'bad_metric':
                    file['scf/mo_coeff'][0, 0] += .1
                elif kind == 'truncated':
                    c = file['scf/mo_coeff'][()][:, :-1]
                    del file['scf/mo_coeff']
                    file['scf/mo_coeff'] = c
                elif kind == 'wrong_model':
                    del file['cube_transfer/model']
                    file['cube_transfer/model'] = 'hf'
                elif kind == 'bad_fock':
                    file['cube_transfer/fock_hartree'][0, 0] += .1
        control = write_control(output/f'reject-{kind}.toml', [bad, *inputs[1:]], assume=True)
        control.with_suffix('.json').write_text('{"status": "completed", "stale": true}\n')
        log = command([ROOT/'scripts/calculate_transfer.py', control], expected=1)
        failure = json.loads(control.with_suffix('.json').read_text())
        assert failure['status'] == 'failed'
        report['rejections'][kind] = {'error': failure['error'], **log}
    control = write_control(output/'partial-esid.toml', inputs, equivalent=False)
    command([ROOT/'scripts/calculate_transfer.py', control], expected=2)
    assert json.loads(control.with_suffix('.json').read_text())['status'] == 'partial'
    report['inapplicable_esid_exit_2'] = 'passed'
    # Nonconvergence must remove stale checkpoint output from stage 1.
    failed_dir = output/'nonconverged'
    failed_dir.mkdir(exist_ok=True)
    stale = failed_dir/'ethylene_a.chk'
    shutil.copyfile(inputs[0], stale)
    report['rejections']['stage1_nonconvergence'] = command([ROOT/'scripts/run_pyscf.py', xyz[0],
        '--output-dir', failed_dir, '--max-cycle', '1'], expected=1)
    assert not stale.exists()
    assert json.loads((failed_dir/'ethylene_a.status.json').read_text())['status'] == 'failed'
    report['status'] = 'passed'
    report['reference_couplings_ev_order'] = 'hole dipro/esid/orbital_gmh; electron dipro/esid/orbital_gmh'
    report['reference_couplings_ev'] = reference_values.tolist()
    (output/'summary.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
