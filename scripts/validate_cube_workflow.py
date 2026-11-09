#!/usr/bin/env python3
"""Real CPU/GPU checks for cube input dispatch, DIPRO and vacuum four-point RE."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import time

import numpy as np
from universal_kmc.quantum.geometry import read_xyz
from universal_kmc.workflow.dimers import fingerprint
from universal_kmc.workflow.control import load_input, scf_config

ROOT = Path(__file__).resolve().parents[1]


def check_gradients(settings, directory, *, embedding=None):
    """Check neutral/cation/anion forces against actual central energy differences."""
    from pyscf.data.nist import BOHR
    from universal_kmc.quantum.coupling import host
    from universal_kmc.workflow.scf import make_mean_field, scf_record
    directory.mkdir(parents=True,exist_ok=True)
    config=scf_config(settings)
    results={}
    h=1e-3  # Bohr
    for charge in (0,1,-1):
        energies=[]
        analytic=None
        for offset in (0,h,-h):
            atoms=[('F',(-.75,0.,0.)),('F',(.75+offset*BOHR,0.,0.))]
            mf=make_mean_field(atoms,config,charge=charge,spin=0 if charge==0 else 1,
                               log=directory/f'q{charge}_{offset:+g}.log', embedding=embedding)
            try:
                mf.kernel()
                scf_record(mf,config.backend)
                lj_energy, lj_gradient = embedding.classical(atoms) if embedding is not None else (0., 0.)
                energies.append(float(mf.e_tot)+lj_energy)
                if offset==0:
                    analytic=float((host(mf.nuc_grad_method().kernel())+lj_gradient)[1,0])
            finally:
                mf.mol.stdout.close()
        finite=(energies[1]-energies[2])/(2*h)
        error=abs(analytic-finite)
        assert error<1e-4,(charge,analytic,finite,error)
        results[str(charge)]={'analytic_hartree_per_bohr':analytic,'finite_difference_hartree_per_bohr':finite,
                              'absolute_error':error,'step_bohr':h}
    return results


def run(control, output, *, resume=False, expected=0):
    args = [sys.executable,'-m','universal_kmc.workflow',str(control)]
    if resume:
        args.append('--resume')
    start = time.perf_counter()
    process = subprocess.run(args,cwd='/tmp',capture_output=True,text=True)
    output.write_text(process.stdout+'\nSTDERR\n'+process.stderr)
    if process.returncode != expected:
        raise AssertionError(f'{args}: exit {process.returncode}, expected {expected}; see {output}\n{process.stderr[-3000:]}')
    return time.perf_counter()-start


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',required=True,type=Path)
    parser.add_argument('--gpu',action='store_true')
    args = parser.parse_args()
    directory=args.output.resolve()
    directory.mkdir(parents=True,exist_ok=True)
    a,b = (read_xyz(ROOT/f'examples/quantum/ethylene_{s}.xyz') for s in ('a','b'))
    coordinates = [np.array([r for e,r in atoms]) for atoms in (a,b)]
    elements = [[e for e,r in atoms] for atoms in (a,b)]
    all_coords = np.concatenate(coordinates)
    monomers = [{'id':str(i+1),'elements':elements[i],
                 'coordinates_angstrom':(c-c.mean(axis=0)).tolist()} for i,c in enumerate(coordinates)]
    dimer = {'id':'1_2','mol_i':'1','mol_j':'2','elements':elements[0]+elements[1],
              'coordinates_angstrom':(all_coords-all_coords.mean(axis=0)).tolist(),
              'atom_count_a':len(a),'core_min_angstrom':3.5,'image':[0,0,0]}
    manifest={'schema':'cube-dimers-1','coordinate_unit':'Angstrom','monomers':monomers,'dimers':[dimer]}
    manifest['fingerprint']=fingerprint(manifest)
    (directory/'mini_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    (directory/'F2.xyz').write_text('2\nF2 starting geometry, vacuum\nF -.75 0 0\nF .75 0 0\n')
    results={}
    for backend in (('cpu','gpu') if args.gpu else ('cpu',)):
        here=directory/backend
        here.mkdir(exist_ok=True)
        engine='gpu-pyscf' if backend=='gpu' else 'pyscf'
        control=here/'input'
        control.write_text(f'''$two-mol
input = ../mini_manifest.json
output = two
method = b3lyp
basis-set = 6-31G(d,p)
transfer-integrals-method = DIPRO
$end
$one-mol
input = ../mini_manifest.json
output = one
qc-engine = {engine}
method = b3lyp
basis-set = 6-31G(d,p)
$end
$reorganization-energy
xyz = ../F2.xyz
output = reorganization
qc-engine = {engine}
method = b3lyp
basis-set = 6-31G(d,p)
re-method = 4-points
environment = vacuum
carriers = hole,electron
opt-maxsteps = 60
$end
''')
        elapsed=run(control,here/'console.log')
        transfer=json.loads((here/'two/1_2/transfer.json').read_text())
        reorg=json.loads((here/'reorganization/F2/reorganization.json').read_text())
        assert transfer['status']=='completed' and transfer['scf_executed'] is False
        assert reorg['status']=='completed' and reorg['environment']=='vacuum'
        assert all(v['reorganization_energy_ev']>0 for v in reorg['results'].values())
        assert all(v['charged_optimization']['converged'] for v in reorg['results'].values())
        if backend=='gpu':
            assert transfer['gpu_runtime']['projection_array_module']=='cupy'
            assert reorg['neutral_optimization']['gradient_class'].startswith('gpu4pyscf.')
            for result in reorg['results'].values():
                assert result['charged_optimization']['gradient_class'].startswith('gpu4pyscf.')
                assert all(p['backend']=='gpu' and p['device'] for p in result['single_points'].values())
        files=[*here.rglob('*.chk'),*here.rglob('*.trace.json'),here/'two/1_2/transfer.json']
        timestamps={str(p):p.stat().st_mtime_ns for p in files}
        resume_elapsed=run(control,here/'resume-console.log',resume=True)
        assert all(Path(p).stat().st_mtime_ns==stamp for p,stamp in timestamps.items())
        summary=json.loads((here/'cube-results/workflow.json').read_text())
        assert all(v['reused'] for v in summary['tasks']['one-mol']['outputs'].values())
        assert summary['tasks']['reorganization-energy']['outputs']['F2']['reused']
        results[backend]={'elapsed_seconds_observation':elapsed,'resume_seconds_observation':resume_elapsed,
            'couplings_ev':{c:transfer['results'][c]['methods']['dipro']['scalar_magnitude_ev'] for c in ('hole','electron')},
            'reorganization_ev':{c:v['reorganization_energy_ev'] for c,v in reorg['results'].items()},
            'neutral_optimization':reorg['neutral_optimization'],'resume_no_recomputation':'passed',
            'gradients':check_gradients(load_input(control).sections['reorganization-energy'],here/'gradient-check')}
    if args.gpu:
        for carrier in ('hole','electron'):
            np.testing.assert_allclose(results['cpu']['couplings_ev'][carrier],results['gpu']['couplings_ev'][carrier],atol=2e-5,rtol=2e-4)
            assert abs(results['cpu']['reorganization_ev'][carrier]-results['gpu']['reorganization_ev'][carrier])<5e-3
    # Force a genuinely unconverged geometry task; do not emit a successful RE.
    failure=directory/'failure'
    failure.mkdir(exist_ok=True)
    (failure/'input').write_text('''$reorganization-energy
xyz = ../F2.xyz
output = reorg
qc-engine = pyscf
method = b3lyp
basis-set = 6-31G(d,p)
re-method = 4-points
environment = vacuum
opt-maxsteps = 1
carriers = hole
$end
''')
    run(failure/'input',failure/'console.log',expected=1)
    failed=json.loads((failure/'reorg/F2/reorganization.json').read_text())
    assert failed['status']=='failed' and not failed['results']
    report={'status':'passed','backends':results,'failed_optimization_no_reorganization':'passed',
            'coupling_tolerance':{'atol_ev':2e-5,'rtol':2e-4},'reorganization_cpu_gpu_atol_ev':5e-3,
            'scope':'ethylene DIPRO and F2 vacuum four-point; timings are observations, not a general speedup benchmark'}
    (directory/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    main()
