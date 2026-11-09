#!/usr/bin/env python3
"""Compare G16 and PySCF energy/gradient providers with the same frozen MM+LJ.

This benchmark is separate from cube's qc-engine controls. Gaussian provides
QM energy/gradients with fixed external charges; cube supplies the same TPR LJ
and geomeTRIC optimizer used for the PySCF provider. This is not native ONIOM.
Gaussian conventional integrals and GPU4PySCF density fitting remain different
algorithms; the report records these settings and numerical discrepancies.
"""
import argparse
from dataclasses import replace
import json
import os
from pathlib import Path
import re
import subprocess
import time

import numpy as np

from universal_kmc.quantum.config import SCFConfig
from universal_kmc.quantum.geometry import read_xyz
from universal_kmc.workflow.dimers import file_hash, xyz_text
from universal_kmc.workflow.embedding import read_charges, snapshot
from universal_kmc.workflow.lj import read_lj_file
from universal_kmc.workflow.reorganization import four_point_terms


def write_json(path, value):
    path = Path(path)
    temporary = path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')
    temporary.replace(path)


def fchk_field(text, name):
    match = re.search(r'^'+re.escape(name)+r'\s+[IR]\s+(?:N=\s*(\d+)|(\S+))\s*$', text, re.M)
    if match is None:
        raise ValueError(f'Gaussian formatted checkpoint lacks {name}')
    if match[2] is not None:
        return float(match[2].replace('D','E'))
    count = int(match[1])
    tokens = text[match.end():].split()[:count]
    if len(tokens) != count:
        raise ValueError(f'incomplete Gaussian {name}')
    try:
        result = np.array([float(s.replace('D','E')) for s in tokens])
    except ValueError as exc:
        raise ValueError(f'incomplete Gaussian {name}') from exc
    if not np.isfinite(result).all():
        raise ValueError(f'nonfinite Gaussian {name}')
    return result


def gaussian_result(log, formatted_checkpoint, count, *, gradient):
    if 'Normal termination of Gaussian' not in log or 'Error termination' in log:
        raise RuntimeError('Gaussian did not terminate normally')
    scf = re.findall(r'SCF Done:\s+E\([^)]*\)\s*=\s*([-\d.DE+]+).*?after\s+(\d+)\s+cycles', log)
    self_energy = re.findall(r'Self energy of the charges\s*=\s*([-\d.DE+]+)', log)
    if not scf or not self_energy:
        raise RuntimeError('Gaussian lacks SCF convergence or external-charge self energy')
    energy = fchk_field(formatted_checkpoint, 'Total Energy')
    constant = float(self_energy[-1].replace('D','E'))
    if not np.isfinite([energy,constant]).all():
        raise RuntimeError('nonfinite Gaussian energy')
    result = {'electronic_energy_hartree': float(energy-constant),
              'g16_raw_energy_hartree': energy, 'g16_mm_self_energy_hartree': constant,
              'scf_cycles': int(scf[-1][1]),
              'basis_functions': int(fchk_field(formatted_checkpoint, 'Number of basis functions'))}
    if gradient:
        g = fchk_field(formatted_checkpoint, 'Cartesian Gradient')
        if g.shape != (count*3,):
            raise RuntimeError('Gaussian QM gradient dimension mismatch')
        result['gradient_hartree_per_bohr'] = g.reshape(count,3).tolist()
    return result


class Calculator:
    def __init__(self, atoms, charge, environment, args, directory):
        begin=time.perf_counter()
        self.args, self.charge, self.environment = args, charge, environment
        self.directory = Path(directory)
        self.directory.mkdir(parents=True, exist_ok=False)
        self.calls = []
        self.scanner = None
        if args.engine != 'g16':
            from universal_kmc.workflow.scf import make_mean_field
            config = SCFConfig(backend=args.engine, basis=args.basis, xc=args.method,
                               auxbasis='def2-universal-jkfit',threads=args.threads,
                               conv_tol=args.conv_tol,max_cycle=100,max_memory_mb=args.memory_mb)
            self.mf = make_mean_field(atoms, config, charge=charge, spin=int(charge!=0),
                                     log=self.directory/'electronic.log',embedding=environment)
            if hasattr(self.mf,'grids'):
                self.mf.grids.atom_grid=(99,590)
            self.scanner=self.mf.nuc_grad_method().as_scanner()
        self.initialization_seconds=time.perf_counter()-begin

    def evaluate(self, atoms, *, gradient=True):
        from pyscf import gto
        self.environment.validate_qm(atoms)
        begin = time.perf_counter()
        if self.args.engine == 'g16':
            result = self._gaussian(atoms,gradient)
        else:
            from universal_kmc.quantum.coupling import host
            from universal_kmc.workflow.scf import scf_record
            mol=gto.M(atom=atoms,basis=self.args.basis,charge=self.charge,spin=int(self.charge!=0),
                      unit='Angstrom',cart=False,verbose=0)
            if gradient:
                energy,g=self.scanner(mol)
            else:
                energy=self.scanner.base(mol)
            record=scf_record(self.scanner.base,self.args.engine)
            self.environment.assert_attached(self.scanner.base)
            result={'electronic_energy_hartree':float(energy),'device':record['device'],
                    'basis_functions':mol.nao,'scf_converged':True}
            if gradient:
                result['gradient_hartree_per_bohr']=np.asarray(host(g),dtype=float).tolist()
        electronic_seconds=time.perf_counter()-begin
        lj_begin=time.perf_counter()
        energy_lj,g_lj=self.environment.classical(atoms)
        result.update(qm_mm_vdw_energy_hartree=energy_lj,
                      energy_hartree=result['electronic_energy_hartree']+energy_lj,
                      electronic_wall_seconds=electronic_seconds,
                      lj_wall_seconds=time.perf_counter()-lj_begin,
                      charge=self.charge,gradient_requested=gradient)
        if gradient:
            g=np.array(result['gradient_hartree_per_bohr'])+g_lj
            if g.shape!=(len(atoms),3) or not np.isfinite(g).all():
                raise RuntimeError('invalid total QM/MM gradient')
            result['gradient_hartree_per_bohr']=g.tolist()
            result.update(gradient_rms=float(np.sqrt(np.mean(g*g))),gradient_max=float(abs(g).max()))
        result['elapsed_seconds']=time.perf_counter()-begin
        result['initialization_seconds']=self.initialization_seconds if not self.calls else 0.
        result['including_initialization_seconds']=result['elapsed_seconds']+result['initialization_seconds']
        self.calls.append(result)
        write_json(self.directory/'evaluations.json',self.calls)
        return result

    def _gaussian(self, atoms, gradient):
        n=len(self.calls)+1
        prefix=self.directory/f'eval{n:04d}'
        checkpoint=self.directory/'state.chk'
        scratch=self.directory/'scratch'
        scratch.mkdir(exist_ok=True)
        guess=' Guess=Read' if checkpoint.is_file() else ''
        integral='' if self.args.method.lower()=='hf' else ' Integral=UltraFine'
        route=(f'#p {self.args.method}/{self.args.basis} '+('Force' if gradient else 'SP')+
               f' Charge NoSymm 5D 7F SCF=(XQC,Conver=10,MaxCycle=100){integral}{guess}')
        text=(f'%chk={checkpoint}\n%nprocshared={self.args.threads}\n%mem={self.args.memory_mb}MB\n'
              f'{route}\n\nCentral QM; frozen MM charges; same external cube LJ added separately\n\n'
              f'{self.charge} {1 if self.charge==0 else 2}\n')
        # All coordinates must be floating-point tokens: Gaussian distinguishes optional flags.
        text+=''.join(f'{e} {r[0]:.15f} {r[1]:.15f} {r[2]:.15f}\n' for e,r in atoms)+'\n'
        text+=''.join(f'{r[0]:.15f} {r[1]:.15f} {r[2]:.15f} {q:.15f}\n'
                      for r,q in zip(self.environment.coordinates_angstrom,self.environment.charges_e))+'\n'
        prefix.with_suffix('.gjf').write_text(text)
        binary=Path(self.args.g16).resolve()
        env=dict(os.environ,GAUSS_EXEDIR=f'{binary.parent}/bsd:{binary.parent}',
                 GAUSS_SCRDIR=str(scratch),OMP_NUM_THREADS=str(self.args.threads))
        with prefix.with_suffix('.gjf').open('rb') as src,prefix.with_suffix('.log').open('wb') as dst:
            process=subprocess.run([str(binary)],cwd=self.directory,env=env,stdin=src,stdout=dst,
                                   stderr=subprocess.STDOUT,timeout=self.args.evaluation_timeout)
        log=prefix.with_suffix('.log').read_text(errors='replace')
        if process.returncode or 'Normal termination of Gaussian' not in log:
            raise RuntimeError(f'Gaussian failed; inspect {prefix.with_suffix(".log")}')
        formatted=prefix.with_suffix('.fchk')
        process=subprocess.run([str(binary.parent/'formchk'),str(checkpoint),str(formatted)],
                               cwd=self.directory,env=env,capture_output=True,text=True,timeout=120)
        if process.returncode:
            raise RuntimeError(f'Gaussian formchk failed: {process.stderr[-1000:]}')
        contents=formatted.read_text()
        from pyscf.data.nist import BOHR
        from pyscf.data.elements import charge
        coordinates=fchk_field(contents,'Current cartesian coordinates').reshape(-1,3)*BOHR
        numbers=fchk_field(contents,'Atomic numbers')
        if (coordinates.shape!=(len(atoms),3) or not np.allclose(coordinates,[r for e,r in atoms],rtol=0,atol=1e-7)
                or not np.array_equal(numbers,[charge(e) for e,r in atoms])):
            raise RuntimeError('Gaussian changed QM coordinate frame or atom ordering')
        result=gaussian_result(log,contents,len(atoms),gradient=gradient)
        result['qm_coordinate_frame_and_order_verified']=True
        return result

    def close(self):
        if self.scanner is not None:
            self.mf.mol.stdout.close()


def optimize(atoms,charge,environment,args,directory):
    from pyscf import gto
    from pyscf.geomopt import as_pyscf_method, geometric_solver
    directory=Path(directory)
    calc=Calculator(atoms,charge,environment,args,directory)
    mol=gto.M(atom=atoms,basis=args.basis,charge=charge,spin=int(charge!=0),unit='Angstrom',cart=False,verbose=0)
    trace=[]
    def evaluate(mol):
        current=[(mol.atom_symbol(i),tuple(r)) for i,r in enumerate(mol.atom_coords(unit='Angstrom'))]
        value=calc.evaluate(current)
        trace.append({'coordinates_angstrom':mol.atom_coords(unit='Angstrom').tolist(),**value})
        write_json(directory/'trace.json',trace)
        return value['energy_hartree'],np.array(value['gradient_hartree_per_bohr'])
    begin=time.perf_counter()
    try:
        converged,result=geometric_solver.kernel(as_pyscf_method(mol,evaluate),maxsteps=args.maxsteps,
                                                 assert_convergence=True)
        if not converged or not trace:
            raise RuntimeError(f'charge {charge} optimization not converged')
        optimized=[(result.atom_symbol(i),tuple(map(float,r))) for i,r in enumerate(result.atom_coords(unit='Angstrom'))]
        (directory/'optimized.xyz').write_text(xyz_text([e for e,r in optimized],[r for e,r in optimized],'optimized QM; fixed MM'))
        report={'converged':True,'evaluations':len(trace),'elapsed_seconds':time.perf_counter()-begin,
                'final_gradient_rms':trace[-1]['gradient_rms'],'final_gradient_max':trace[-1]['gradient_max']}
        write_json(directory/'optimization.json',report)
        return optimized,report
    finally:
        calc.close()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--engine',choices=['g16','gpu','cpu'],required=True)
    parser.add_argument('--xyz',type=Path,required=True)
    parser.add_argument('--mm-charges',type=Path,required=True)
    parser.add_argument('--lj',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--task',choices=['force','four-point'],default='force')
    parser.add_argument('--method',default='b3lyp')
    parser.add_argument('--basis',default='6-31G(d,p)')
    parser.add_argument('--threads',type=int,default=4)
    parser.add_argument('--memory-mb',type=int,default=24000)
    parser.add_argument('--conv-tol',type=float,default=1e-10)
    parser.add_argument('--maxsteps',type=int,default=100)
    parser.add_argument('--evaluation-timeout',type=int,default=7200)
    parser.add_argument('--g16',default='/opt/g16_C01_AVX2/g16/g16')
    args=parser.parse_args()
    if args.threads<1 or args.memory_mb<1 or args.conv_tol<=0:
        parser.error('threads, memory and SCF tolerance must be positive')
    if args.engine=='g16' and args.conv_tol!=1e-10:
        parser.error('G16 comparison currently requires conv-tol=1e-10 (SCF Conver=10)')
    args.output=args.output.resolve()
    args.output.mkdir(parents=True,exist_ok=False)
    atoms=read_xyz(args.xyz)
    environment=read_charges(args.mm_charges)
    environment=replace(environment,vdw=read_lj_file(args.lj,atoms,environment))
    snapshot(environment,args.output/'mm-environment.json')
    report={'schema':'cube-qmmm-engine-benchmark-1','status':'running','engine':args.engine,'task':args.task,
            'settings':{k:str(v) if isinstance(v,Path) else v for k,v in vars(args).items()},
            'input_sha256':{k:file_hash(v) for k,v in [('xyz',args.xyz),('mm_charges',args.mm_charges),('lj',args.lj)]},
            'qm_atoms':len(atoms),'mm_atoms':len(environment.charges_e),'mm_frozen':True,
            'embedding_fingerprint':environment.fingerprint,'results':{},
            'scf_algorithm':'Gaussian_conventional_integrals' if args.engine=='g16' else 'PySCF_density_fitting_def2_universal_jkfit',
            'grid':'nominal_99_radial_590_angular_native_pruning_and_quadrature_differ_between_engines',
            'optimizer':'same_geomeTRIC_QM_only','gaussian_scope':'external_charge_QM_plus_cube_LJ_not_native_ONIOM'}
    write_json(args.output/'result.json',report)
    begin=time.perf_counter()
    def single(geometry,charge,label,gradient=False):
        calculator=Calculator(geometry,charge,environment,args,args.output/label)
        try:return calculator.evaluate(geometry,gradient=gradient)
        finally:calculator.close()
    try:
        if args.task=='force':
            for charge in [0,1,-1]:
                report['results'][str(charge)]=single(atoms,charge,f'q{charge:+d}',gradient=True)
                write_json(args.output/'result.json',report)
        else:
            from pyscf.data.nist import HARTREE2EV
            r0,opt0=optimize(atoms,0,environment,args,args.output/'neutral')
            report['neutral_optimization']=opt0
            e0r0=single(r0,0,'E0_R0')
            for carrier,charge in [('hole',1),('electron',-1)]:
                rq,optq=optimize(r0,charge,environment,args,args.output/carrier)
                points={'E0_R0':e0r0,'Eq_R0':single(r0,charge,f'{carrier}_Eq_R0'),
                        'Eq_Rq':single(rq,charge,f'{carrier}_Eq_Rq'),'E0_Rq':single(rq,0,f'{carrier}_E0_Rq')}
                terms=four_point_terms(*(points[k]['energy_hartree'] for k in ['E0_R0','Eq_R0','Eq_Rq','E0_Rq']))
                if min(terms['charged_relaxation'],terms['neutral_relaxation'])*HARTREE2EV < -1e-5:
                    raise RuntimeError('negative relaxation energy; inspect optimized states')
                report['results'][carrier]={'optimization':optq,'single_points':points,
                    'relaxation_hartree':terms,'reorganization_energy_ev':terms['total']*HARTREE2EV}
                write_json(args.output/'result.json',report)
        report.update(status='completed',elapsed_seconds=time.perf_counter()-begin)
        write_json(args.output/'result.json',report)
    except Exception as exc:
        report.update(status='failed',elapsed_seconds=time.perf_counter()-begin,error=f'{type(exc).__name__}: {exc}')
        write_json(args.output/'result.json',report)
        raise


if __name__=='__main__':
    main()
