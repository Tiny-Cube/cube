#!/usr/bin/env python3
"""Validate real molecular couplings and exercise the universal KMC engines.

Requires numpy, optionally matplotlib for plots. Does not modify source data.
Run with the system Python or another environment providing these analysis deps.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import random
import re
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
import numpy as np
from universal_kmc import Edge, Engine, Graph, marcus
from universal_kmc.rates import KB, HBAR

MASSES = {'H': 1.008, 'C': 12.011, 'N': 14.007, 'O': 15.999, 'S': 32.06}


def check(ok, message):
    if not ok:
        raise ValueError(message)


def rows(path):
    with Path(path).open(newline='') as f:
        return list(csv.DictReader(f))


def write_csv(path, header, data):
    with Path(path).open('w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(data)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_geometry(source):
    lines = (source / 'anneal_a_fix.gro').read_text().splitlines()
    n_atoms = int(lines[1])
    check(len(lines) == n_atoms + 3, 'GRO length mismatch')
    box = np.array([float(v) for v in lines[-1].split()])
    check(box.shape == (3,) and np.all(np.isfinite(box)) and np.all(box > 0), 'orthorhombic GRO required')
    molecules = {}
    for line in lines[2:-1]:
        resid = int(line[:5])
        element = re.match(r'[A-Za-z]+', line[10:15].strip())[0]
        check(element in MASSES, f'unsupported element {element}')
        xyz = [float(line[start:start+8]) for start in (20, 28, 36)]
        molecules.setdefault(resid, []).append((element, xyz))
    ids = sorted(molecules)
    check(ids == list(range(1, len(ids)+1)), 'expected consecutive 1-based molecule IDs')
    core = np.array([int(v)-1 for v in (source/'core.txt').read_text().split()])
    check(len(core) > 0 and len(set(core)) == len(core) and np.min(core) >= 0, 'invalid 1-based core IDs')
    centers, positions, elements, weights = [], [], [], []
    for resid in ids:
        element = [a[0] for a in molecules[resid]]
        xyz = np.array([a[1] for a in molecules[resid]])
        check(core.max() < len(xyz) and np.isfinite(xyz).all(), 'invalid molecular geometry')
        pts = xyz[core]
        pts = pts[0] + (pts - pts[0] - box*np.rint((pts-pts[0])/box))
        w = np.array([MASSES[element[i]] for i in core])
        centers.append(np.average(pts, axis=0, weights=w))
        positions.append(xyz)
        elements.append(element)
        weights.append(w)
    return box, core, np.array(centers), positions, elements, weights, n_atoms


def gjf_coordinates(path):
    lines = path.read_text().splitlines()
    start = next(i for i, line in enumerate(lines) if line.strip() == '0 1') + 1
    atoms = []
    for line in lines[start:]:
        if not line.strip():
            break
        fields = line.split()
        check(len(fields) == 4, f'unsupported GJF coordinate row: {path}')
        atoms.append((fields[0], [float(v)/10 for v in fields[1:]]))
    return [a[0] for a in atoms], np.array([a[1] for a in atoms])


def prepare(source, output, temperature):
    box, core, centers, positions, elements, weights, n_atoms = load_geometry(source)
    table = rows(source/'transfer_integrals.csv')
    validated = {(int(r['mol_i']), int(r['mol_j'])): r for r in rows(source/'momap/transfer_integrals.csv')}
    energy_rows = [r for r in rows(source/'e_006.csv') if r['name'] == 'e_006_ether_12-17_0']
    check(len(energy_rows) == 1, 'missing/ambiguous reorganization-energy row')
    energies = energy_rows[0]
    lambdas = {'hole': float(energies['re_cat']), 'electron': float(energies['re_ani'])}
    check(all(math.isfinite(v) and v > 0 for v in lambdas.values()), 'invalid reorganization energy')
    # Verify that the supplied reorganization energies use the stated Hartree->eV conversion.
    expected_h = 27.211386245981*(abs(float(energies['neu_me'])-float(energies['cat'])) + abs(float(energies['cat_pe'])-float(energies['neu'])))
    expected_e = 27.211386245981*(abs(float(energies['ani_me'])-float(energies['neu'])) + abs(float(energies['neu_pe'])-float(energies['ani'])))
    check(abs(expected_h-lambdas['hole']) < 1e-8 and abs(expected_e-lambdas['electron']) < 1e-8, 'lambda unit/formula mismatch')
    checked = []
    seen = set()
    files = {name: sha(source/name) for name in [
        'transfer_integrals.csv', 'momap/transfer_integrals.csv', 'anneal_a_fix.gro', 'core.txt', 'e_006.csv']}
    coupling_count = 0
    for row in table:
        i, j = int(row['mol_i'])-1, int(row['mol_j'])-1
        check(0 <= i < j < len(centers), 'pair IDs out of bounds or not canonical')
        check((i, j) not in seen, 'duplicate dimer')
        seen.add((i, j))
        name = f'{i+1}_{j+1}'
        ref = validated[(i+1, j+1)]
        check(ref['2mol'] == name, 'validated pair label mismatch')
        for channel in ['VH', 'VL']:
            value = float(row[channel])
            check(math.isfinite(value) and value == float(ref[channel]), 'coupling table mismatch')
            path = source/f'momap/{name}_{channel}.inf'
            fields = path.read_text().split()
            check(len(fields) == 2 and fields[1] == 'meV' and float(fields[0]) == value, f'MOMAP result mismatch: {path}')
            files[path.relative_to(source).as_posix()] = sha(path)
            coupling_count += 1
        path = source/f'mol/frame2000/{name}.gjf'
        labels, xyz = gjf_coordinates(path)
        files[path.relative_to(source).as_posix()] = sha(path)
        n_i, n_j = len(positions[i]), len(positions[j])
        check(len(xyz) == n_i+n_j and labels == elements[i]+elements[j], 'GJF molecule/atom mapping mismatch')
        # Check every atom up to independent molecular translations, not just COM distance.
        shape_errors = []
        for mol, coords in [(i, xyz[:n_i]), (j, xyz[n_i:])]:
            delta = (coords - coords[0]) - (positions[mol] - positions[mol][0])
            delta -= box*np.rint(delta/box)
            shape_errors.append(float(np.linalg.norm(delta, axis=1).max()))
        check(max(shape_errors) < .002, f'{name}: GJF/GRO atom geometry mismatch')
        ref_vector = np.average(xyz[n_i+core],axis=0,weights=weights[j])-np.average(xyz[core],axis=0,weights=weights[i])
        raw = centers[j]-centers[i]
        image = np.rint((ref_vector-raw)/box).astype(int)
        vector = raw + image*box
        gjf_error = float(np.linalg.norm(vector-ref_vector))
        distance = float(np.linalg.norm(vector))
        check(gjf_error < .002, f'{name}: wrong periodic image')
        check(abs(distance-float(ref['com'])) < .002, f'{name}: validated COM mismatch')
        rates = {carrier: marcus(float(row['VH' if carrier=='hole' else 'VL'])/1000, lam, 0.0, temperature)
                 for carrier, lam in lambdas.items()}
        checked.append(dict(mol_i=i+1,mol_j=j+1,VH_meV=float(row['VH']),VL_meV=float(row['VL']),
                            input_com_nm=float(row['com']), corrected_com_nm=distance,
                            image_a=int(image[0]),image_b=int(image[1]),image_c=int(image[2]),
                            dx_nm=float(vector[0]),dy_nm=float(vector[1]),dz_nm=float(vector[2]),
                            gjf_vector_error_nm=gjf_error, max_atom_shape_error_nm=max(shape_errors),
                            hole_rate_s_1=rates['hole'],electron_rate_s_1=rates['electron']))
    check({(i+1,j+1) for i,j in seen} == validated.keys(), 'pair completeness mismatch')
    write_csv(output/'rates_and_geometry.csv', list(checked[0]), [list(r.values()) for r in checked])
    write_csv(output/'sites.csv', ['site','mol_id','core_x_nm','core_y_nm','core_z_nm'],
              [[i,i+1,*center] for i,center in enumerate(centers)])
    errors = [abs(r['input_com_nm']-r['corrected_com_nm']) for r in checked]
    audit = dict(atoms=n_atoms, molecules=len(centers), pairs=len(checked),box_nm=box.tolist(),
                 core_atom_ids_1based=(core+1).tolist(),couplings_verified_against_raw_momap=coupling_count,
                 gjf_pairs_verified=len(checked),root_label_values=sorted(set(r['2mol'] for r in table)),
                 input_com_errors_gt_0_002nm=sum(d>.002 for d in errors),max_input_com_error_nm=max(errors),
                 max_gjf_vector_error_nm=max(r['gjf_vector_error_nm'] for r in checked),
                 max_atom_shape_error_nm=max(r['max_atom_shape_error_nm'] for r in checked),
                 lambda_eV=lambdas,temperature_K=temperature,delta_G_eV=0.0,
                 assumptions=['static periodic morphology','zero field','no site-energy disorder',
                              'classical Marcus','one carrier','550 supplied dimer channels; no invented missing couplings'],
                 input_sha256=files)
    (output/'input_audit.json').write_text(json.dumps(audit,indent=2))
    return checked, audit


def make_edges(table, carrier):
    edges = []
    for row in table:
        i,j = row['mol_i']-1,row['mol_j']-1
        vector=tuple(row[k]*1e-9 for k in ['dx_nm','dy_nm','dz_nm'])
        rate=row[f'{carrier}_rate_s_1']
        edges.extend([Edge(i,j,rate,vector), Edge(j,i,rate,tuple(-x for x in vector))])
    return edges


def components(n, edges):
    adjacency=[set() for _ in range(n)]
    for edge in edges:
        if edge.rate>0: adjacency[edge.source].add(edge.target)
    unseen=set(range(n)); groups=[]
    while unseen:
        todo=[min(unseen)]; visited=set(todo)
        while todo:
            for j in adjacency[todo.pop()]-visited:
                visited.add(j); todo.append(j)
        unseen-=visited; groups.append(sorted(visited))
    return sorted(groups,key=len,reverse=True)


def save_model(path,n,edges):
    with path.open('w') as f:
        f.write('# Actual couplings; Marcus 300 K unless overridden. SI: rate s^-1, displacement m.\n')
        f.write(f'KMC1 graph\n{n} {len(edges)} 0\n')
        for e in edges:
            f.write(f'{e.source} {e.target} {e.rate:.17g} '+ ' '.join(f'{x:.17g}' for x in e.displacement)+'\n')


def cross_language(output,carrier,n,edges):
    model_path=output/f'{carrier}.kmc'
    save_model(model_path,n,edges)
    commands={'python':[sys.executable,'-m','universal_kmc'],
              'cpp':[str(ROOT/'build/ukmc-cpp')], 'rust':[str(ROOT/'rust/target/release/ukmc-rust')],
              'fortran':[str(ROOT/'build/ukmc-fortran')]}
    env=dict(os.environ,PYTHONPATH=str(ROOT/'src'))
    results=[]
    for kind,args in [('steps',['10000','42']),('time',['100000','20261102','1e-9'])]:
        reference=None
        for language,command in commands.items():
            path=output/f'{carrier}_{kind}_{language}.csv'
            began=time.perf_counter()
            with path.open('w') as f:
                proc=subprocess.run([*command,str(model_path),*args],stdout=f,stderr=subprocess.PIPE,text=True,
                                    cwd=ROOT,env=env,timeout=120)
            check(proc.returncode==0,f'{language}: {proc.stderr}')
            data=rows(path)
            if reference is None: reference=data
            check(len(data)==len(reference),f'{language}: trajectory length differs')
            max_time_error=0.
            for row,ref in zip(data,reference):
                for key in ['step','event','state0']:
                    check(int(row[key])==int(ref[key]),f'{language}: event/state differs')
                for key in ['time','x','y','z']:
                    actual,expected=float(row[key]),float(ref[key])
                    check(math.isclose(actual,expected,rel_tol=3e-13,abs_tol=1e-22),f'{language}: numeric mismatch')
                max_time_error=max(max_time_error,abs(float(row['time'])-float(ref['time'])))
            results.append(dict(mode=kind,language=language,events=int(data[-1]['step']),time_s=float(data[-1]['time']),
                                visited_sites=len({r['state0'] for r in data}),rows=len(data),
                                max_time_difference_vs_python_s=max_time_error,elapsed_seconds=time.perf_counter()-began,
                                status=proc.stderr.strip(),passed=True))
        print(carrier,kind,'all four implementations agree',flush=True)
    return results


class SparseHopping:
    """Application adapter over arbitrary outgoing events; uses the unchanged Engine."""
    def __init__(self,adjacency,site):
        self.adjacency=adjacency; self.site=site; self.xyz=(0.,0.,0.)
        self.rates=[[e.rate for e in a] for a in adjacency]
    def propensities(self): return self.rates[self.site]
    def fire(self,event):
        edge=self.adjacency[self.site][event]
        self.site=edge.target
        self.xyz=tuple(a+b for a,b in zip(self.xyz,edge.displacement))
    def state(self): return (self.site,)
    def position(self): return self.xyz


class PythonRNG:
    def __init__(self,seed): self.rng=random.Random(seed)
    def uniform(self):
        while True:
            u=self.rng.random()
            if u>0: return u


def exact_benchmark(n,edges,grid):
    """Symmetric periodic CTMC: corrector diffusivity and stationary finite-time MSD."""
    L=np.zeros((n,n)); b=np.zeros((n,3))
    for e in edges:
        L[e.source,e.source]+=e.rate; L[e.source,e.target]-=e.rate
        b[e.source]+=e.rate*np.array(e.displacement)*1e9
    scale=float(np.max(np.diag(L)))
    eigen,vectors=np.linalg.eigh(L/scale)
    check(eigen[1]>1e-14,'network must be connected for this benchmark')
    eigen=eigen[1:]; modes=vectors[:,1:].T@(b/scale)
    chi_modes=modes/eigen[:,None]
    chi=vectors[:,1:]@chi_modes
    residual=float(np.max(np.abs((L/scale)@chi-b/scale)))
    check(residual<1e-8,'corrector residual too large')
    tensor=np.zeros((3,3))
    for e in edges:
        d=np.array(e.displacement)*1e9+chi[e.target]-chi[e.source]
        tensor+=e.rate*np.outer(d,d)/(2*n)
    transient=2/n*((-np.expm1(-np.outer(grid,eigen*scale)))@np.sum(chi_modes**2,axis=1))
    msd=2*np.trace(tensor)*grid+transient
    return tensor,msd,residual


def ensemble(output,carrier,n,edges,temperature,trajectories,duration):
    adjacency=[[] for _ in range(n)]
    for e in edges: adjacency[e.source].append(e)
    # Check the sparse adapter has the same actual trajectory as the generic full graph.
    sparse=SparseHopping(adjacency,0); dense=Graph(n,edges)
    a,b=Engine(sparse,42),Engine(dense,42)
    for _ in range(2000):
        ra,rb=a.step(),b.step()
        check(ra.state==rb.state and ra.position==rb.position and ra.time==rb.time,'sparse adapter changes trajectory')
    grid=np.linspace(0,duration,101)
    rng=PythonRNG(20261102+(0 if carrier=='hole' else 1))
    squared=np.zeros((trajectories,len(grid)))
    endpoints=np.zeros((trajectories,3))
    hops=[]
    tensor,exact_msd,residual=exact_benchmark(n,edges,grid)
    for trajectory in range(trajectories):
        initial=rng.rng.randrange(n)
        model=SparseHopping(adjacency,initial)
        engine=Engine(model,rng=rng)
        pointer=1
        trace=[]
        if trajectory<3: trace.append([0,0,-1,initial,0.,0.,0.])
        while engine.time<duration:
            old_position=np.array(model.xyz)*1e9
            record=engine.step(duration)
            # Waiting interval has the PRE-event position; no linear interpolation.
            while pointer<len(grid) and grid[pointer]<engine.time:
                squared[trajectory,pointer]=float(old_position@old_position)
                pointer+=1
            if trajectory<3:
                trace.append([engine.steps,engine.time,-1 if record is None else record.event,
                              model.site,*[x*1e9 for x in model.xyz]])
            if record is None: break
            check(engine.steps<1000000,'ensemble hop cap reached')
        final=np.array(model.xyz)*1e9
        squared[trajectory,pointer:]=float(final@final)
        endpoints[trajectory]=final
        hops.append(engine.steps)
        if trajectory<3:
            write_csv(output/f'{carrier}_trace_{trajectory+1}.csv',
                      ['step','time_s','local_event','site','dx_nm','dy_nm','dz_nm'],trace)
        if (trajectory+1)%32==0:
            print(f'{carrier}: {trajectory+1}/{trajectories} fixed-time trajectories',flush=True)
    mean=squared.mean(axis=0); sem=squared.std(axis=0,ddof=1)/math.sqrt(trajectories)
    z=float((mean[-1]-exact_msd[-1])/sem[-1])
    check(abs(z)<4,'ensemble MSD differs from exact CTMC by >=4 standard errors')
    # Late-window per-trajectory slopes give an uncertainty estimate without treating times as independent.
    start=len(grid)//2
    x=grid[start:]-grid[start:].mean()
    per_d=(squared[:,start:]@x/(x@x))/6
    mu=per_d*1e-14/(KB*temperature)
    mu_mean=float(mu.mean()); mu_sem=float(mu.std(ddof=1)/math.sqrt(trajectories))
    write_csv(output/f'{carrier}_msd.csv',['time_s','mean_r2_nm2','sem_r2_nm2','exact_finite_time_r2_nm2'],
              zip(grid,mean,sem,exact_msd))
    exact_mu=float(np.trace(tensor)/3*1e-14/(KB*temperature))
    result=dict(trajectories=trajectories,duration_s=duration,total_events=sum(hops),
                hops_min_median_max=[min(hops),float(np.median(hops)),max(hops)],
                final_msd_nm2=float(mean[-1]),final_msd_sem_nm2=float(sem[-1]),
                exact_final_msd_nm2=float(exact_msd[-1]),final_msd_z_score=z,
                mean_displacement_nm=endpoints.mean(axis=0).tolist(),
                exact_diffusion_tensor_nm2_s=tensor.tolist(),exact_model_mobility_cm2_Vs=exact_mu,
                estimated_model_mobility_cm2_Vs=mu_mean,mobility_95ci_cm2_Vs=[mu_mean-1.96*mu_sem,mu_mean+1.96*mu_sem],
                corrector_residual_nm=residual,sparse_adapter_2000_events_verified=True,
                sampler='Python standard random.Random; shared stream; universal_kmc.Engine',passed=True)
    return result,grid,mean,sem,exact_msd


def plot(output,curves):
    os.environ.setdefault('MPLCONFIGDIR','/tmp/ukmc-mpl-cache')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(11,4),constrained_layout=True)
    for ax,(carrier,grid,mean,sem,exact) in zip(axes,curves):
        ax.plot(grid*1e9,mean,label='Universal KMC ensemble',lw=2)
        ax.fill_between(grid*1e9,mean-1.96*sem,mean+1.96*sem,alpha=.18,label='95% MC interval')
        ax.plot(grid*1e9,exact,'--',color='black',label='Exact finite-time CTMC')
        ax.set(xlabel='Time (ns)',ylabel='Mean squared displacement (nm²)',title=carrier.capitalize())
        ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.suptitle('e_006_ether_12-17_0 · 300 K · zero field · static periodic network')
    fig.savefig(output/'kmc_validation.png',dpi=170)
    fig.savefig(output/'kmc_validation.svg')
    plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--temperature',type=float,default=300.)
    parser.add_argument('--trajectories',type=int,default=256)
    parser.add_argument('--time-ns',type=float,default=10.)
    args=parser.parse_args()
    check(args.trajectories>=2 and args.time_ns>0 and math.isfinite(args.time_ns),'invalid ensemble parameters')
    output=args.output.resolve();output.mkdir(parents=True,exist_ok=True)
    check(not (output/'summary.json').exists(),'output already completed; choose a new directory')
    began=time.perf_counter()
    table,audit=prepare(args.source.resolve(),output,args.temperature)
    summary={'system':'e_006_ether_12-17_0','input_audit':'input_audit.json','carriers':{}}
    curves=[]
    for carrier in ['hole','electron']:
        edges=make_edges(table,carrier)
        groups=components(audit['molecules'],edges)
        check(len(groups)==1,'disconnected coupling graph')
        rates=[e.rate for e in edges]
        print(carrier,'positive directed edges',sum(r>0 for r in rates),'rate range',min(rates),max(rates),flush=True)
        parity=cross_language(output,carrier,audit['molecules'],edges)
        result,grid,mean,sem,exact=ensemble(output,carrier,audit['molecules'],edges,args.temperature,args.trajectories,args.time_ns*1e-9)
        summary['carriers'][carrier]={'lambda_eV':audit['lambda_eV'][carrier],
                                     'directed_edges':len(edges),'positive_directed_edges':sum(r>0 for r in rates),
                                     'rate_min_max_s_1':[min(rates),max(rates)],
                                     'components':[len(c) for c in groups],'cross_language':parity,'ensemble':result}
        curves.append((carrier,grid,mean,sem,exact))
        (output/'summary.partial.json').write_text(json.dumps(summary,indent=2))
    plot(output,curves)
    summary['elapsed_seconds']=time.perf_counter()-began
    summary['passed']=True
    (output/'summary.json').write_text(json.dumps(summary,indent=2))
    (output/'summary.partial.json').unlink(missing_ok=True)
    print(json.dumps({k:v['ensemble'] for k,v in summary['carriers'].items()},indent=2),flush=True)


if __name__=='__main__':
    main()
