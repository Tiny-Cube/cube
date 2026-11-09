#!/usr/bin/env python3
"""Independent split-operator quantum reference vs Wigner-sampled Tully-1 FSSH.

Not a unit test: the stochastic estimator has sampling error and FSSH is an
approximation to quantum dynamics. No external published numbers are assumed.
"""
import argparse
import json
from pathlib import Path
import time
import numpy as np

from universal_kmc.quantum.cli import _write_result
from universal_kmc.surface_hopping.engine import FSSH
from universal_kmc.surface_hopping.models import TullyModel


def quantum_scattering(points=2048, dt=1., duration=2200.):
    # Minimum uncertainty wavepacket: sigma_x=1, sigma_p=0.5 (atomic units).
    x = np.linspace(-35.,35.,points,endpoint=False)
    dx = x[1]-x[0]
    # Construct original diabatic potentials directly, independently of provider.
    v = np.sign(x)*.01*(1-np.exp(-1.6*abs(x)))
    coupling = .005*np.exp(-x*x)
    h = np.zeros((points,2,2))
    h[:,0,0], h[:,1,1] = v,-v
    h[:,0,1] = h[:,1,0] = coupling
    energies,vectors = np.linalg.eigh(h)
    half = np.einsum('xai,xi,xbi->xab',vectors,np.exp(-.5j*dt*energies),vectors)
    momentum = 2*np.pi*np.fft.fftfreq(points,d=dx)
    kinetic = np.exp(-1j*momentum*momentum/(2*2000.)*dt)
    psi = np.zeros((points,2),dtype=complex)
    psi[:,0] = np.exp(-(x+8.)**2/4+20j*(x+8.))
    psi /= np.sqrt(np.sum(abs(psi)**2)*dx)
    for _ in range(round(duration/dt)):
        psi = np.einsum('xab,xb->xa',half,psi)
        psi = np.fft.ifft(kinetic[:,None]*np.fft.fft(psi,axis=0),axis=0)
        psi = np.einsum('xab,xb->xa',half,psi)
    adiabatic = np.einsum('xai,xa->xi',vectors,psi)
    return {'right':(np.sum(abs(adiabatic[x>0])**2,axis=0)*dx).tolist(),
            'left':(np.sum(abs(adiabatic[x<0])**2,axis=0)*dx).tolist(),
            'central_probability':float(np.sum(abs(psi[abs(x)<3])**2)*dx),
            'norm':float(np.sum(abs(psi)**2)*dx), 'points':points,'dt_au':dt,'duration_au':duration}


def fssh_scattering(trajectories,seed,dt):
    engine=FSSH(TullyModel(1),[2000.],dt=dt,electronic_substeps=4)
    children=np.random.SeedSequence(seed).spawn(trajectories)
    counts={'left':[0,0],'right':[0,0],'unfinished':[0,0]}
    max_drift=max_jump=max_norm=0.
    for index,child in enumerate(children):
        rng=np.random.Generator(np.random.PCG64(child))
        state=engine.initialize([rng.normal(-8.,1.)],[rng.normal(20.,.5)])
        e0=engine.energy(state)
        reason='unfinished'
        for _ in range(round(2200/dt)):
            state,event=engine.step(state,rng)
            max_drift=max(max_drift,abs(engine.energy(state)-e0))
            max_jump=max(max_jump,abs(event['jump_energy_error']))
            max_norm=max(max_norm,abs(float(np.vdot(state.c,state.c).real)-1))
            if state.q[0]>8.:
                reason='right';break
            if state.q[0]<-12.:
                reason='left';break
        counts[reason][state.active]+=1
        if (index+1)%max(1,trajectories//8)==0:
            print(f'reference validation: {index+1}/{trajectories}',flush=True)
    return {'trajectories':trajectories,'seed':seed,'dt_au':dt,'counts':counts,
            'fractions':{k:(np.array(v)/trajectories).tolist() for k,v in counts.items()},
            'max_energy_drift_hartree':max_drift,'max_hop_energy_error_hartree':max_jump,'max_norm_error':max_norm}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--trajectories',type=int,default=128)
    parser.add_argument('--seed',type=int,default=29101)
    parser.add_argument('--dt',type=float,default=1.)
    parser.add_argument('--output',type=Path,default=Path('build/surface-hopping/reference.json'))
    args=parser.parse_args()
    if args.trajectories<1 or args.dt<=0:parser.error('positive trajectories and dt required')
    start=time.perf_counter()
    coarse=quantum_scattering()
    refined=quantum_scattering(points=4096,dt=.5)
    fssh=fssh_scattering(args.trajectories,args.seed,args.dt)
    p=np.array(refined['right']);observed=np.array(fssh['fractions']['right'])
    se=np.sqrt(p*(1-p)/args.trajectories)
    report={'schema':'cube-fssh-validation-1','model':'Tully 1','scope':'single high-momentum scattering benchmark',
            'reference_doi':'10.1063/1.459170','initial':'Gaussian pure diabatic 0; q=-8 sigma_q=1, p=20 sigma_p=0.5; FSSH uses matching Wigner distribution',
            'quantum_coarse':coarse,'quantum_refined':refined,'fssh':fssh,
            'quantum_refinement_max_difference':float(np.max(abs(np.array(coarse['right'])-p))),
            'right_fraction_difference':(observed-p).tolist(),'binomial_standard_errors_using_quantum_p':se.tolist(),
            'right_difference_in_standard_errors':((observed-p)/se).tolist(),
            'elapsed_seconds':time.perf_counter()-start,
            'interpretation':'Stochastic and approximate FSSH; one benchmark does not establish accuracy for material transport or all Tully regimes.'}
    _write_result(args.output,report)
    print(json.dumps({k:report[k] for k in ('right_fraction_difference','right_difference_in_standard_errors','quantum_refinement_max_difference','elapsed_seconds')},indent=2))


if __name__=='__main__':main()
