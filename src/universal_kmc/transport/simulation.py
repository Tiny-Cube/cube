"""Independent single-carrier trajectories observed at fixed physical times."""
import hashlib
import math
import random
import numpy as np
from universal_kmc import Engine
from universal_kmc.rates import KB


class OutgoingGraph:
    def __init__(self,adjacency,initial):
        self.adjacency=adjacency;self.site=initial;self.xyz=(0.,0.,0.)
    def propensities(self):return [edge.rate for _,edge in self.adjacency[self.site]]
    def fire(self,event):
        _,edge=self.adjacency[self.site][event]
        self.site=edge.target
        self.xyz=tuple(a+b for a,b in zip(self.xyz,edge.displacement))
    def state(self):return (self.site,)
    def position(self):return self.xyz


class IndependentRNG:
    """MT19937 with separately hashed seeds; uniform values strictly in (0,1)."""
    def __init__(self,seed):self.generator=random.Random(seed)
    def uniform(self):
        while True:
            value=self.generator.random()
            if value>0:return value


def trajectory_seed(seed,carrier,index,kind):
    return int.from_bytes(hashlib.sha256(f'{seed}:{carrier}:{index}:{kind}'.encode()).digest()[:8],'little')


def simulate(sites,edges,pi,settings,carrier):
    adjacency=[[] for _ in range(sites)]
    for index,edge in enumerate(edges):
        if edge.rate>0:adjacency[edge.source].append((index,edge))
    grid=np.linspace(0,settings['duration'],settings['samples'])
    positions=np.zeros((settings['trajectories'],len(grid),3))
    traces=[];counts=[];capped=[]
    for trajectory in range(settings['trajectories']):
        seed=trajectory_seed(settings['seed'],carrier,trajectory,'events')
        draw=IndependentRNG(trajectory_seed(settings['seed'],carrier,trajectory,'initial')).uniform()
        initial=min(int(np.searchsorted(np.cumsum(pi),draw,side='right')),sites-1)
        model=OutgoingGraph(adjacency,initial);engine=Engine(model,rng=IndependentRNG(seed))
        pointer=1
        trace=[{'step':0,'time_s':0.,'site':initial,'event':None,'displacement_m':[0.,0.,0.]}]
        while engine.time<settings['duration']:
            if engine.steps>=settings['max-events']:
                capped.append(trajectory);break
            before=model.position();old_site=model.site
            outgoing=adjacency[old_site]
            record=engine.step(settings['duration'])
            # A waiting interval holds the PRE-event coordinate. No interpolation.
            while pointer<len(grid) and grid[pointer]<engine.time:
                positions[trajectory,pointer]=before;pointer+=1
            if trajectory<settings['trace-trajectories']:
                trace.append({'step':engine.steps,'time_s':engine.time,'site':model.site,
                              'event':None if record is None else outgoing[record.event][0],
                              'displacement_m':list(model.position())})
            if record is None:break
        if trajectory in capped:
            positions[trajectory,pointer:]=np.nan
        else:positions[trajectory,pointer:]=model.position()
        counts.append(engine.steps)
        if trajectory<settings['trace-trajectories']:traces.append({'trajectory':trajectory,'initial_site':initial,'event_seed':seed,'records':trace})
        if (trajectory+1)%max(1,settings['trajectories']//8)==0:
            print(f'cube: {carrier} KMC {trajectory+1}/{settings["trajectories"]}',flush=True)
    return grid,positions,traces,counts,capped


def _fit(grid,covariance,start):
    x=grid[start:]-grid[start:].mean()
    return .5*np.einsum('t,tij->ij',x,covariance[start:])/float(x@x)


def analyze(grid,positions,temperature,fit_start,relative_target,rank):
    n=len(positions)
    if n<4:return {'status':'partial','reason':'at_least_four_trajectories_required_for_uncertainty'}
    mean=positions.mean(axis=0)
    outer=np.einsum('nti,ntj->tij',positions,positions)
    sums=positions.sum(axis=0)
    covariance=(outer-n*np.einsum('ti,tj->tij',mean,mean))/(n-1)
    start=min(int((len(grid)-1)*fit_start),len(grid)-3)
    tensor=_fit(grid,covariance,start)
    jack=[]
    for sample in positions:
        leave_sum=sums-sample
        leave_cov=(outer-np.einsum('ti,tj->tij',sample,sample)
                   -np.einsum('ti,tj->tij',leave_sum,leave_sum)/(n-1))/(n-2)
        jack.append(_fit(grid,leave_cov,start))
    jack=np.asarray(jack)
    sem=np.sqrt((n-1)/n*np.sum((jack-jack.mean(axis=0))**2,axis=0))
    later=_fit(grid,covariance,min(int((len(grid)-1)*(.5+.5*fit_start)),len(grid)-3))
    norm=max(float(np.linalg.norm(tensor)),1e-300)
    relative_uncertainty=float(1.96*np.linalg.norm(sem)/norm)
    window_change=float(np.linalg.norm(later-tensor)/norm)
    converged=relative_uncertainty<=relative_target and window_change<=relative_target
    if np.linalg.eigvalsh(tensor).min() < -3*float(np.linalg.norm(sem)):converged=False
    reason='within_requested_diagnostics' if converged else 'increase_duration_or_trajectories'
    if rank==0:converged=True;reason='bounded_network_no_unbounded_transport'
    conversion=1e4/(KB*temperature)
    return {'status':'completed' if converged else 'partial','reason':reason,
            'diffusion_tensor_m2_s':tensor.tolist(),'diffusion_tensor_sem_m2_s':sem.tolist(),
            'einstein_mobility_tensor_cm2_Vs':(tensor*conversion).tolist(),
            'mobility_tensor_sem_cm2_Vs':(sem*conversion).tolist(),
            'uncertainty_method':'trajectory jackknife; componentwise 1.96*SEM, norm diagnostic is not a simultaneous confidence region',
            'mean_displacement_m':mean.tolist(),'covariance_m2':covariance.tolist(),
            'fit_start_fraction':fit_start,'fit_start_index':start,
            'relative_95ci_norm':relative_uncertainty,'late_window_relative_change':window_change,
            'relative_target':relative_target,'periodic_transport_rank':rank,
            'note':'finite-window zero-field estimator; diagnostics do not prove material convergence'}
