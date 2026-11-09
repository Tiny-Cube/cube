"""Periodic molecular graphs: explicit eV/Angstrom data, no invented couplings."""
import math
import re
import numpy as np

from universal_kmc.models import Edge
from universal_kmc.rates import KB, HBAR, marcus


def _identity(value):
    text = str(value)
    if not re.fullmatch(r'[A-Za-z0-9_-]+', text):
        raise ValueError('site/pair IDs must be nonempty safe names')
    return text


def _real(value, name, positive=False):
    if type(value) not in (int,float) or not math.isfinite(value) or (positive and value <= 0):
        raise ValueError(f'{name} must be finite'+(' and positive' if positive else ''))
    return float(value)


def validate_network(data, carriers):
    if not isinstance(data,dict) or data.get('schema') != 'cube-charge-network-1':
        raise ValueError('KMC input must have schema cube-charge-network-1')
    if set(data)-{'schema','lattice_angstrom','sites','pairs','provenance'}:
        raise ValueError('unknown charge network fields')
    cell=np.asarray(data['lattice_angstrom'],dtype=float)
    if cell.shape!=(3,3) or not np.isfinite(cell).all() or np.linalg.det(cell)<=1e-10:
        raise ValueError('a finite right-handed 3x3 lattice in Angstrom is required')
    sites, pairs=data['sites'],data['pairs']
    if not isinstance(sites,list) or not sites or not isinstance(pairs,list):
        raise ValueError('network requires nonempty sites and a list of pairs')
    ids=[];positions=[];energies={c:[] for c in carriers}
    for site in sites:
        if set(site)!={'id','position_angstrom','free_energy_ev'}:
            raise ValueError('sites require id, position_angstrom, free_energy_ev')
        ids.append(_identity(site['id']))
        position=np.asarray(site['position_angstrom'],dtype=float)
        if position.shape!=(3,) or not np.isfinite(position).all():
            raise ValueError('site position must have 3 finite Angstrom coordinates')
        positions.append(position)
        for carrier in carriers:
            if carrier not in site['free_energy_ev']:
                raise ValueError(f'missing {carrier} site free energy; uniform energies must be explicit')
            energies[carrier].append(_real(site['free_energy_ev'][carrier],'site free energy'))
    if len(set(ids))!=len(ids):raise ValueError('duplicate site IDs')
    index={name:i for i,name in enumerate(ids)}
    seen=set();pair_ids=set();records=[]
    for pair in pairs:
        if set(pair)!={'id','source','target','image','coupling_ev','reorganization_ev'}:
            raise ValueError('pairs require id, source, target, image, coupling_ev, reorganization_ev')
        identity=_identity(pair['id'])
        if identity in pair_ids:raise ValueError('duplicate pair ID')
        pair_ids.add(identity)
        try:i,j=index[str(pair['source'])],index[str(pair['target'])]
        except KeyError as exc:raise ValueError('pair references an absent site') from exc
        image=pair['image']
        if not isinstance(image,list) or len(image)!=3 or any(type(v) is not int for v in image):
            raise ValueError('image must have 3 integer lattice translations')
        if i==j and not any(image):raise ValueError('a self hop requires a nonzero periodic image')
        forward=(i,j,*image);reverse=(j,i,*(-v for v in image))
        canonical=min(forward,reverse)
        if canonical in seen:raise ValueError('duplicate/reversed periodic pair channel')
        seen.add(canonical)
        displacement=(positions[j]-positions[i]+np.asarray(image)@cell)*1e-10
        if not np.isfinite(displacement).all():raise ValueError('non-finite hop displacement')
        for carrier in carriers:
            if carrier not in pair['coupling_ev'] or carrier not in pair['reorganization_ev']:
                raise ValueError(f'missing {carrier} coupling/reorganization for pair {identity}; missing is not zero')
            _real(pair['coupling_ev'][carrier],'coupling')
            _real(pair['reorganization_ev'][carrier],'pair reorganization',positive=True)
        records.append((pair,i,j,displacement))
    return ids,np.asarray(positions),energies,records


def build_rates(data,carrier,temperature):
    ids,positions,energies,pairs=validate_network(data,[carrier])
    if not math.isfinite(temperature) or temperature<=0:raise ValueError('positive temperature required')
    energy=np.asarray(energies[carrier])
    kt=KB*temperature
    weights=np.exp(-(energy-energy.min())/kt)
    equilibrium=weights/weights.sum()
    if np.any(equilibrium==0):raise ValueError('Boltzmann probabilities underflow; this energy range is not supported in double precision')
    edges=[];rows=[];images=[]
    for pair,i,j,d in pairs:
        coupling=float(pair['coupling_ev'][carrier]);lam=float(pair['reorganization_ev'][carrier])
        delta=float(energy[j]-energy[i])
        logs=[]
        for sign,source,target in ((1,i,j),(-1,j,i)):
            try:
                rate=marcus(coupling,lam,sign*delta,temperature)
                log_rate=(-math.inf if coupling==0 else math.log(2*math.pi/HBAR)+2*math.log(abs(coupling))
                          -.5*math.log(4*math.pi*lam*kt)-(sign*delta+lam)**2/(4*lam*kt))
            except OverflowError as exc:
                raise ValueError('Marcus parameters/rates exceed floating-point range') from exc
            logs.append(log_rate)
            edges.append(Edge(source,target,rate,tuple(float(v) for v in sign*d)))
            image=[sign*v for v in pair['image']];images.append(image)
            rows.append({'pair':str(pair['id']),'source':ids[source],'target':ids[target],
                         'image':image,'displacement_m':list(edges[-1].displacement),
                         'coupling_ev':coupling,'reorganization_ev':lam,'delta_g_ev':sign*delta,
                         'rate_s_1':rate,'log_rate_s_1':None if not math.isfinite(log_rate) else log_rate,
                         'zero_coupling':coupling==0,'underflow':coupling!=0 and rate==0})
        if coupling and abs((logs[0]-logs[1])+delta/kt)>1e-8*max(1.,abs(delta/kt)):
            raise RuntimeError('Marcus detailed-balance check failed')
    adjacency=[[] for _ in ids]
    for n,edge in enumerate(edges):
        if edge.rate>0:adjacency[edge.source].append((n,edge))
    labels={0:np.zeros(3,dtype=np.int64)};todo=[0];cycles=[]
    while todo:
        i=todo.pop()
        for e_index,e in adjacency[i]:
            candidate=labels[i]+images[e_index]
            if e.target not in labels:
                labels[e.target]=candidate;todo.append(e.target)
            else:cycles.append(candidate-labels[e.target])
    reverse_adjacency=[[] for _ in ids]
    for e in edges:
        if e.rate>0:reverse_adjacency[e.target].append(e.source)
    backward={0};todo=[0]
    while todo:
        for target in reverse_adjacency[todo.pop()]:
            if target not in backward:backward.add(target);todo.append(target)
    if len(labels)!=len(ids) or len(backward)!=len(ids):
        raise ValueError('positive-rate network is disconnected; choose a connected component explicitly')
    rank=int(np.linalg.matrix_rank(np.asarray(cycles,dtype=float))) if cycles and np.any(cycles) else 0
    return ids,edges,rows,equilibrium,rank


def exact_diffusion(sites,edges,pi,max_sites=512):
    """Zero-field reversible CTMC corrector, arbitrary site energies/parallel edges.

    Q chi = v-b, pi.chi=0; D = 1/2 sum pi_i k_ij (d+chi_j-chi_i)^2.
    Self-image hops cancel from Q while remaining in drift and diffusivity.
    """
    if sites>max_sites:return {'status':'skipped','reason':'dense_reference_site_limit','max_sites':max_sites}
    q=np.zeros((sites,sites));b=np.zeros((sites,3))
    for e in edges:
        if e.source!=e.target:
            q[e.source,e.target]+=e.rate;q[e.source,e.source]-=e.rate
        b[e.source]+=e.rate*np.asarray(e.displacement)
    scale=max(float(np.max(abs(q))),max((e.rate for e in edges),default=0.),1e-300)
    stationary=float(np.max(abs(pi @ q))/scale)
    if stationary>1e-9:raise ValueError('provided equilibrium does not solve the network generator')
    drift=pi @ b
    a=-q/scale+np.ones((sites,1))*pi
    condition=float(np.linalg.cond(a))
    if not math.isfinite(condition) or condition>1e14:
        return {'status':'skipped','reason':'ill_conditioned_dense_reference','condition_number':condition}
    chi=np.linalg.solve(a,(b-drift)/scale)
    residual=float(np.max(abs(a@chi-(b-drift)/scale)))
    relative=residual/max(float(np.max(abs((b-drift)/scale))),1e-300)
    if relative>1e-8:raise RuntimeError('corrector solve did not reach the required residual')
    tensor=np.zeros((3,3))
    for e in edges:
        d=np.asarray(e.displacement)+chi[e.target]-chi[e.source]
        tensor+=.5*pi[e.source]*e.rate*np.outer(d,d)
    return {'status':'completed','diffusion_tensor_m2_s':tensor.tolist(),
            'stationary_residual_scaled':stationary,'corrector_residual_m':residual,
            'corrector_relative_residual':relative,'condition_number':condition,
            'equilibrium':pi.tolist(),'drift_velocity_m_s':drift.tolist()}
