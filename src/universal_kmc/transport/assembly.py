"""Connect a single extraction frame + DIPRO CSV + completed four-point results.

No orbital energy is silently interpreted as a free energy. Pair lambda is
accepted from computed monomers only under an explicit equivalent-site model.
"""
import csv
import json
from pathlib import Path
import numpy as np

HARTREE_EV=27.211386245981


def assemble(manifest_path,couplings_path,reorganization,settings,sources):
    from .runner import read_source
    manifest_path=Path(manifest_path)
    manifest=read_source(manifest_path,sources)
    if manifest.get('schema')!='cube-dimers-1' or manifest.get('coordinate_unit')!='Angstrom':
        raise ValueError('KMC assembly requires one cube-dimers-1 Angstrom manifest per calculation')
    ids=[str(m['id']) for m in manifest['monomers']]
    if len(ids)!=len(set(ids)):raise ValueError('duplicate monomer IDs in KMC manifest')
    positions={}
    for m in manifest['monomers']:
        if 'transport_position_angstrom' in m:positions[str(m['id'])]=m['transport_position_angstrom']
    if len(positions)!=len(ids):
        # Old extraction records store centered isolated monomers; their centers
        # are not material positions. Reconstruct from verified MD provenance.
        from universal_kmc.workflow.embedding import load_gromacs_environment
        frame=load_gromacs_environment(manifest)
        if len(frame['fragments'])!=len(ids):raise ValueError('MD/frame molecule count changed')
        positions={identity:np.asarray(fragment['coordinates']).mean(axis=0).tolist()
                   for identity,fragment in zip(ids,frame['fragments'])}
        sources.update({manifest['provenance']['settings'][key]:value
                        for key,value in manifest['provenance']['sha256'].items()})
    if isinstance(reorganization,(str,Path)):
        index_path=Path(reorganization)
        reorganization=read_source(index_path,sources)
        reorganization={str(k):str((index_path.parent/v).resolve()) for k,v in reorganization.items()}
    records={}
    for identity in ids:
        if identity not in reorganization:raise ValueError(f'missing reorganization result for site {identity}')
        path=reorganization[identity]
        if isinstance(path,dict):path=path['output']
        result=read_source(path,sources)
        if result.get('schema')!='cube-reorganization-1' or result.get('status')!='completed':
            raise ValueError('KMC requires completed cube four-point results')
        records[identity]=result
    sites=[]
    for identity in ids:
        energy={}
        for carrier in settings['carriers']:
            row=records[identity]['results'][carrier]
            energy[carrier]=(0. if settings['site-energy-model']=='uniform' else
                (row['energies_hartree']['Eq_Rq']-row['energies_hartree']['E0_R0'])*HARTREE_EV)
        sites.append({'id':identity,'position_angstrom':positions[identity],'free_energy_ev':energy})
    read_source(couplings_path,sources,json_data=False)
    with Path(couplings_path).open(newline='') as f:table=list(csv.DictReader(f))
    by_pair={row['dimer']:row for row in table}
    if len(by_pair)!=len(table):raise ValueError('duplicate coupling CSV dimer IDs')
    pairs=[]
    for pair in manifest['dimers']:
        identity=str(pair['id']);a,b=str(pair['mol_i']),str(pair['mol_j'])
        if identity not in by_pair:raise ValueError(f'missing coupling for dimer {identity}')
        row=by_pair[identity]
        if str(row['mol_i'])!=a or str(row['mol_j'])!=b:raise ValueError('coupling CSV and manifest molecule IDs disagree')
        coupling={};lambdas={}
        for carrier,column in (('hole','VH_meV'),('electron','VL_meV')):
            if carrier not in settings['carriers']:continue
            if not row.get(column):raise ValueError(f'missing {column} for {identity}; missing is not zero')
            coupling[carrier]=float(row[column])/1000.
            ra=records[a]['results'][carrier]['relaxation_ev'];rb=records[b]['results'][carrier]['relaxation_ev']
            if min(ra['charged_relaxation'],ra['neutral_relaxation'],rb['charged_relaxation'],rb['neutral_relaxation']) < -1e-5:
                raise ValueError('negative computed relaxation energy cannot define pair lambda')
            forward=ra['charged_relaxation']+rb['neutral_relaxation']
            backward=rb['charged_relaxation']+ra['neutral_relaxation']
            if not np.isclose(forward,backward,rtol=1e-4,atol=1e-6):
                raise ValueError(f'inequivalent forward/reverse four-point lambdas for {identity}; provide physically defined pair lambdas in an explicit network')
            lambdas[carrier]=float(forward)
        pairs.append({'id':identity,'source':a,'target':b,'image':pair['image'],
                      'coupling_ev':coupling,'reorganization_ev':lambdas})
    if set(by_pair)!=set(str(p['id']) for p in manifest['dimers']):raise ValueError('coupling CSV contains pairs outside the extraction manifest')
    return {'schema':'cube-charge-network-1','lattice_angstrom':manifest['lattice_angstrom'],
            'sites':sites,'pairs':pairs,'provenance':{'scope':'static_snapshot_single_carrier',
            'frame':manifest.get('frame'),'position_definition':'whole_molecule_geometric_center_in_material',
            'site_energy_model':settings['site-energy-model'],
            'pair_reorganization_model':'charged_relaxation_donor_plus_neutral_relaxation_acceptor_equivalent_sites',
            'outer_sphere_relaxation':False,'source_sha256':dict(sources)}}
