"""Analytical periodic graphs, parameter integrity and charge workflow wiring."""
from contextlib import redirect_stdout
import copy
import importlib.util
import io
import json
import math
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from universal_kmc.workflow.control import load_input
from universal_kmc.transport.config import normalize

HAS_NUMPY=importlib.util.find_spec('numpy') is not None
ROOT=Path(__file__).resolve().parents[1]


class CarrierControlTests(unittest.TestCase):
    def test_requested_parent_syntax_and_early_missing_data_checks(self):
        with tempfile.TemporaryDirectory() as directory:
            p=Path(directory);path=p/'input'
            path.write_text('$charge-carrier\ndo=0\nget_dimers_monomers=1\none-mol=1\ntwo-mol=1\nreorganization-energy=1\nkmc=1\n$end\n')
            self.assertEqual(load_input(path).plan()['tasks'],[])
            path.write_text(path.read_text().replace('do=0','do=1'))
            with self.assertRaisesRegex(ValueError,'settings are missing'):load_input(path).execution_plan()
            path.write_text('$charge-carrier\ndo=1\nkmc=1\n$end\n$kmc\n$end\n')
            with self.assertRaisesRegex(ValueError,'needs manifest'):load_input(path).execution_plan()

    def test_closed_kmc_ignores_unused_paths_and_strict_settings(self):
        with tempfile.TemporaryDirectory() as directory:
            p=Path(directory);path=p/'input'
            path.write_text('$charge-carrier\ndo=0\nkmc=0\n$end\n$kmc\ndo=0\ninput=missing.json\nbackend=unused\n$end\n')
            self.assertEqual(load_input(path).plan()['tasks'],[])
            for settings in ({'samples':'4'},{'duration':'nan'},{'fit-start':'1'}, {'backend':'gpu'},{'seed':'-1'}):
                with self.assertRaises(ValueError):normalize(settings,p)


@unittest.skipUnless(HAS_NUMPY,'install transport extra for KMC physics tests')
class CarrierPhysicsTests(unittest.TestCase):
    def setUp(self):
        import numpy as np
        self.np=np
        self.data=json.loads((ROOT/'examples/charge-carrier/chain.json').read_text())

    def test_periodic_self_image_diffusion_matches_analytic_and_marcus(self):
        from universal_kmc.transport.network import build_rates,exact_diffusion
        from universal_kmc.rates import marcus
        ids,edges,rows,pi,rank=build_rates(self.data,'hole',300.)
        reference=exact_diffusion(len(ids),edges,pi)
        expected=marcus(.01,.3,0.,300.)*1e-18
        self.np.testing.assert_allclose(reference['diffusion_tensor_m2_s'],self.np.diag([expected,0,0]),rtol=1e-13,atol=1e-30)
        self.assertEqual(rank,1)
        self.np.testing.assert_allclose(edges[0].displacement,[1e-9,0,0],atol=1e-25)
        self.np.testing.assert_array_equal(edges[1].displacement,[-1e-9,0,0])

    def dimer(self,periodic=False):
        data=copy.deepcopy(self.data)
        data['sites'].append({'id':'B','position_angstrom':[5,0,0],'free_energy_ev':{'hole':.02,'electron':-.01}})
        pair=data['pairs'][0];pair.update(id='AB',target='B',image=[0,0,0])
        if periodic:
            other=copy.deepcopy(pair);other.update(id='AB_img_minus_a',image=[-1,0,0]);data['pairs'].append(other)
        return data

    def test_disordered_closed_dimer_has_detailed_balance_and_zero_diffusion(self):
        from universal_kmc.transport.network import build_rates,exact_diffusion
        from universal_kmc.rates import KB
        ids,edges,rows,pi,rank=build_rates(self.dimer(),'hole',300.)
        self.assertAlmostEqual(edges[0].rate/edges[1].rate,math.exp(-.02/(KB*300)),places=13)
        self.assertAlmostEqual(pi[1]/pi[0],math.exp(-.02/(KB*300)),places=13)
        reference=exact_diffusion(2,edges,pi)
        self.np.testing.assert_allclose(reference['diffusion_tensor_m2_s'],0,atol=1e-34)
        self.assertEqual(rank,0)

    def test_parallel_periodic_channels_and_triclinic_vectors_are_preserved(self):
        from universal_kmc.transport.network import build_rates
        data=self.dimer(periodic=True)
        data['lattice_angstrom']=[[10,2,0],[0,10,1],[0,0,10]]
        ids,edges,_,_,rank=build_rates(data,'hole',300.)
        self.assertEqual(len(edges),4);self.assertEqual(rank,1)
        self.np.testing.assert_allclose(edges[2].displacement,[-5e-10,-2e-10,0],atol=1e-25)

    def test_missing_parameters_duplicates_and_disconnection_are_rejected(self):
        from universal_kmc.transport.network import build_rates
        cases=[]
        data=copy.deepcopy(self.data);del data['pairs'][0]['coupling_ev']['hole'];cases.append(data)
        data=copy.deepcopy(self.data);p=copy.deepcopy(data['pairs'][0]);p['id']='reverse';p['image']=[-1,0,0];data['pairs'].append(p);cases.append(data)
        data=self.dimer();data['pairs'][0]['coupling_ev']['hole']=0.;cases.append(data)
        data=copy.deepcopy(self.data);data['pairs'][0]['reorganization_ev']['hole']=-1;cases.append(data)
        for data in cases:
            with self.assertRaises(ValueError):build_rates(data,'hole',300.)

    def test_piecewise_constant_observations_do_not_interpolate_future_hops(self):
        from universal_kmc.models import Edge
        from universal_kmc.transport.simulation import simulate
        class Fixed:
            def __init__(self,values):self.values=iter(values)
            def uniform(self):return next(self.values)
        settings=normalize({'trajectories':'1','duration':'1','samples':'5','max-events':'20','trace-trajectories':'1'},ROOT)
        with patch('universal_kmc.transport.simulation.IndependentRNG',side_effect=[Fixed([.5]),Fixed([math.exp(-.3),.2,math.exp(-2.)])]),redirect_stdout(io.StringIO()):
            grid,positions,traces,counts,capped=simulate(1,[Edge(0,0,1.,(1.,0,0))],[1.],settings,'hole')
        self.np.testing.assert_array_equal(positions[0,:,0],[0,0,1,1,1])
        self.assertEqual(counts,[1]);self.assertEqual(capped,[])

    def test_sparse_outgoing_adapter_preserves_generic_engine_trajectory(self):
        from universal_kmc import Engine,Graph
        from universal_kmc.transport.network import build_rates
        from universal_kmc.transport.simulation import OutgoingGraph,IndependentRNG
        ids,edges,_,_,_=build_rates(self.dimer(periodic=True),'hole',300.)
        adjacency=[[] for _ in ids]
        for i,e in enumerate(edges):adjacency[e.source].append((i,e))
        sparse=OutgoingGraph(adjacency,0);dense=Graph(len(ids),edges)
        a=Engine(sparse,rng=IndependentRNG(2345));b=Engine(dense,rng=IndependentRNG(2345))
        for _ in range(300):
            source=sparse.site
            x,y=a.step(),b.step()
            self.assertEqual(adjacency[source][x.event][0],y.event)
            self.assertEqual(x.state,y.state);self.assertEqual(x.position,y.position);self.assertEqual(x.time,y.time)

    def test_ensemble_diffusion_agrees_with_analytic_within_uncertainty(self):
        from universal_kmc.transport.network import build_rates,exact_diffusion
        from universal_kmc.transport.simulation import simulate,analyze
        ids,edges,_,pi,rank=build_rates(self.data,'hole',300.)
        settings=normalize({'trajectories':'512','duration':'2e-10','samples':'31','seed':'37821','trace-trajectories':'0'},ROOT)
        with redirect_stdout(io.StringIO()):grid,positions,_,_,capped=simulate(len(ids),edges,pi,settings,'hole')
        self.assertEqual(capped,[])
        result=analyze(grid,positions,300.,.5,.5,rank)
        actual=result['diffusion_tensor_m2_s'][0][0];sem=result['diffusion_tensor_sem_m2_s'][0][0]
        exact=exact_diffusion(1,edges,pi)['diffusion_tensor_m2_s'][0][0]
        self.assertLess(abs(actual-exact),4*sem)
        self.assertGreater(sem,0)


@unittest.skipUnless(HAS_NUMPY,'install transport extra for task integration tests')
class CarrierWorkflowTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
        (self.root/'network.json').write_text((ROOT/'examples/charge-carrier/chain.json').read_text())
        self.path=self.root/'input'

    def control(self,parent=True,standalone=False,cap=1000000):
        self.path.write_text('$charge-carrier\ndo='+str(int(parent))+'\nget_dimers_monomers=0\none-mol=0\ntwo-mol=0\nreorganization-energy=0\nkmc=1\n$end\n'
            +'$kmc\n'+('do=1\n' if standalone else '')+'input=network.json\noutput=result\ncarriers=hole\ntrajectories=32\nduration=1e-10\nsamples=11\nrelative-tolerance=10\nmax-events='+str(cap)+'\n$end\n')
        return load_input(self.path)

    def test_parent_and_standalone_kmc_use_same_task_and_audited_resume(self):
        from universal_kmc.workflow.runner import run_workflow
        control=self.control(parent=True,standalone=True)
        self.assertEqual(control.plan()['tasks'],['kmc'])
        self.assertEqual(control.plan()['module_owners']['kmc'],['charge-carrier','standalone'])
        with redirect_stdout(io.StringIO()):
            self.assertEqual(run_workflow(control),0)
            self.assertEqual(run_workflow(control,resume=True),0)
        report=json.loads((self.root/'cube-results/workflow.json').read_text())
        self.assertEqual(report['calculations']['charge-carrier']['status'],'completed')
        self.assertTrue(report['tasks']['kmc']['outputs']['network']['reused'])
        (self.root/'result/hole.kmc').write_text('tampered')
        with self.assertRaisesRegex(ValueError,'missing/changed'):run_workflow(control,resume=True)

    def test_standalone_without_active_parent_and_event_cap_partial(self):
        from universal_kmc.workflow.runner import run_workflow
        control=self.control(parent=False,standalone=True,cap=1)
        with redirect_stdout(io.StringIO()):self.assertEqual(run_workflow(control),2)
        report=json.loads((self.root/'cube-results/workflow.json').read_text())
        self.assertEqual(report['calculations']['charge-carrier']['status'],'disabled')
        result=json.loads((self.root/'result/result.json').read_text())
        self.assertEqual(result['carriers']['hole']['estimate']['reason'],'event_cap_reached')
        self.assertNotIn('diffusion_tensor_m2_s',result['carriers']['hole']['estimate'])

    def test_output_directory_cannot_overwrite_network_input(self):
        from universal_kmc.workflow.runner import run_workflow
        self.control()
        original=(self.root/'network.json').read_bytes()
        self.path.write_text(self.path.read_text().replace('output=result','output=.'))
        with self.assertRaisesRegex(ValueError,'overwrite a source file'):
            run_workflow(load_input(self.path))
        self.assertEqual((self.root/'network.json').read_bytes(),original)

    def test_independent_kmc_does_not_complete_a_parent_that_deselected_it(self):
        from universal_kmc.workflow.runner import run_workflow
        self.control(parent=True,standalone=True,cap=1)
        self.path.write_text(self.path.read_text().replace('kmc=1','kmc=0'))
        with redirect_stdout(io.StringIO()):self.assertEqual(run_workflow(load_input(self.path)),2)
        report=json.loads((self.root/'cube-results/workflow.json').read_text())
        self.assertEqual(report['calculations']['charge-carrier']['status'],'prepared')
        self.assertEqual(report['tasks']['kmc']['status'],'partial')

    def test_qc_assembly_uses_material_origins_and_checks_pair_reorganization(self):
        from universal_kmc.transport.assembly import assemble
        manifest={'schema':'cube-dimers-1','coordinate_unit':'Angstrom','lattice_angstrom':[[10,0,0],[0,10,0],[0,0,10]],
            'monomers':[{'id':str(i),'coordinates_angstrom':[[0,0,0]],'transport_position_angstrom':[5*(i-1),0,0]} for i in (1,2)],
            'dimers':[{'id':'1_2','mol_i':'1','mol_j':'2','image':[0,0,0]}]}
        (self.root/'manifest.json').write_text(json.dumps(manifest))
        (self.root/'couplings.csv').write_text('dimer,mol_i,mol_j,VH_meV\n1_2,1,2,10\n')
        paths={}
        for i in (1,2):
            data={'schema':'cube-reorganization-1','status':'completed','results':{'hole':{
                'relaxation_ev':{'charged_relaxation':.18,'neutral_relaxation':.12},
                'energies_hartree':{'Eq_Rq':-9.,'E0_R0':-10.}}}}
            path=self.root/f'reorg{i}.json';path.write_text(json.dumps(data));paths[str(i)]=str(path)
        settings={'carriers':['hole'],'site-energy-model':'uniform'};sources={}
        data=assemble(self.root/'manifest.json',self.root/'couplings.csv',paths,settings,sources)
        self.assertEqual(data['sites'][1]['position_angstrom'],[5,0,0])
        self.assertAlmostEqual(data['pairs'][0]['coupling_ev']['hole'],.01)
        self.assertAlmostEqual(data['pairs'][0]['reorganization_ev']['hole'],.3)
        # Exercise the public control path with all expensive preparation stages
        # disabled and existing quantum outputs, rather than a prebuilt network.
        manifest['dimers'].append({'id':'1_2_img_m1','mol_i':'1','mol_j':'2','image':[-1,0,0]})
        (self.root/'manifest.json').write_text(json.dumps(manifest))
        with (self.root/'couplings.csv').open('a') as f:f.write('1_2_img_m1,1,2,10\n')
        (self.root/'reorg-index.json').write_text(json.dumps(paths))
        self.path.write_text('$charge-carrier\ndo=1\nget_dimers_monomers=0\none-mol=0\ntwo-mol=0\nreorganization-energy=0\nkmc=1\n$end\n'
            '$kmc\nmanifest=manifest.json\ncouplings=couplings.csv\nreorganization=reorg-index.json\nsite-energy-model=uniform\ncarriers=hole\n'
            'output=assembled\ntrajectories=32\nduration=1e-10\nsamples=11\nrelative-tolerance=10\n$end\n')
        from universal_kmc.workflow.runner import run_workflow
        with redirect_stdout(io.StringIO()):self.assertEqual(run_workflow(load_input(self.path)),0)
        assembled=json.loads((self.root/'assembled/network.json').read_text())
        self.assertEqual(assembled['sites'][1]['position_angstrom'],[5,0,0])
        self.assertEqual(len(assembled['pairs']),2)
        path=self.root/'reorg2.json';bad=json.loads(path.read_text());bad['results']['hole']['relaxation_ev']['charged_relaxation']=.3;path.write_text(json.dumps(bad))
        with self.assertRaisesRegex(ValueError,'inequivalent'):assemble(self.root/'manifest.json',self.root/'couplings.csv',paths,settings,{})


if __name__=='__main__':unittest.main()
