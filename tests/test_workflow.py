"""Section controls and scientifically meaningful extraction/workflow checks."""
from contextlib import redirect_stdout
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

from universal_kmc.workflow.control import load_input

try:
    import numpy as np
except ImportError:
    np = None


TEMPLATE = '''$get_dimers_monomers
tpr = a.tpr
gro = a.gro
xtc = a.xtc
frame =
dimers_cutoff = 6 # Angstrom
fragment-file = core.txt
output = mol
format = xyz,gjf
$end
$two-mol
method = b3lyp
basis-set = 6-31G(d,p)
transfer-integrals-method = DIPRO
$end
$one-mol
qc-engine = gpu-pyscf
method = b3lyp
basis-set = 6-31G(d,p)
$end
$reorganization-energy
method = b3lyp
basis-set = 6-31G(d,p)
re-method = 4-points
$end
'''


class SectionControlTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        for name in ('a.tpr','a.gro','a.xtc','core.txt'):
            (self.root/name).write_text('1 2 3\n')
        self.path = self.root/'input'
        self.path.write_text(TEMPLATE)

    def test_exact_user_syntax_defaults_engine_inheritance_and_order(self):
        control = load_input(self.path)
        self.assertEqual(control.plan()['tasks'], ['get_dimers_monomers','one-mol','two-mol','reorganization-energy'])
        self.assertEqual(control.sections['get_dimers_monomers']['frame'], -1)
        self.assertEqual(control.sections['get_dimers_monomers']['format'], ['xyz','gjf'])
        self.assertEqual(control.sections['get_dimers_monomers']['dimers_cutoff'], 6.)
        self.assertEqual(control.sections['two-mol']['basis-set'], '6-31G(d,p)')
        self.assertEqual(control.sections['two-mol']['qc-engine'], 'gpu-pyscf')
        self.assertEqual(control.sections['reorganization-energy']['environment'], 'vacuum')
        self.assertEqual(control.sections['reorganization-energy']['qc-engine'], 'gpu-pyscf')

    def test_plan_is_dependency_free_and_resolves_paths_against_input(self):
        script = '''import sys
from universal_kmc.workflow.cli import main
assert main([sys.argv[1],'--plan']) == 0
assert all(x not in sys.modules for x in ('numpy','pyscf','cupy','MDAnalysis'))
'''
        env = dict(os.environ, PYTHONPATH=str(Path(__file__).resolve().parents[1]/'src'))
        p = subprocess.run([sys.executable,'-c',script,str(self.path)], cwd='/tmp', env=env, text=True,capture_output=True)
        self.assertEqual(p.returncode, 0, p.stderr)
        self.assertEqual(json.loads(p.stdout)['sections']['get_dimers_monomers']['tpr'], str(self.root/'a.tpr'))
        self.assertFalse((self.root/'mol').exists())

    def test_invalid_or_ambiguous_controls_rejected(self):
        changes = [('$two-mol','$typo'), ('dimers_cutoff = 6','cutoff = nan'), ('frame =\n','frame = -2\n'),
                   ('format = xyz,gjf','format = xyz,XYZ'), ('format = xyz,gjf','format = pdb'),
                   ('qc-engine = gpu-pyscf','qc-engine = gaussian'),
                   ('re-method = 4-points','re-method = approximate'),
                   ('re-method = 4-points','re-method = 4-points\nenvironment = qmmm\nmm-frozen = false'),
                   ('tpr = a.tpr','tpr ='), ('dimers_cutoff = 6','dimers_cutoff = 6\ndimers_cutoff = 7'),
                   ('$get_dimers_monomers','$get_dimers_monomers\nunknown = 3'),
                   ('basis-set = 6-31G(d,p)\ntransfer', 'basis-set = sto-3g\ntransfer')]
        for old,new in changes:
            with self.subTest(new=new):
                self.path.write_text(TEMPLATE.replace(old,new))
                with self.assertRaises(ValueError):
                    load_input(self.path)
        for content in (TEMPLATE+'$one-mol\n$end\n', TEMPLATE.rsplit('$end',1)[0],
                        '$end\n', 'dimers_cutoff = 6\n', '$get_dimers_monomers\n$one-mol\n$end\n'):
            self.path.write_text(content)
            with self.assertRaises(ValueError):
                load_input(self.path)

    def test_hash_in_quoted_path_and_standalone_xyz(self):
        xyz = self.root/'a # molecule.xyz'
        xyz.write_text('2\ncomment\nH 0 0 0\nH .7 0 0\n')
        self.path.write_text('$one-mol\nxyz = "a # molecule.xyz" # comment\nqc-engine = pyscf\n$end\n')
        self.assertEqual(load_input(self.path).sections['one-mol']['xyz'], [str(xyz)])

    def test_sampling_cluster_defaults_and_legacy_aliases(self):
        source = TEMPLATE.replace('frame =\n', 'frame = 5\nframe_step = 20\nmonomers-qmmm-cutoff =\n')
        self.path.write_text(source)
        s = load_input(self.path).sections['get_dimers_monomers']
        self.assertEqual((s['frame'], s['frame_step'], s['monomers-qmmm-cutoff']), (5, 20, 12.))
        for key, value in [('frame_step','0'), ('frame_step','-2'), ('frame_step','1.5'),
                           ('monomers-qmmm-cutoff','nan'), ('monomers-qmmm-cutoff','0')]:
            self.path.write_text(TEMPLATE.replace('frame =\n', f'frame = -1\n{key} = {value}\n'))
            with self.subTest(key=key, value=value), self.assertRaises(ValueError):
                load_input(self.path)
        self.path.write_text(TEMPLATE.replace('$get_dimers_monomers', '$get_dimer').replace('dimers_cutoff', 'cutoff'))
        c = load_input(self.path)
        self.assertEqual(c.plan(['get_dimer'])['tasks'], ['get_dimers_monomers'])
        self.assertEqual(c.sections['get_dimers_monomers']['dimers_cutoff'], 6.)
        self.path.write_text(TEMPLATE.replace('re-method = 4-points', 're-method = 4-points\nenvironment = qmmm')
                             .replace('frame =\n', 'frame = -1\nmonomers-qmmm-cutoff = 9\n'))
        self.assertEqual(load_input(self.path).sections['reorganization-energy']['mm-cutoff'], 9.)


@unittest.skipIf(np is None, 'NumPy optional dependency required')
class PeriodicCoreTests(unittest.TestCase):
    def test_frame_selection_starts_at_requested_frame_and_does_not_append_last(self):
        from universal_kmc.workflow.dimers import selected_frames
        self.assertEqual(selected_frames(-1, 20, 46), [45])
        self.assertEqual(selected_frames(5, 20, 46), [5, 25, 45])
        self.assertEqual(selected_frames(0, 20, 46), [0, 20, 40])
        with self.assertRaises(ValueError):
            selected_frames(46, 20, 46)

    def test_boundary_contact_whole_image_and_strict_cutoff(self):
        from universal_kmc.workflow.dimers import contact_images
        cell = np.diag([10.,10.,10.])
        cores = [np.array([[9.8,0,0],[9.9,.1,0]]), np.array([[.2,0,0],[.3,.1,0]])]
        contacts = list(contact_images(cores,cell,.5))
        self.assertEqual(len(contacts),1)
        self.assertEqual(contacts[0][2].tolist(),[1,0,0])
        self.assertAlmostEqual(contacts[0][3],np.sqrt(.1**2+.3**2))
        self.assertFalse(list(contact_images([np.array([[0.,0,0]]),np.array([[1.,0,0]])],cell,1.)))

    def test_triclinic_images_match_exhaustive_search(self):
        from universal_kmc.workflow.dimers import contact_images
        cell = np.array([[8.,0,0],[3.,7.,0],[2.,1.,9.]])
        cores = [np.array([[.1,.2,.3],[.5,.3,.2]]),np.array([[7.8,6.8,8.7],[8.1,6.7,8.5]])]
        import itertools
        expected = {}
        for shift in itertools.product(range(-3,4),repeat=3):
            distances = np.linalg.norm(cores[0][:,None,:]-(cores[1]+np.array(shift)@cell)[None,:,:],axis=2)
            if distances.min()<5:
                expected[shift] = float(distances.min())
        actual = {tuple(shift):d for i,j,shift,d in contact_images(cores,cell,5)}
        self.assertEqual(actual.keys(),expected.keys())
        for shift in actual:
            self.assertAlmostEqual(actual[shift],expected[shift])

    def test_multiple_images_are_preserved(self):
        from universal_kmc.workflow.dimers import contact_images
        contacts = list(contact_images([np.array([[0.,0,0]]),np.array([[1.,0,0]])],np.diag([2.,5.,5.]),1.1))
        self.assertEqual(len(contacts),2)
        self.assertEqual({tuple(c[2]) for c in contacts},{(0,0,0),(-1,0,0)})

    def test_four_point_formula_preserves_relaxation_signs(self):
        from universal_kmc.workflow.reorganization import four_point_terms
        terms = four_point_terms(-10.,-9.,-9.2,-9.8)
        self.assertAlmostEqual(terms['total'],.4)
        self.assertLess(four_point_terms(-10.,-9.3,-9.2,-9.8)['charged_relaxation'],0)
        with self.assertRaises(ValueError):
            four_point_terms(float('nan'),0,0,0)


try:
    import MDAnalysis as mda
except ImportError:
    mda = None

try:
    import pyscf
except ImportError:
    pyscf = None


@unittest.skipIf(pyscf is None or np is None, 'quantum dependencies required')
class CheckpointMappingTests(unittest.TestCase):
    def test_multiframe_quantum_dispatch_separates_outputs_and_reuses_checkpoints(self):
        from universal_kmc.workflow.dimers import file_hash
        from universal_kmc.workflow.runner import run_workflow, _dataset_sources
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            entries = []
            for frame, separation in [(0, 3.), (2, 4.)]:
                here = root/f'frame{frame:04d}'
                here.mkdir()
                monomers = [{'id': str(i), 'elements': ['H','H'],
                    'coordinates_angstrom': [[-.37,0,0],[.37,0,0]]} for i in (1,2)]
                pair = {'id':'1_2','mol_i':'1','mol_j':'2','elements':['H']*4, 'atom_count_a':2,
                        'coordinates_angstrom':[[-.37,0,-separation/2],[.37,0,-separation/2],
                                               [-.37,0,separation/2],[.37,0,separation/2]]}
                manifest = here/'manifest.json'
                manifest.write_text(json.dumps({'schema':'cube-dimers-1','coordinate_unit':'Angstrom',
                    'frame':frame,'monomers':monomers,'dimers':[pair]}))
                entries.append({'frame':frame,'manifest':str(manifest.relative_to(root)),'sha256':file_hash(manifest)})
            index = root/'frames.json'
            index.write_text(json.dumps({'schema':'cube-frames-1','frame_count':2,'frames':entries}))
            control = root/'input'
            common = 'input = frames.json\nmethod = hf\nbasis-set = sto-3g\nthreads = 1\n'
            control.write_text('$one-mol\n'+common+'output = one\n$end\n$two-mol\n'+common+'output = two\n$end\n')
            with redirect_stdout(io.StringIO()):
                self.assertEqual(run_workflow(load_input(control)), 0)
            report = json.loads((root/'cube-results/workflow.json').read_text())
            self.assertEqual(report['tasks']['two-mol']['frame_count'], 2)
            self.assertEqual(report['tasks']['one-mol']['count'], 4)
            magnitudes = []
            for label in ('frame0000','frame0002'):
                transfer = json.loads((root/'two'/label/'1_2/transfer.json').read_text())
                self.assertEqual(transfer['status'], 'completed')
                magnitudes.append(transfer['results']['hole']['methods']['dipro']['scalar_magnitude_ev'])
            self.assertGreater(magnitudes[0], magnitudes[1])
            checkpoints = {p:p.stat().st_mtime_ns for p in root.rglob('*.chk')}
            with redirect_stdout(io.StringIO()):
                self.assertEqual(run_workflow(load_input(control), resume=True), 0)
            self.assertTrue(all(p.stat().st_mtime_ns == t for p,t in checkpoints.items()))
            (root/'frame0000/manifest.json').write_text('{}')
            with self.assertRaisesRegex(ValueError, 'missing/changed'):
                _dataset_sources({'input':str(index)}, None)

    def test_translation_preserves_complete_mo_and_fock_and_rejects_rotation(self):
        from pyscf import gto
        from pyscf.scf import chkfile
        from universal_kmc.workflow.scf import translated_checkpoint
        from universal_kmc.workflow.dimers import file_hash
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            mol=gto.M(atom='H 0 0 0; H .74 0 0',basis='sto-3g',verbose=0)
            s=mol.intor_symmetric('int1e_ovlp')
            eigenvectors, u=np.linalg.eigh(s)
            c=(u/np.sqrt(eigenvectors))@u.T
            e=np.array([-.5,.2])
            f=(s@c*e)@c.T@s
            source=root/'source.chk'
            chkfile.dump_scf(mol,str(source),-1.,e,c,np.array([2.,0.]))
            chkfile.save(str(source),'cube_transfer',{'schema_version':1,'converged':True,'model':'hf','fock_hartree':f})
            original_hash=file_hash(source)
            shift=np.array([10.,-4.,3.])
            atoms=[(mol.atom_symbol(i),tuple(mol.atom_coord(i,unit='Angstrom')+shift)) for i in range(2)]
            target=translated_checkpoint(source,atoms,root/'translated.chk')
            moved,data=chkfile.load_scf(str(target))
            np.testing.assert_allclose(moved.atom_coords(unit='Angstrom'),mol.atom_coords(unit='Angstrom')+shift,atol=1e-12)
            np.testing.assert_array_equal(data['mo_coeff'],c)
            np.testing.assert_array_equal(chkfile.load(str(target),'cube_transfer')['fock_hartree'],f)
            self.assertEqual(file_hash(source),original_hash)
            with self.assertRaisesRegex(ValueError,'pure translation'):
                translated_checkpoint(source,[('H',(0.,0.,0.)),('H',(0.,.74,0.))],root/'rotated.chk')

    def test_failed_optimization_never_reports_successful_lambda(self):
        from universal_kmc.workflow.reorganization import calculate_reorganization
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            (root/'H2.xyz').write_text('2\ncomment\nH 0 0 0\nH .74 0 0\n')
            (root/'input').write_text('$reorganization-energy\nxyz = H2.xyz\nqc-engine = pyscf\n$end\n')
            settings=load_input(root/'input').sections['reorganization-energy']
            with patch('universal_kmc.workflow.reorganization._optimize',side_effect=RuntimeError('not converged')):
                with self.assertRaises(RuntimeError):
                    calculate_reorganization([('H',(0.,0.,0.)),('H',(.74,0.,0.))],settings,root/'reorg')
            result=json.loads((root/'reorg/reorganization.json').read_text())
            self.assertEqual(result['status'],'failed')
            self.assertFalse(result['results'])


@unittest.skipIf(mda is None, 'workflow MDAnalysis dependency required')
class ExtractionTests(unittest.TestCase):
    def universe(self):
        u = mda.Universe.empty(6,n_residues=2,atom_resindex=[0,0,0,1,1,1],trajectory=True)
        for name,values in (('names',['C','H','H']*2),('elements',['C','H','H']*2),
                            ('masses',[12.,1.,1.]*2),('charges',[0.]*6)):
            u.add_TopologyAttr(name,values)
        u.add_TopologyAttr('bonds',[(0,1),(0,2),(3,4),(3,5)])
        u.dimensions=[10.,10.,10.,90.,90.,90.]
        u.atoms.positions=[[9.8,0,0],[.1,.1,0],[9.7,.2,0],[.2,0,0],[.4,.1,0],[.2,.2,0]]
        return u

    def test_format_selection_pbc_and_verified_resume(self):
        from universal_kmc.workflow.dimers import extract_dimers
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            for name in ('a.tpr','a.gro','a.xtc'):
                (root/name).write_bytes(b'fixture')
            (root/'core.txt').write_text('1 # local carbon\n')
            settings={k:str(root/v) for k,v in {'tpr':'a.tpr','gro':'a.gro','xtc':'a.xtc','fragment-file':'core.txt','output':'mol'}.items()}
            settings.update(frame=0,cutoff=1.,format=['gjf'])
            with patch.object(mda,'Universe',side_effect=[self.universe(),self.universe()]):
                manifest_path=extract_dimers(settings,one_settings={'method':'hf','basis-set':'sto-3g'})
            manifest=json.loads(manifest_path.read_text())
            self.assertEqual(manifest['molecule_count'],2)
            self.assertEqual(manifest['dimer_count'],1)
            self.assertEqual(manifest['dimers'][0]['image'],[1,0,0])
            self.assertFalse(list(manifest_path.parent.rglob('*.xyz')))
            self.assertEqual(len(list(manifest_path.parent.rglob('*.gjf'))),3)
            self.assertIn('hf/sto-3g',(manifest_path.parent/'monomers/1.gjf').read_text())
            with patch.object(mda,'Universe',side_effect=[self.universe(),self.universe()]):
                self.assertEqual(extract_dimers(settings,one_settings={'method':'hf','basis-set':'sto-3g','molecules':['1']},resume=True),manifest_path)
            (manifest_path.parent/'monomers/1.gjf').write_text('changed')
            with patch.object(mda,'Universe',side_effect=[self.universe(),self.universe()]),self.assertRaisesRegex(ValueError,'changed'):
                extract_dimers(settings,one_settings={'method':'hf','basis-set':'sto-3g'},resume=True)


    def test_multiple_frames_full_cluster_mapping_formats_and_resume_integrity(self):
        from universal_kmc.workflow.dimers import extract_dimers_monomers
        from universal_kmc.workflow.embedding import read_charges
        def trajectory():
            u = self.universe()
            xyz = u.atoms.positions.copy()
            frames = np.stack([xyz + [0, i*.1, 0] for i in range(5)])
            frames[:, 3:, 0] += 1.+np.arange(5)[:,None]*.1
            u.load_new(frames, dimensions=np.tile([10.,10.,10.,90.,90.,90.], (5,1)))
            return u
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for name in ('a.tpr','a.gro','a.xtc'):
                (root/name).write_bytes(b'fixture')
            (root/'core.txt').write_text('1\n')
            settings = {k:str(root/v) for k,v in {'tpr':'a.tpr','gro':'a.gro','xtc':'a.xtc',
                        'fragment-file':'core.txt','output':'mol'}.items()}
            settings.update(frame=1, frame_step=2, dimers_cutoff=2., format=['xyz','gjf'], **{'monomers-qmmm-cutoff': 2.})
            with patch.object(mda,'Universe',side_effect=[trajectory(), self.universe()]):
                index_path = extract_dimers_monomers(settings)
            index = json.loads(index_path.read_text())
            self.assertEqual([f['frame'] for f in index['frames']], [1,3])
            self.assertFalse((root/'mol/frame0004').exists())
            cluster_files = []
            for f in index['frames']:
                frame_dir = index_path.parent/Path(f['manifest']).parent
                manifest = json.loads((frame_dir/'manifest.json').read_text())
                self.assertEqual(manifest['frame'], f['frame'])
                cluster = json.loads((frame_dir/'monomers/1_cluster.json').read_text())
                self.assertEqual((cluster['qm_atom_count'], cluster['mm_atom_count']), (3,3))
                mapping = cluster['mm']['provenance']['molecular_images']
                self.assertEqual(mapping[0]['global_atom_indices_zero_based'], [3,4,5])
                self.assertEqual(mapping[0]['image'], [1,0,0])
                self.assertEqual(mapping[0]['molecule_id'], '2')
                qm = np.asarray(cluster['qm']['coordinates_angstrom'])
                mm = np.asarray(cluster['mm']['coordinates_angstrom'])
                pc = read_charges(frame_dir/'monomers/1_mm.pc')
                np.testing.assert_allclose(pc.coordinates_angstrom, mm, atol=1e-14)
                self.assertEqual(int((frame_dir/'monomers/1.xyz').read_text().splitlines()[0]), 3)
                self.assertEqual(int((frame_dir/'monomers/1_cluster.xyz').read_text().splitlines()[0]), 6)
                gjf = (frame_dir/'monomers/1_cluster.gjf').read_text()
                self.assertIn(' SP Charge ', gjf)
                self.assertEqual(sum(line.startswith(('C ', 'H ')) for line in gjf.splitlines()), 3)
                self.assertAlmostEqual(mm[0,0]-qm[0,0], 1.4+f['frame']*.1, places=5)
                cluster_files.append(frame_dir/'monomers/1_cluster.xyz')
            self.assertNotEqual(cluster_files[0].read_bytes(), cluster_files[1].read_bytes())
            times = {p:p.stat().st_mtime_ns for p in cluster_files}
            with patch.object(mda,'Universe',side_effect=[trajectory(), self.universe()]):
                self.assertEqual(extract_dimers_monomers(settings,resume=True), index_path)
            self.assertTrue(all(p.stat().st_mtime_ns == t for p,t in times.items()))
            (cluster_files[0].parent/'1_mm.pc').write_text('0 0 0 99\n')
            with patch.object(mda,'Universe',side_effect=[trajectory(), self.universe()]), self.assertRaisesRegex(ValueError,'changed'):
                extract_dimers_monomers(settings,resume=True)
            settings.update(frame=-1, frame_step=20, output=str(root/'last'), format=['gjf'], **{'monomers-qmmm-cutoff': .01})
            with patch.object(mda,'Universe',side_effect=[trajectory(), self.universe()]):
                last = extract_dimers_monomers(settings)
            self.assertEqual(json.loads(last.read_text())['frame'], 4)
            self.assertFalse(list(last.parent.rglob('*.xyz')))
            empty = json.loads((last.parent/'monomers/1_cluster.json').read_text())
            self.assertEqual(empty['mm_atom_count'], 0)
            self.assertNotIn(' SP Charge ', (last.parent/'monomers/1_cluster.gjf').read_text())


if __name__ == '__main__':
    unittest.main()
