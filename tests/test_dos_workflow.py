from contextlib import redirect_stdout
import csv
import importlib.util
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from wren.dos.control import normalize


class DosControlTests(unittest.TestCase):
    @unittest.skipUnless(importlib.util.find_spec('numpy'), 'DOS math needs NumPy')
    def test_mesh_study_edges_and_full_spectrum_change(self):
        import numpy as np
        from wren.dos.mesh_study import shape_change
        from wren.dos.periodic import band_edges
        bands = np.array([[-2., 2.], [-1., 3.]])
        occ = np.array([[2., 0.], [2., 0.]])
        edges = band_edges(bands, occ)
        self.assertEqual((edges['vbm_k_index'], edges['cbm_k_index']), (1, 0))
        self.assertEqual((edges['sampled_gap_ev'], edges['minimum_sampled_direct_gap_ev']), (3., 4.))
        self.assertFalse(edges['partial_occupation'])
        self.assertIsNone(band_edges(bands, np.full_like(occ, 2.))['sampled_gap_ev'])
        self.assertEqual(band_edges(bands, np.ones_like(occ))['sampled_gap_ev'], 0.)
        base = {'orbital_energies_ev': [[0., 1.]], 'kpoint_weights': [1.]}
        doubled = {'orbital_energies_ev': [[0., 1.]]*2, 'kpoint_weights': [.5, .5]}
        self.assertAlmostEqual(shape_change(base, doubled, .1), 0.)
        shifted = dict(base, orbital_energies_ev=[[3., 4.]])
        self.assertAlmostEqual(shape_change(base, shifted, .1), 1., places=8)

    def test_check_is_dependency_free_and_paths_are_relative_to_control(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root/'h2.xyz').write_text('2\nH2\nH 0 0 0\nH .74 0 0\n')
            (root/'input').write_text('$dos\ndo=1\nxyz=h2.xyz\n$end\n')
            script = '''import sys
from wren.workflow.cli import main
assert main([sys.argv[1], '--plan']) == 0
assert all(name not in sys.modules for name in ('numpy','pyscf','cupy','MDAnalysis'))
'''
            env = dict(os.environ, PYTHONPATH=str(Path(__file__).resolve().parents[1]/'src'))
            process = subprocess.run([sys.executable, '-c', script, str(root/'input')], cwd='/tmp',
                                     env=env, text=True, capture_output=True)
            self.assertEqual(process.returncode, 0, process.stderr)
            self.assertEqual(json.loads(process.stdout)['sections']['dos']['xyz'], str(root/'h2.xyz'))
            with self.assertRaises(ValueError):
                normalize({'xyz': 'h2.xyz', 'sigma': 0}, root)
            periodic = normalize({'mode': 'periodic', 'input': 'h2.xyz', 'qc-engine': 'gpu-pyscf'}, root)
            self.assertEqual(periodic['qc-engine'], 'gpu-pyscf')
            for mode in ('vertical-dscf', 'tight-binding'):
                with self.assertRaisesRegex(ValueError, 'element-projections'):
                    normalize({'mode': mode, 'input': 'h2.xyz', 'element-projections': True}, root)


@unittest.skipUnless(importlib.util.find_spec('numpy'), 'DOS NumPy dependency required')
class DosWorkflowTests(unittest.TestCase):
    @unittest.skipUnless(importlib.util.find_spec('pyscf'), 'PySCF required')
    def test_vertical_dscf_reads_md_frames_keeps_geometry_and_separate_line_densities(self):
        import numpy as np
        from pyscf.data.nist import HARTREE2EV
        from wren.workflow.dimers import file_hash
        from wren.dos.runner import run
        from wren.workflow.cli import main
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            frames = []
            for frame, length in ((2, .74), (4, .8)):
                path = root/f'f{frame}.json'
                path.write_text(json.dumps({'schema': 'cube-dimers-1', 'coordinate_unit': 'Angstrom',
                    'frame': frame, 'time_ps': frame/2, 'monomers': [{'id': 'A', 'elements': ['H', 'H'],
                        'coordinates_angstrom': [[0., 0., 0.], [length, 0., 0.]]}]}))
                frames.append({'frame': frame, 'manifest': path.name, 'sha256': file_hash(path)})
            source = root/'frames.json'
            source.write_text(json.dumps({'schema': 'cube-frames-1', 'frame_count': 2, 'frames': frames}))
            control = root/'input'
            control.write_text('$dos\ndo=1\nmode=vertical-dscf\ninput=frames.json\nmolecules=A\n'
                               'output=out\nmethod=hf\nbasis-set=aug-cc-pvdz\nauxbasis=weigend\nthreads=1\n$end\n')
            with redirect_stdout(io.StringIO()):
                self.assertEqual(main([str(control)]), 0)
            directory = root/'out'
            result = json.loads((directory/'result.json').read_text())
            self.assertEqual(result['request']['source_sha256']['input'], file_hash(source))
            for path, digest in result['request']['implementation_sha256'].items():
                self.assertEqual(digest, file_hash(path))
            self.assertEqual(result['nonpositive_vertical_ea_count'], 2)
            self.assertFalse((directory/'dos.csv').exists())
            with (directory/'vertical-energies.csv').open() as file:
                raw = list(csv.DictReader(file))
            self.assertEqual([row['frame'] for row in raw], ['2', '4'])
            self.assertEqual([row['mol'] for row in raw], ['A', 'A'])
            for row in raw:
                neutral = float(row['neutral_hartree'])
                self.assertAlmostEqual(float(row['vertical_ip_ev']), (float(row['cation_hartree'])-neutral)*HARTREE2EV)
                self.assertAlmostEqual(float(row['vertical_ea_ev']), (neutral-float(row['anion_hartree']))*HARTREE2EV)
            self.assertAlmostEqual(float(raw[0]['vertical_ip_ev']), 15.3375, places=3)
            self.assertAlmostEqual(float(raw[0]['vertical_ea_ev']), -1.6347, places=3)
            for channel, property_name in (('removal', 'vertical_ip_ev'), ('addition', 'vertical_ea_ev')):
                data = np.genfromtxt(directory/(channel+'.csv'), delimiter=',', names=True)
                self.assertAlmostEqual(result['channels'][channel]['finite_grid_integral'], 1., places=10)
                for index, row in enumerate(raw):
                    moment = np.trapezoid(data['energy_ev']*data[f'frame_{index}'], data['energy_ev'])
                    self.assertAlmostEqual(moment, -float(row[property_name]), places=10)
                np.testing.assert_allclose(data['density_per_ev'], (data['frame_0']+data['frame_1'])/2)
            molecule = json.loads((directory/'spectra.json').read_text())['frames'][0]['molecules'][0]
            for state, charge, spin in (('neutral', 0, 0), ('cation', 1, 1), ('anion', -1, 1)):
                status = json.loads((directory/molecule['scf'][state]['status_file']).read_text())
                self.assertEqual((status['charge'], status['spin_2s']), (charge, spin))
                self.assertEqual(status['request']['atoms'], [['H', [0., 0., 0.]], ['H', [.74, 0., 0.]]])
                self.assertEqual(file_hash(status['checkpoint']), molecule['scf'][state]['checkpoint_sha256'])
            self.assertTrue(run(result['request']['settings'], resume=True)['reused'])
            (directory/'addition.csv').write_text('changed')
            with self.assertRaisesRegex(ValueError, 'changed'):
                run(result['request']['settings'], resume=True)

    def test_vertical_dscf_rejects_k_weights_and_spin_multipliers(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root/'h2.xyz').write_text('2\n\nH 0 0 0\nH .74 0 0')
            for option in ({'spin-degeneracy': 2}, {'kmesh': [2, 1, 1]}):
                with self.assertRaisesRegex(ValueError, 'one transition'):
                    normalize({'mode': 'vertical-dscf', 'xyz': 'h2.xyz', **option}, root)

    def test_site_projection_control_outputs_and_resume(self):
        import numpy as np
        from wren.dos.cli import main
        from wren.dos.runner import run
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root/'model.json').write_text(json.dumps({'schema': 'wren-dos-tb-1',
                'coupling_convention': 'signed_real_consistent_orbital_phases',
                'site_energies_ev': [0., 0.], 'pairs': [[0, 1, [0, 0, 0], -.1], [0, 1, [1, 0, 0], -.1]]}))
            control = root/'tb.toml'
            control.write_text('[dos]\nmode="tight-binding"\ninput="model.json"\noutput="out"\n'
                               'spin-degeneracy=2\nkmesh=[4,1,1]\nsite-projections=true\n')
            with redirect_stdout(io.StringIO()):
                self.assertEqual(main([str(control)]), 0)
                self.assertEqual(main([str(control), '--resume']), 0)
            result = json.loads((root/'out/result.json').read_text())
            self.assertAlmostEqual(result['finite_grid_integral'], 4.)
            info = result['site_projections']
            np.testing.assert_allclose(info['finite_grid_integral_per_site'], [2., 2.], atol=1e-12)
            total = np.genfromtxt(root/'out/dos.csv', delimiter=',', names=True)
            sites = np.genfromtxt(root/'out/pdos.csv', delimiter=',', names=True)
            np.testing.assert_array_equal(total['energy_ev'], sites['energy_ev'])
            np.testing.assert_allclose(sites['site_0'] + sites['site_1'], total['dos'], atol=1e-12)
            with (root/'out/participation.csv').open() as file:
                groups = list(csv.DictReader(file))
            self.assertEqual(len(groups), 7)
            degenerate, = [g for g in groups if g['rank'] == '2']
            self.assertEqual((degenerate['k_index'], degenerate['ipr']), ('2', ''))
            self.assertAlmostEqual(float(degenerate['subspace_profile_ipr']), .5)
            (root/'out/pdos.csv').write_text('tampered')
            with self.assertRaisesRegex(ValueError, 'missing or changed'):
                run(result['request']['settings'], resume=True)
            for bad in ({'site-projections': 'maybe'}, {'degeneracy-tol-ev': 1.},
                        {'site-projections': True, 'degeneracy-tol-ev': -1}):
                with self.assertRaises(ValueError):
                    normalize({'mode': 'tight-binding', 'input': 'model.json', **bad}, root)
            with self.assertRaisesRegex(ValueError, 'require tight-binding'):
                normalize({'input': 'model.json', 'site-projections': 'true'}, root)

    def test_bloch_control_routing_and_output_hash_resume(self):
        from wren.workflow.cli import main
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root/'model.json').write_text(json.dumps({'schema': 'wren-dos-tb-1',
                'coupling_convention': 'signed_real_consistent_orbital_phases',
                'site_energies_ev': [0.], 'pairs': [[0, 0, [1, 0, 0], -.1]]}))
            p = root/'input'
            p.write_text('$density-of-states\ndo=1\n$end\n$dos\nmode=tight-binding\n'
                         'input=model.json\nkmesh=32,1,1\noutput=results\n$end\n')
            with redirect_stdout(io.StringIO()):
                self.assertEqual(main([str(p)]), 0)
                self.assertEqual(main([str(p), '--resume']), 0)
            result = json.loads((root/'results/result.json').read_text())
            self.assertAlmostEqual(result['finite_grid_integral'], 1., places=9)
            report = json.loads((root/'cube-results/workflow.json').read_text())
            self.assertTrue(report['tasks']['dos']['reused'])
            (root/'results/dos.csv').write_text('tampered')
            from wren.dos.runner import run
            with self.assertRaisesRegex(ValueError, 'missing or changed'):
                run(normalize({'mode': 'tight-binding', 'input': 'model.json', 'output': 'results', 'kmesh': '32,1,1'}, root), resume=True)

    def test_raw_network_schema_and_sum_rule(self):
        from wren.dos.runner import run, _network_terms
        network = Path(__file__).resolve().parents[1]/'examples/dos/nonorthogonal-chain.json'
        with tempfile.TemporaryDirectory() as temp:
            settings = normalize({'mode': 'tight-binding', 'input': network, 'output': temp,
                                  'kmesh': [9, 1, 1], 'site-projections': True}, Path(temp))
            result = run(settings)
            self.assertEqual(result['site_projections']['basis'], 'symmetric_lowdin_sites_within_model_cell')
            self.assertAlmostEqual(result['finite_grid_integral'], 2., places=10)
        model = json.loads(network.read_text())
        with self.assertRaisesRegex(ValueError, 'require schema'):
            _network_terms(dict(model, schema='wren-dos-tb-1', coupling_convention='signed_real_consistent_orbital_phases'))
        del model['overlap_pairs']
        with self.assertRaisesRegex(ValueError, 'overlap_pairs'):
            _network_terms(model)

    @unittest.skipUnless(importlib.util.find_spec('pyscf'), 'periodic PySCF required')
    def test_element_resolved_periodic_lih(self):
        import numpy as np
        from wren.dos.runner import run
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            source = root/'cell.json'
            source.write_text(json.dumps({'schema': 'cube-dimers-1', 'coordinate_unit': 'Angstrom',
                'lattice_angstrom': np.diag([5., 5., 5.]).tolist(), 'frame': 0,
                'monomers': [{'id': '1', 'elements': ['Li', 'H'], 'transport_position_angstrom': [2.5]*3,
                              'coordinates_angstrom': [[-.8, 0., 0.], [.8, 0., 0.]]}]}))
            settings = normalize({'mode': 'periodic', 'input': source, 'output': 'out', 'method': 'hf',
                                  'basis-set': 'gth-szv', 'pseudo': 'gth-pade', 'threads': 1,
                                  'kmesh': [3, 1, 1], 'precision': 1e-6, 'element-projections': True}, root)
            with redirect_stdout(io.StringIO()):
                result = run(settings)
            projected = result['element_projections']
            frame, = json.loads((root/'out/spectra.json').read_text())['frames']
            self.assertTrue(frame['scf_converged'])
            self.assertEqual(frame['scf']['modeled_electrons_per_cell'], 4)
            self.assertAlmostEqual(frame['scf']['occupation_integral_per_cell'], 4.)
            self.assertIsNone(frame['scf']['xc_grid'])
            self.assertLessEqual(frame['band_edges']['sampled_gap_ev'], frame['band_edges']['minimum_sampled_direct_gap_ev'])
            self.assertEqual(projected['elements'], ['H', 'Li'])
            for symbol, count in projected['expected_integral_per_element'].items():
                self.assertAlmostEqual(projected['finite_grid_integral_per_element'][symbol], count, places=8)
            self.assertAlmostEqual(projected['expected_integral_per_element']['H'], 2., places=8)
            self.assertLess(projected['max_abs_element_sum_residual_per_ev'], 1e-8)

    @unittest.skipUnless(importlib.util.find_spec('pyscf'), 'PySCF required')
    def test_real_monomer_multiframe_scf_and_raw_frontier_statistics(self):
        from wren.workflow.dimers import file_hash
        from wren.dos.runner import run
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            frames = []
            for i, length in enumerate((.74, .80)):
                here = root/f'frame{i:04d}'
                here.mkdir()
                p = here/'manifest.json'
                p.write_text(json.dumps({'schema': 'cube-dimers-1', 'coordinate_unit': 'Angstrom',
                    'frame': i, 'monomers': [{'id': '1', 'elements': ['H', 'H'],
                        'coordinates_angstrom': [[0., 0., 0.], [length, 0., 0.]]},
                        {'id': '2', 'elements': ['H', 'F'],
                         'coordinates_angstrom': [[0., 0., 0.], [length+.16, 0., 0.]]}], 'dimers': []}))
                frames.append({'frame': i, 'manifest': str(p.relative_to(root)), 'sha256': file_hash(p)})
            (root/'frames.json').write_text(json.dumps({'schema': 'cube-frames-1', 'frame_count': 2, 'frames': frames}))
            settings = normalize({'input': 'frames.json', 'output': 'results', 'method': 'hf', 'element-projections': True,
                                  'basis-set': 'sto-3g', 'threads': '1'}, root)
            with redirect_stdout(io.StringIO()):
                result = run(settings)
            self.assertEqual(result['frame_count'], 2)
            self.assertAlmostEqual(result['finite_grid_integral'], 8., places=10)
            integrals = result['element_projections']['finite_grid_integral_per_element']
            self.assertAlmostEqual(integrals['H'], 3., places=10)
            self.assertAlmostEqual(integrals['F'], 5., places=10)
            self.assertGreater(result['frontier_disorder']['homo']['variance_dynamic_ev2'], 0.)
            raw = json.loads((root/'results/spectra.json').read_text())
            self.assertNotEqual(raw['frames'][0]['molecules'][0]['homo_ev'], raw['frames'][1]['molecules'][0]['homo_ev'])
            with self.assertRaisesRegex(ValueError, 'potential alignment'):
                run(normalize({'mode': 'periodic', 'input': 'frames.json', 'output': 'periodic'}, root))


if __name__ == '__main__':
    unittest.main()
