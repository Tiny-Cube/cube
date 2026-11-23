from contextlib import redirect_stdout
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
            with self.assertRaisesRegex(ValueError, 'periodic GPU'):
                normalize({'mode': 'periodic', 'input': 'h2.xyz', 'qc-engine': 'gpu-pyscf'}, root)


@unittest.skipUnless(importlib.util.find_spec('numpy'), 'DOS NumPy dependency required')
class DosWorkflowTests(unittest.TestCase):
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
                        'coordinates_angstrom': [[0., 0., 0.], [length, 0., 0.]]}], 'dimers': []}))
                frames.append({'frame': i, 'manifest': str(p.relative_to(root)), 'sha256': file_hash(p)})
            (root/'frames.json').write_text(json.dumps({'schema': 'cube-frames-1', 'frame_count': 2, 'frames': frames}))
            settings = normalize({'input': 'frames.json', 'output': 'results', 'method': 'hf',
                                  'basis-set': 'sto-3g', 'threads': '1'}, root)
            with redirect_stdout(io.StringIO()):
                result = run(settings)
            self.assertEqual(result['frame_count'], 2)
            self.assertAlmostEqual(result['finite_grid_integral'], 4., places=10)
            self.assertGreater(result['frontier_disorder']['homo']['variance_dynamic_ev2'], 0.)
            raw = json.loads((root/'results/spectra.json').read_text())
            self.assertNotEqual(raw['frames'][0]['molecules'][0]['homo_ev'], raw['frames'][1]['molecules'][0]['homo_ev'])
            with self.assertRaisesRegex(ValueError, 'potential alignment'):
                run(normalize({'mode': 'periodic', 'input': 'frames.json', 'output': 'periodic'}, root))


if __name__ == '__main__':
    unittest.main()
