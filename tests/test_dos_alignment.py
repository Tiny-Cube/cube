from contextlib import redirect_stdout
import csv
import importlib.util
import io
import json
from pathlib import Path
import tempfile
import unittest

from wren.dos.control import normalize
from wren.workflow.dimers import file_hash


@unittest.skipUnless(importlib.util.find_spec('numpy'), 'NumPy required')
class AlignmentTests(unittest.TestCase):
    def test_vbm_alignment_preserves_gap_and_rejects_partial_filling(self):
        import numpy as np
        from wren.dos.alignment import aligned_levels
        edges = {'vbm_ev': -2., 'cbm_ev': 1., 'sampled_gap_ev': 3., 'partial_occupation': False}
        frame = {'band_edges': edges}
        settings = {'periodic-reference': 'frame-vbm'}
        np.testing.assert_array_equal(aligned_levels(settings, frame, np.array([-2., 1.]), {}), [0., 3.])
        self.assertEqual(edges['vbm_ev'], -2.)
        self.assertEqual(frame['energy_alignment']['aligned_cbm_ev'], 3.)
        edges['partial_occupation'] = True
        with self.assertRaisesRegex(ValueError, 'sampled gap'):
            aligned_levels(settings, frame, np.array([-2., 1.]), {})

    def test_external_reference_gauge_invariance_and_hash_binding(self):
        import numpy as np
        from wren.dos.alignment import aligned_levels, reference_values
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root/'frames.json'
            source.write_text('{}')
            table = root/'reference.csv'
            table.write_text('frame,reference_ev,manifest_sha256\n17,-1.75,abc\n3,-2,def\n')
            settings = normalize({'mode': 'periodic', 'input': source, 'reference-file': table}, root)
            datasets = [{'monomers': [], 'frame': index, '_source_sha256': digest}
                        for index, digest in [(3, 'def'), (17, 'abc')]]
            refs = reference_values(settings, datasets)
            energies = np.array([-4., 1.])
            frame = {'frame': 3, 'band_edges': {'vbm_ev': -4., 'cbm_ev': 1.}}
            a = aligned_levels(settings, frame, energies, refs)
            shifted = {'frame': 17, 'band_edges': {'vbm_ev': -3.75, 'cbm_ev': 1.25}}
            np.testing.assert_array_equal(a, aligned_levels(settings, shifted, energies+.25, refs))
            self.assertEqual(shifted['energy_alignment']['aligned_cbm_ev'], 3.)
            for content in ('17,-1,abc\n', '17,nan,abc\n3,-2,def\n', '17,-1,abc\n17,2,abc\n',
                            '17,-1,bad\n3,-2,def\n', '17,-1,abc,wrong\n3,-2,def\n'):
                table.write_text('frame,reference_ev,manifest_sha256\n'+content)
                with self.assertRaises(ValueError):
                    reference_values(settings, datasets)
            changed = [{**datasets[0], 'monomers': [{'id': '2', 'elements': ['He']}]}, datasets[1]]
            with self.assertRaisesRegex(ValueError, 'species'):
                reference_values(settings, changed)

    @unittest.skipUnless(importlib.util.find_spec('pyscf'), 'periodic PySCF required')
    def test_gth_cell_mean_matches_independent_radial_integral(self):
        import math
        import numpy as np
        from scipy.integrate import quad
        from pyscf.pbc import gto
        from pyscf.data.nist import HARTREE2EV
        from wren.dos.alignment import mean_local_potential_ev
        for edge in (6., 7.):
            cell = gto.Cell(atom='C 0 0 0; H 2 0 0; H 0 2 0', a=np.eye(3)*edge,
                            basis='gth-szv', pseudo='gth-pbe', unit='Angstrom', verbose=0).build()
            # Integrate V_local(r) + Z/r for each isolated ion. Subtracting the
            # point Coulomb kernel isolates the finite neutral-cell G=0 coefficient.
            total = 0.
            for i in range(cell.natm):
                z = cell.atom_charge(i)
                radius, _, coefficients = cell._pseudo[cell.atom_symbol(i)][1:4]
                def integrand(r):
                    u = r/radius
                    short = math.exp(-.5*u*u)*sum(c*u**(2*j) for j, c in enumerate(coefficients))
                    return 4*math.pi*(z*r*math.erfc(u/math.sqrt(2))+r*r*short)
                total += quad(integrand, 0., 20*radius)[0]
            self.assertAlmostEqual(mean_local_potential_ev(cell), total*HARTREE2EV/cell.vol, places=11)

    @unittest.skipUnless(importlib.util.find_spec('pyscf'), 'periodic PySCF required')
    def test_multiframe_periodic_dos_pdos_blocks_and_raw_spectrum(self):
        import numpy as np
        from wren.dos.runner import run
        template = Path(__file__).resolve().parents[1]/'examples/dos/h2-cell.json'
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            entries = []
            for frame, length in [(10, .74), (20, .84)]:
                cell = json.loads(template.read_text())
                cell.update(frame=frame, time_ps=frame/10)
                cell['monomers'][0]['coordinates_angstrom'] = [[-length/2, 0, 0], [length/2, 0, 0]]
                path = root/f'cell{frame}.json'
                path.write_text(json.dumps(cell))
                entries.append({'frame': frame, 'manifest': path.name, 'sha256': file_hash(path)})
            (root/'frames.json').write_text(json.dumps({'schema': 'cube-frames-1', 'frame_count': 2, 'frames': entries}))
            control = {'mode': 'periodic', 'input': 'frames.json', 'output': 'bands', 'method': 'hf',
                       'pseudo': 'gth-pade', 'threads': 1, 'precision': 1e-6, 'element-projections': True,
                       'block-size-frames': 1, 'periodic-reference': 'cell-mean-local', 'sigma': .2}
            with redirect_stdout(io.StringIO()):
                result = run(normalize(control, root))
            self.assertEqual(result['energy_reference'], 'per_cell_mean_hartree_plus_local_gth')
            self.assertTrue(result['sampling_uncertainty_estimated'])
            self.assertEqual(result['sampling']['full_blocks'], 2)
            frames = json.loads((root/'bands/spectra.json').read_text())['frames']
            refs = [f['cell_mean_local_ev'] for f in frames]
            self.assertEqual(refs[0], refs[1])
            self.assertNotEqual(frames[0]['band_edges']['vbm_ev'], frames[1]['band_edges']['vbm_ev'])
            dos = np.genfromtxt(root/'bands/dos.csv', delimiter=',', names=True)
            pdos = np.genfromtxt(root/'bands/pdos.csv', delimiter=',', names=True)
            se = np.genfromtxt(root/'bands/dos-blocks.csv', delimiter=',', names=True)
            raw = np.array([f['orbital_energies_ev'][0] for f in frames])
            x = dos['energy_ev']
            sigma = control['sigma']
            direct = 2*np.exp(-.5*((x[:, None, None]-(raw-refs[0]))/sigma)**2).sum(axis=2).mean(axis=1)
            np.testing.assert_allclose(dos['dos'], direct/(np.sqrt(2*np.pi)*sigma), atol=1e-12)
            np.testing.assert_allclose(dos['dos'], pdos['H'], atol=1e-12)
            np.testing.assert_allclose(se['standard_error'], abs(dos['frame_0']-dos['frame_1'])/2, atol=1e-12)
            self.assertAlmostEqual(result['finite_grid_integral'], 4., places=9)
            with (root/'reference.csv').open('w') as stream:
                writer = csv.writer(stream)
                writer.writerow(['frame', 'reference_ev'])
                writer.writerows((f['frame'], ref) for f, ref in zip(frames, refs))
            external = normalize(dict(control, output='external', **{'periodic-reference': 'external',
                                                                      'reference-file': 'reference.csv'}), root)
            with redirect_stdout(io.StringIO()):
                repeat = run(external)
            self.assertEqual(result['output_sha256']['dos.csv'], repeat['output_sha256']['dos.csv'])
            (root/'reference.csv').write_text('frame,reference_ev\n10,1\n20,2\n')
            with self.assertRaisesRegex(ValueError, 'identical'):
                run(external, resume=True)
