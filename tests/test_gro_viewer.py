"""Standalone geometry viewer: image selection, GRO identities and HTTP export."""
import importlib.util
import itertools
import json
from pathlib import Path
import sys
import tempfile
import threading
import unittest
from urllib.error import HTTPError
from urllib.request import urlopen

try:
    import numpy as np
    spec = importlib.util.spec_from_file_location('gro_viewer', Path(__file__).resolve().parents[1] / 'scripts/view_clusters.py')
    viewer = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = viewer
    spec.loader.exec_module(viewer)
except ImportError:
    np = None


def write_gro(path, records, box='1 1 1', precision=3, velocities=False):
    width = precision + 5
    rows = ['synthetic geometry', str(len(records))]
    for serial, (resid, name, xyz) in enumerate(records, 1):
        row = f'{resid:5d}{"MOL":<5}{name:>5}{serial:5d}'
        row += ''.join(f'{x / 10:{width}.{precision}f}' for x in xyz)
        if velocities:
            row += ''.join(f'{x:8.4f}' for x in (.1, -.2, .3))
        rows.append(row)
    path.write_text('\n'.join([*rows, box, '']))
    return path


@unittest.skipIf(np is None, 'NumPy required')
class GROViewerTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.path = Path(self.temp.name) / 'test.gro'

    def tearDown(self):
        self.temp.cleanup()

    def test_fixed_width_precision_velocity_and_residue_wrap(self):
        write_gro(self.path, [(99999, 'C', (1.2345, -2.3456, 3.4567)),
                              (0, 'H', (4, 5, 6))], precision=5, velocities=True,
                  box='1 2 3 0 0 .2 0 .3 .4')
        data = viewer.read_gro(self.path)
        np.testing.assert_allclose(data['atoms'][0]['coordinates'], [1.2345, -2.3456, 3.4567])
        np.testing.assert_allclose(data['cell'], [[10, 0, 0], [2, 20, 0], [3, 4, 30]])
        self.assertEqual(data['blocks'], [[0], [1]])
        self.assertEqual(data['atoms'][1]['residue'], 0)

    def test_whole_molecule_unwrap_and_strict_all_atom_boundary(self):
        # Carbon-carbon distance is 5 A, but hydrogen makes minimum distance 3 A.
        write_gro(self.path, [(1, 'C', (9, 0, 0)), (1, 'H', (0, 0, 0)),
                              (2, 'C', (4, 0, 0)), (2, 'H', (3, 0, 0))])
        structure = viewer.Structure(self.path, max_cutoff=6)
        np.testing.assert_allclose(structure.molecules[0].coordinates, [[9, 0, 0], [10, 0, 0]])
        at_boundary = structure.cluster(1, 3)
        self.assertEqual(at_boundary['neighbor_images'], 0)
        included = structure.cluster(1, 3.01)
        self.assertEqual(included['neighbor_images'], 1)
        self.assertEqual(included['neighbor_atoms'], 2)
        self.assertEqual(included['groups'][1]['image'], (1, 0, 0))
        self.assertEqual(included['global_indices'], [1, 2, 3, 4])
        self.assertEqual(structure.cluster(1, 0)['neighbor_atoms'], 0)
        self.assertGreater(structure.cluster(1, 6)['neighbor_images'], included['neighbor_images'])
        nonperiodic = viewer.Structure(self.path, pbc=False, max_cutoff=6)
        self.assertEqual(nonperiodic.cluster(1, 3.01)['neighbor_images'], 1)
        self.assertEqual(nonperiodic.cluster(1, 3.01)['groups'][1]['image'], (0, 0, 0))

    def test_triclinic_neighbors_match_independent_brute_force(self):
        write_gro(self.path, [(1, 'C', (1, 1, 1)), (1, 'H', (2, 1, 1)),
                              (2, 'N', (8, 7, 1)), (2, 'H', (9, 7, 1)),
                              (3, 'O', (4, 4, 7))], box='1 .9 .8 0 0 .4 0 .2 .3')
        structure = viewer.Structure(self.path, max_cutoff=12)
        target = structure.molecules[0].coordinates
        expected = set()
        for i, molecule in enumerate(structure.molecules, 1):
            for image in itertools.product(range(-3, 4), repeat=3):
                if i == 1 and image == (0, 0, 0):
                    continue
                moved = molecule.coordinates + np.array(image) @ structure.cell
                if np.linalg.norm(target[:, None] - moved[None, :], axis=-1).min() < 6:
                    expected.add((i, image))
        result = structure.cluster(1, 6)
        actual = {(g['molecule'], g['image']) for g in result['groups'][1:]}
        self.assertEqual(actual, expected)
        larger = {(g['molecule'], g['image']) for g in structure.cluster(1, 12)['groups'][1:]}
        self.assertTrue(any(i == 1 for i, image in larger))  # central periodic copies retained
        self.assertLess(actual, larger)

    def test_skew_cell_closest_vector_and_bond_graph(self):
        cell = np.array([[10., 0, 0], [9, 2, 0], [0, 0, 10]])
        delta = np.array([[4.8, 1.1, 0], [-4.8, -1.1, 0]])
        actual = viewer.minimum_images(delta, cell)
        candidates = np.array([delta + np.array(i) @ cell for i in itertools.product(range(-8, 9), repeat=3)])
        np.testing.assert_allclose(np.linalg.norm(actual, axis=1), np.linalg.norm(candidates, axis=2).min(axis=0))
        points = np.array([[9., 0, 0], [0, 0, 0], [1, 0, 0]])
        np.testing.assert_allclose(viewer.unwrap_molecule(points, np.eye(3) * 10, [[0, 1], [1, 2]]),
                                   [[9, 0, 0], [10, 0, 0], [11, 0, 0]])

    def test_invalid_frames_boxes_and_thresholds_rejected(self):
        write_gro(self.path, [(1, 'C', (0, 0, 0))], box='0 0 0')
        structure = viewer.Structure(self.path)
        self.assertFalse(structure.metadata()['periodic'])
        for center, cutoff in [(0, 6), (2, 6), (1, -1), (1, float('nan')), (1, 31)]:
            with self.assertRaises(ValueError):
                structure.cluster(center, cutoff)
        self.path.write_text(self.path.read_text() + 'another frame\n')
        with self.assertRaisesRegex(ValueError, 'single GRO frame'):
            viewer.read_gro(self.path)
        write_gro(self.path, [(1, 'C', (0, 0, 0))], box='1 0 0')
        with self.assertRaisesRegex(ValueError, 'nondegenerate'):
            viewer.read_gro(self.path)

    def test_http_cluster_xyz_and_invalid_queries(self):
        write_gro(self.path, [(1, 'C', (0, 0, 0)), (2, 'N', (4, 0, 0))])
        server = viewer.make_server(viewer.Structure(self.path, max_cutoff=6))
        worker = threading.Thread(target=server.serve_forever, daemon=True)
        worker.start()
        base = f'http://127.0.0.1:{server.server_port}'
        try:
            with urlopen(base + '/api/structure') as response:
                self.assertEqual(json.load(response)['molecule_count'], 2)
            with urlopen(base + '/api/cluster?center=1&cutoff=5') as response:
                self.assertEqual(json.load(response)['neighbor_atoms'], 1)
            with urlopen(base + '/api/cluster.xyz?center=1&cutoff=5') as response:
                lines = response.read().decode().splitlines()
                self.assertEqual(int(lines[0]), len(lines) - 2)
                self.assertIn('attachment;', response.headers['Content-Disposition'])
            for query in ['center=0', 'cutoff=nan', 'cutoff=7', 'center=1&center=2', 'other=1']:
                with self.assertRaises(HTTPError) as error:
                    urlopen(base + '/api/cluster?' + query)
                self.assertEqual(error.exception.code, 400)
        finally:
            server.shutdown()
            server.server_close()
            worker.join()


if __name__ == '__main__':
    unittest.main()
