"""Round-trip the production extractor on a triclinic, two-frame LAMMPS trajectory."""
import importlib.util
import itertools
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from wren.workflow.control import load_input
from wren.workflow.frame_geometry import indexed_xyz, load_frame
from wren.workflow.runner import _dataset_sources

EXAMPLE = Path(__file__).resolve().parents[1] / 'examples/lammps/input'


@unittest.skipUnless(all(importlib.util.find_spec(m) for m in ('h5py', 'MDAnalysis')), 'workflow extra required')
class ArchiveTests(unittest.TestCase):
    def setUp(self):
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.output = Path(directory.name)
        self.settings = dict(load_input(EXAMPLE).sections['get_dimers_monomers'], output=str(self.output))
        self.settings.pop('format')  # exercise the default public entry point

    def test_indices_and_distances_match_exhaustive_periodic_search(self):
        import h5py
        import numpy as np
        from wren.workflow.dimers import extract_dimers_monomers
        self.settings['monomers-qmmm-cutoff'] = 11.
        path = extract_dimers_monomers(self.settings)
        self.assertEqual(path, self.output/'geometry.wren.h5')
        sources = _dataset_sources({'input': str(path)}, None)
        self.assertEqual([load_frame(s)['lammps_timestep'] for s in sources], [0,200])
        with self.assertRaisesRegex(ValueError, 'choose a frame'):
            load_frame(path)
        with h5py.File(path) as h:
            self.assertEqual(h.attrs['schema'], 'wren-geometry/1')
            self.assertEqual(h['frames/id'].id, h['neighbors/frame_id'].id)
            np.testing.assert_array_equal(h['topology/molecule_id'][:], [7,19])
            for frame in (0,1):
                data = load_frame(path, frame=frame)
                self.assertEqual(data['dimers'][0]['id'], '7_19_img_0_-1_0')
                xyz = h['geometry/positions'][frame].reshape(2,2,3)
                cell = h['frames/cell'][frame]
                lines = indexed_xyz(self.output, f'frame{frame}.mol19').splitlines()[2:]
                np.testing.assert_array_equal(xyz[1], [list(map(float,s.split()[1:])) for s in lines])
                for kind in ('monomer','dimer'):
                    g = h['neighbors/'+kind]
                    a,b = g['frame_ptr'][frame:frame+2]
                    self.assertEqual(b-a, 2 if kind=='monomer' else 1)
                    for row in range(a,b):
                        ids = np.atleast_1d(g['qm_mol_id'][row]).tolist()
                        qm = [(ids[0], (0,0,0))]
                        if kind=='dimer':
                            qm.append((ids[1],tuple(g['qm_b_image'][row])))
                        coords = np.concatenate([xyz[[7,19].index(i)] + np.array(img)@cell for i,img in qm])
                        expected = {}
                        for k, identity in enumerate((7,19)):
                            for shift in itertools.product(range(-3,4), repeat=3):
                                if (identity,shift) in qm:
                                    continue
                                delta = coords[:,None,:] - (xyz[k]+np.array(shift)@cell)
                                d = np.linalg.norm(delta,axis=2).min()
                                if d < 11:
                                    expected[identity,shift] = d
                        lo,hi = g['mm_ptr'][row:row+2]
                        actual = {(i,tuple(img)):dist for i,img,dist in zip(g['mm_mol_id'][lo:hi],
                                  g['mm_image'][lo:hi],g['mm_min_angstrom'][lo:hi], strict=True)}
                        self.assertEqual(len(actual),hi-lo)
                        self.assertEqual(actual.keys(),expected.keys())
                        for site in expected:
                            self.assertAlmostEqual(actual[site],expected[site])
        self.assertEqual(extract_dimers_monomers(self.settings, resume=True), path)
        with h5py.File(path,'r+') as h:
            h['neighbors/dimer/mm_mol_id'][0] = 123
        with self.assertRaisesRegex(ValueError,'changed'):
            extract_dimers_monomers(self.settings, resume=True)

    def test_empty_tables_and_atomic_failure(self):
        import h5py
        from wren.workflow.dimers import extract_dimers_monomers
        self.settings.update(dimers_cutoff=.01, **{'monomers-qmmm-cutoff':.01})
        with patch('wren.workflow.hdf_archive._write_frame', side_effect=ValueError('interrupted')):
            with self.assertRaisesRegex(ValueError,'interrupted'):
                extract_dimers_monomers(self.settings)
        self.assertEqual(list(self.output.iterdir()), [])
        path = extract_dimers_monomers(self.settings)
        with h5py.File(path) as h:
            self.assertEqual(h['neighbors/monomer/frame_ptr'][:].tolist(), [0,2,4])
            self.assertEqual(h['neighbors/dimer/frame_ptr'][:].tolist(), [0,0,0])
            self.assertEqual(h['neighbors/monomer/mm_ptr'][:].tolist(), [0]*5)
            self.assertEqual(h['neighbors/dimer/mm_ptr'][:].tolist(), [0])
        with self.assertRaisesRegex(ValueError,'output exists'):
            extract_dimers_monomers(self.settings)
