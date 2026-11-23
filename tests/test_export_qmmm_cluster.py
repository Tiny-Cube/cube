"""Exercise real export/selection/serialization with a small in-memory topology."""
import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch


@unittest.skipIf(importlib.util.find_spec('scipy') is None, 'export needs molecular workflow dependencies')
class ClusterExportTests(unittest.TestCase):
    def setUp(self):
        source = Path(__file__).resolve().parents[1]/'scripts/export_qmmm_cluster.py'
        spec = importlib.util.spec_from_file_location('wren_export_test', source)
        self.exporter = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(self.exporter)
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        root = Path(tmp.name)
        self.dest, self.path = root/'selected', root/'frame.json'
        tpr = root/'synthetic.tpr'
        tpr.write_text('synthetic topology fixture')
        fingerprint = hashlib.sha256(tpr.read_bytes()).hexdigest()
        positions = [[0., 0., 0.], [1., 0., 0.]]
        self.frame = {'lattice': [[18.,0.,0.],[0.,18.,0.],[0.,0.,18.]], 'fragments': [
            {'indices':[0,1], 'elements':['H','F'], 'coordinates':positions, 'charges':[.2,-.2]},
            {'indices':[2], 'elements':['Na'], 'coordinates':[[4.,1.,0.]], 'charges':[1.]}
        ], 'provenance': {'sha256': {'tpr': fingerprint}}}
        self.manifest = {'schema':'cube-dimers-1', 'coordinate_unit':'Angstrom', 'frame':0,
                         'provenance':{'settings':{'tpr':str(tpr)}},
                         'monomers':[{'id':'1', 'global_atom_indices_zero_based':[0,1], 'elements':['H','F'],
                                      'coordinates_angstrom': positions, 'cluster':{'cutoff_angstrom':4.}},
                                     {'id':'2', 'global_atom_indices_zero_based':[2]}]}
        self.path.write_text(json.dumps(self.manifest))
        self.top = (['H','F','Na'], ['h','f','na'], [.2,-.2,1.], [[1],[0],[]])

    def test_explicit_bonds_original_ids_cutoff_and_overwrite_refusal(self):
        from wren.workflow.separated_embedding import read_separated_embedding
        with patch('wren.workflow.embedding.load_gromacs_environment', return_value=self.frame), \
             patch.object(self.exporter, '_read_topology', return_value=self.top):
            result = self.exporter.export_cluster(self.path, 1, self.dest)
        self.assertEqual(result['cutoff_angstrom'], 4.)
        atoms, mm, identity = read_separated_embedding(self.dest/'coordinates.extxyz', self.dest/'topology.json')
        self.assertEqual([x[0] for x in atoms], ['H','F'])
        self.assertEqual((mm.charges_e, identity), ((1.,), '1'))
        data = json.loads((self.dest/'topology.json').read_text())
        self.assertEqual(data['bonds'], [[1,2]])
        self.assertEqual([a['source_atom_id'] for a in data['atoms']], [1,2,3])
        self.assertFalse((self.dest/'forcefield.json').exists())
        marker = self.dest/'my-previous-result'
        marker.write_text('preserve')
        with self.assertRaisesRegex(FileExistsError, 'new directory'):
            self.exporter.export_cluster(self.path, 1, self.dest)
        self.assertEqual(marker.read_text(), 'preserve')

    def test_bond_cut_and_empty_shell_publish_nothing(self):
        with patch('wren.workflow.embedding.load_gromacs_environment', return_value=self.frame), \
             patch.object(self.exporter, '_read_topology', return_value=(*self.top[:3], [[1],[0,2],[1]])):
            with self.assertRaisesRegex(ValueError, 'cuts a TPR bond'):
                self.exporter.export_cluster(self.path, 1, self.dest)
            with self.assertRaisesRegex(ValueError, 'nonempty'):
                self.exporter.export_cluster(self.path, 1, self.dest, cutoff=.1)
            self.assertFalse(self.dest.exists())
