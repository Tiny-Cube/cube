"""Join explicit region flags, stable atom identities, charged MM sites and resolved LJ."""
import copy
import importlib.util
import json
from pathlib import Path
import shutil
import tempfile
import unittest

try:
    import numpy as np
    from wren.workflow.separated_embedding import read_separated_embedding
except ImportError:
    np = None

EXAMPLE = Path(__file__).resolve().parents[1] / 'examples/reorganization/separated'


@unittest.skipIf(np is None, 'install QM/MM workflow dependencies')
class SeparatedEmbeddingTests(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name)
        self.xyz, self.top, self.ff = (self.root/name for name in ('coordinates.extxyz', 'topology.json', 'forcefield.json'))
        for p in (self.xyz, self.top, self.ff):
            shutil.copyfile(EXAMPLE/p.name, p)
        self.topology = json.loads(self.top.read_text())
        self.forcefield = json.loads(self.ff.read_text())

    def test_charged_mm_and_explicit_cross_parameters(self):
        atoms, env, molecule = read_separated_embedding(self.xyz, self.top, self.ff)
        self.assertEqual(molecule, '1')
        self.assertEqual([atom[0] for atom in atoms], ['O', 'H', 'H'])
        self.assertEqual(env.payload()['total_charge_e'], 1.)
        self.assertEqual(env.charges_e, (1.,))
        self.assertTrue(env.payload()['mm_frozen'])
        self.assertEqual(env.provenance['qm_atom_ids'], [1, 2, 3])
        self.assertEqual(env.provenance['mm_atom_ids'], [4])
        lj = env.vdw
        self.assertEqual(lj.c6_kj_mol_nm6[lj.qm_types[0]][lj.mm_types[0]], .0015)
        self.assertEqual(lj.c12_kj_mol_nm12[lj.qm_types[1]][lj.mm_types[0]], 0.)
        self.assertEqual(lj.provenance['type_names'], ['HW', 'Na', 'OW'])
        no_lj = read_separated_embedding(self.xyz, self.top)[1]
        self.assertIsNone(no_lj.vdw)
        self.assertNotEqual(env.fingerprint, no_lj.fingerprint)
        self.ff.write_text(self.ff.read_text() + '\n')
        self.assertNotEqual(read_separated_embedding(self.xyz, self.top, self.ff)[1].fingerprint, env.fingerprint)

    def test_join_by_id_with_reordered_columns_topology_and_extra_properties(self):
        self.topology['atoms'].reverse()
        self.top.write_text(json.dumps(self.topology))
        self.xyz.write_text('4\nProperties=region:S:1:atom_id:I:1:species:S:1:pos:R:3:frozen:L:1:v:R:3 '
                            'coordinate_unit="Angstrom"\nMM 4 Na 0 -2.44 0 True 0 0 0\n'
                            'QM 3 H -0.7572 0.5861 0 False 0 0 0\nQM 1 O 0 0 0 F 0 0 0\n'
                            'QM 2 H 0.7572 0.5861 0 F 0 0 0\n')
        atoms, env, _ = read_separated_embedding(self.xyz, self.top, self.ff)
        self.assertEqual([x[0] for x in atoms], ['H', 'O', 'H'])
        self.assertEqual(env.provenance['qm_atom_ids'], [3, 1, 2])
        names = env.vdw.provenance['type_names']
        self.assertEqual([names[i] for i in env.vdw.qm_types], ['HW', 'OW', 'HW'])

    def test_incomplete_cut_and_inconsistent_id_mapping(self):
        for field, value in [('atom_id', 2), ('atom_id', 7), ('element', 'Ca'), ('molecule_id', 1),
                             ('charge', True), ('image', [0, 0.5, 0])]:
            changed = copy.deepcopy(self.topology)
            changed['atoms'][-1][field] = value
            self.top.write_text(json.dumps(changed))
            with self.subTest(field=field, value=value), self.assertRaises(ValueError):
                read_separated_embedding(self.xyz, self.top)
        for bonds in ([[1, 4]], [[1, 1]], [[1, 2], [2, 1]], [[1, 5]]):
            self.top.write_text(json.dumps({**self.topology, 'bonds': bonds}))
            with self.subTest(bonds=bonds), self.assertRaises(ValueError):
                read_separated_embedding(self.xyz, self.top)

    def test_region_unit_and_geometry_validation(self):
        original = self.xyz.read_text()
        for a, b in [('QM F', 'QM T'), ('MM T', 'MM F'), ('MM T', 'MM X'), ('Angstrom', 'nm'),
                     ('region:S:1', 'zone:S:1'), ('atom_id:I:1', 'atom_id:R:1'),
                     ('MM T 4', 'MM T 1'), ('Angstrom', 'Angstrom pbc="T T T"'),
                     ('-2.4400', '-0.1'), ('-2.4400', '0.0000')]:
            self.xyz.write_text(original.replace(a, b))
            with self.subTest(change=b), self.assertRaises(ValueError):
                read_separated_embedding(self.xyz, self.top)

    def test_periodic_copy_of_central_molecule_is_a_separate_complete_mm_image(self):
        self.xyz.write_text('6\nProperties=species:S:1:pos:R:3:region:S:1:frozen:L:1:atom_id:I:1 coordinate_unit=Angstrom\n'
                            'O 0 0 0 QM F 1\nH .75 .58 0 QM F 2\nH -.75 .58 0 QM F 3\n'
                            'O 5 0 0 MM T 4\nH 5.75 .58 0 MM T 5\nH 4.25 .58 0 MM T 6\n')
        base = self.topology['atoms'][:3]
        image = [{**a, 'atom_id': i+4, 'source_atom_id': i+1, 'image': [1, 0, 0]} for i, a in enumerate(base)]
        self.top.write_text(json.dumps({**self.topology, 'atoms': base + image,
                                        'bonds': [[1, 2], [1, 3], [4, 5], [4, 6]]}))
        _, env, central = read_separated_embedding(self.xyz, self.top)
        self.assertEqual(central, '1')
        self.assertEqual(env.charges_e, (-.834, .417, .417))
        self.assertEqual(env.provenance['mm_atom_ids'], [4, 5, 6])

    def test_no_inferred_lj_or_unsupported_units(self):
        pairs = self.forcefield['pairs']
        invalid = [{'pairs': pairs[:1]}, {'pairs': [*pairs, {**pairs[0], 'types': list(reversed(pairs[0]['types']))}]},
                   {'units': {'energy': 'eV', 'length': 'Angstrom'}}, {'model': 'LJ9-6'},
                   {'pairs': [{**pairs[0], 'c6': True}, pairs[1]]}]
        for update in invalid:
            self.ff.write_text(json.dumps({**self.forcefield, **update}))
            with self.subTest(update=update), self.assertRaises(ValueError):
                read_separated_embedding(self.xyz, self.top, self.ff)

    @unittest.skipIf(importlib.util.find_spec('pyscf') is None, 'PySCF required for actual QM/MM energy')
    def test_physical_single_point_matches_manual_charge_and_lj(self):
        from scipy.constants import Avogadro, physical_constants
        from pyscf import gto, lib, qmmm, scf
        from wren.quantum.config import SCFConfig
        from wren.workflow.scf import make_mean_field
        atoms, env, _ = read_separated_embedding(self.xyz, self.top, self.ff)
        with lib.with_omp_threads(1):
            conf = SCFConfig(xc='hf', basis='sto-3g', threads=1, auxbasis='weigend')
            mf = make_mean_field(atoms, conf, embedding=env)
            mf.kernel()
            self.assertTrue(mf.converged)
            direct = gto.M(atom=atoms, basis='sto-3g', unit='Angstrom', verbose=0)
            reference = qmmm.mm_charge(scf.RHF(direct).density_fit(auxbasis='weigend'), [[0., -2.44, 0.]], [1.], unit='Angstrom')
            reference.conv_tol = conf.conv_tol
            reference.kernel()
            self.assertTrue(reference.converged)
            kj_per_hartree = physical_constants['Hartree energy'][0] * Avogadro / 1000
            analytic_lj = (1e-6/.244**12 - .0015/.244**6) / kj_per_hartree
            energy, _ = env.classical(atoms)
            self.assertAlmostEqual(energy, analytic_lj, places=13)
            self.assertAlmostEqual(mf.e_tot + energy, reference.e_tot + analytic_lj, places=9)
