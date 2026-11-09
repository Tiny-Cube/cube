"""Frozen-MM physical consistency, source integrity and checkpoint isolation."""
import json
import importlib.util
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from universal_kmc.workflow.control import load_input

try:
    import numpy as np
    from universal_kmc.workflow.embedding import FrozenMM, cluster_from_frame, read_charges
except ImportError:
    np = None

try:
    import pyscf
except ImportError:
    pyscf = None


class FrozenMMControlTests(unittest.TestCase):
    def test_explicit_file_and_manifest_controls_are_unambiguous(self):
        with tempfile.TemporaryDirectory() as directory:
            p = Path(directory)
            (p/'qm.xyz').write_text('2\nHF\nH 0 0 0\nF 1 0 0\n')
            (p/'env.pc').write_text('0 3 0 .1\n')
            (p/'manifest.json').write_text('{}\n')
            source = '$reorganization-energy\nxyz = qm.xyz\nenvironment = qmmm\nmm-charges = env.pc\n$end\n'
            (p/'input').write_text(source)
            s = load_input(p/'input').sections['reorganization-energy']
            self.assertTrue(s['mm-frozen'])
            self.assertEqual(s['mm-source'], 'file')
            self.assertEqual(s['mm-charges'], str(p/'env.pc'))
            manifest = source.replace('xyz = qm.xyz', 'input = manifest.json').replace('mm-charges = env.pc', 'mm-source = manifest\nmm-cutoff = 9')
            (p/'input').write_text(manifest)
            self.assertEqual(load_input(p/'input').sections['reorganization-energy']['mm-cutoff'], 9.)
            for invalid in (source.replace('environment = qmmm', 'environment = vacuum'),
                            source.replace('$end', 'mm-frozen = false\n$end'),
                            source.replace('$end', 'mm-cutoff = 9\n$end'),
                            source.replace('mm-charges = env.pc', 'mm-source = manifest'),
                            source.replace('xyz = qm.xyz', 'input = manifest.json'),
                            manifest.replace('mm-cutoff = 9', 'mm-cutoff = nan'),
                            manifest.replace('$end', 'mm-charges = env.pc\n$end')):
                (p/'input').write_text(invalid)
                with self.assertRaises(ValueError):
                    load_input(p/'input')


@unittest.skipIf(np is None, 'NumPy dependency required')
class FrozenMMGeometryTests(unittest.TestCase):
    def test_point_charge_format_finite_data_and_collision(self):
        with tempfile.TemporaryDirectory() as directory:
            p = Path(directory)/'environment.pc'
            p.write_text('# x y z q; Angstrom/e\n1 2 3 .2 # site\n-1 2 3 -.2\n')
            env = read_charges(p)
            self.assertEqual(env.charges_e, (.2, -.2))
            self.assertEqual(env.payload()['total_charge_e'], 0.)
            with self.assertRaisesRegex(ValueError, '0.5 Angstrom'):
                env.validate_qm([('H', (1., 2., 3.))])
            for content in ('', '0 0 nan 1\n', '0 1 2\n', 'q x y z\n'):
                p.write_text(content)
                with self.assertRaises(ValueError):
                    read_charges(p)

    def test_whole_molecule_images_match_brute_force_and_translate_together(self):
        import itertools
        cell = np.array([[8., 0, 0], [2., 7., 0], [1., 1., 9.]])
        fragments = [
            {'indices': [0, 1], 'elements': ['H', 'F'], 'coordinates': [[.1, .2, .3], [.9, .2, .3]], 'charges': [.2, -.2]},
            {'indices': [2, 3], 'elements': ['H', 'F'], 'coordinates': [[7.7, 6.7, 8.7], [8.5, 6.7, 8.7]], 'charges': [.2, -.2]}]
        frame = {'lattice': cell.tolist(), 'fragments': fragments, 'provenance': {'fixture': True}}
        qm = np.array(fragments[0]['coordinates'])
        record = {'elements': ['H', 'F'], 'global_atom_indices_zero_based': [0, 1], 'coordinates_angstrom': qm.tolist()}
        env = cluster_from_frame(frame, record, 9.)
        expected = set()
        for f in fragments:
            for image in itertools.product(range(-3, 4), repeat=3):
                if f is fragments[0] and image == (0, 0, 0):
                    continue
                xyz = np.array(f['coordinates'])+np.array(image)@cell
                if np.linalg.norm(qm[:, None]-xyz[None, :], axis=2).min() < 9:
                    expected.add((tuple(f['indices']), image))
        actual = {(tuple(m['global_atom_indices_zero_based']), tuple(m['image'])) for m in env.provenance['molecular_images']}
        self.assertEqual(actual, expected)
        self.assertNotIn(((0, 1), (0, 0, 0)), actual)
        self.assertEqual(len(env.charges_e), 2*len(expected))
        shift = np.array([5., -4., 3.])
        moved = cluster_from_frame(frame, {**record, 'coordinates_angstrom': (qm+shift).tolist()}, 9.)
        np.testing.assert_allclose(moved.coordinates_angstrom, np.array(env.coordinates_angstrom)+shift, atol=1e-12)
        self.assertEqual(moved.charges_e, env.charges_e)
        with self.assertRaisesRegex(ValueError, 'complete topology fragment'):
            cluster_from_frame(frame, {**record, 'global_atom_indices_zero_based': [0]}, 9.)
        with self.assertRaisesRegex(ValueError, 'common translation'):
            cluster_from_frame(frame, {**record, 'coordinates_angstrom': (qm*2).tolist()}, 9.)

    @unittest.skipIf(importlib.util.find_spec('MDAnalysis') is None, 'MDAnalysis dependency required')
    def test_gromacs_source_hash_change_is_rejected_before_loading(self):
        from universal_kmc.workflow.embedding import load_gromacs_environment
        with tempfile.TemporaryDirectory() as directory:
            p = Path(directory)/'a.tpr'
            p.write_bytes(b'changed')
            with self.assertRaisesRegex(ValueError, 'tpr changed'):
                load_gromacs_environment({'provenance': {'settings': {'tpr': str(p)}, 'sha256': {'tpr': 'invalid'}}})


@unittest.skipIf(pyscf is None or np is None, 'quantum dependencies required')
class FrozenMMPhysicsTests(unittest.TestCase):
    def settings(self, directory):
        p = Path(directory)
        (p/'HF.xyz').write_text('2\nHF\nH -.4 0 0\nF .5 0 0\n')
        (p/'mm.pc').write_text('2 3 0 .2\n-1 -3 0 -.1\n')
        (p/'input').write_text('''$reorganization-energy
xyz = HF.xyz
environment = qmmm
mm-charges = mm.pc
qc-engine = pyscf
method = hf
basis-set = sto-3g
conv-tol = 1e-11
$end
''')
        return load_input(p/'input').sections['reorganization-energy']

    def test_energy_embedding_gradients_frozen_sites_and_translation_invariance(self):
        from pyscf.data.nist import BOHR
        from pyscf.qmmm.itrf import QMMMGrad
        from universal_kmc.workflow.control import scf_config
        from universal_kmc.workflow.scf import make_mean_field, scf_record
        with tempfile.TemporaryDirectory() as directory:
            settings = self.settings(directory)
            config = scf_config(settings)
            env = read_charges(settings['mm-charges'])
            atoms = [('H', (-.4, 0., 0.)), ('F', (.5, 0., 0.))]
            vacuum = make_mean_field(atoms, config)
            vacuum.kernel()
            for charge in (0, 1, -1):
                energies = []
                for offset in (0., 1e-3, -1e-3):
                    moved = [atoms[0], ('F', (.5+offset*BOHR, 0., 0.))]
                    mf = make_mean_field(moved, config, charge=charge, spin=int(charge != 0), embedding=env)
                    mf.kernel()
                    scf_record(mf, 'cpu')
                    env.assert_attached(mf)
                    energies.append(float(mf.e_tot))
                    if offset == 0:
                        gradient = mf.nuc_grad_method()
                        self.assertIsInstance(gradient, QMMMGrad)
                        analytic = gradient.kernel()[1, 0]
                self.assertLess(abs(analytic-(energies[1]-energies[2])/.002), 2e-5)
                if charge == 0:
                    self.assertGreater(abs(energies[0]-vacuum.e_tot), 1e-5)
                    shift = np.array([4., -2., 1.])
                    translated = FrozenMM(np.array(env.coordinates_angstrom)+shift, env.charges_e, {})
                    translated_mf = make_mean_field([(e, np.array(r)+shift) for e, r in atoms], config, embedding=translated)
                    translated_mf.kernel()
                    self.assertAlmostEqual(translated_mf.e_tot, energies[0], places=9)
                    translated_mf.mm_mol.set_geom_(np.array(translated.coordinates_angstrom)+[.01, 0, 0], unit='Angstrom')
                    with self.assertRaisesRegex(RuntimeError, 'frozen MM'):
                        translated.assert_attached(translated_mf)

    def test_optimizer_initialization_failure_is_recorded(self):
        from universal_kmc.workflow.control import scf_config
        from universal_kmc.workflow.reorganization import _optimize
        with tempfile.TemporaryDirectory() as directory:
            settings = self.settings(directory)
            root = Path(directory)
            with patch('universal_kmc.workflow.reorganization.make_mean_field', side_effect=RuntimeError('initialization failed')):
                with self.assertRaisesRegex(RuntimeError, 'initialization failed'):
                    _optimize([('H', (-.4, 0., 0.)), ('F', (.5, 0., 0.))], scf_config(settings), 0, root, 10,
                              embedding=read_charges(settings['mm-charges']))
            self.assertEqual(json.loads((root/'neutral.optimization.json').read_text())['status'], 'failed')

    def test_embedded_checkpoint_isolated_from_vacuum_and_mm_edits_invalidate_resume(self):
        from pyscf.scf import chkfile
        from universal_kmc.workflow.scf import run_single_point, translated_checkpoint
        from universal_kmc.quantum.checkpoints import _load_checkpoint
        from universal_kmc.workflow.reorganization import calculate_reorganization
        with tempfile.TemporaryDirectory() as directory:
            settings = self.settings(directory)
            root = Path(directory)
            atoms = [('H', (-.4, 0., 0.)), ('F', (.5, 0., 0.))]
            env = read_charges(settings['mm-charges'])
            result = run_single_point(atoms, settings, root/'point', embedding=env)
            self.assertEqual(result['embedding_fingerprint'], env.fingerprint)
            self.assertIsNone(chkfile.load(result['checkpoint'], 'cube_transfer'))
            marker = chkfile.load(result['checkpoint'], 'cube_reorganization')
            np.testing.assert_array_equal(marker['mm_charges_e'], env.charges_e)
            for assumed in (False, True):
                with self.assertRaisesRegex(ValueError, 'QM/MM reorganization'):
                    _load_checkpoint(result['checkpoint'], assumed)
            with self.assertRaises(ValueError):
                translated_checkpoint(result['checkpoint'], atoms, root/'translated.chk')
            self.assertTrue(run_single_point(atoms, settings, root/'point', embedding=env, resume=True)['reused'])
            changed = FrozenMM(env.coordinates_angstrom, [.3, -.1], {})
            with self.assertRaisesRegex(ValueError, 'identical settings'):
                run_single_point(atoms, settings, root/'point', embedding=changed, resume=True)
            with patch('universal_kmc.workflow.reorganization._optimize', side_effect=RuntimeError('forced failed optimization')):
                with self.assertRaises(RuntimeError):
                    calculate_reorganization(atoms, settings, root/'reorg', embedding=env)
            snapshot = json.loads((root/'reorg/mm-environment.json').read_text())
            self.assertEqual(snapshot['fingerprint'], env.fingerprint)
            with self.assertRaisesRegex(ValueError, 'identical inputs'):
                calculate_reorganization(atoms, settings, root/'reorg', embedding=changed, resume=True)


if __name__ == '__main__':
    unittest.main()
