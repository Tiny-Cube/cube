"""Finite-cluster material Hamiltonian, forces, invariances and FSSH end-to-end."""
import copy
from contextlib import redirect_stdout
import io
import json
from pathlib import Path
import tempfile
import unittest

from wren.surface_hopping.control import load_control
from wren.surface_hopping.fragment_sites import read_fragment_model

try:
    import numpy as np
except ImportError:
    np = None
if np is not None:
    from wren.surface_hopping.engine import adiabatic_frame
    from wren.surface_hopping.models import FragmentSitesModel
    from wren.surface_hopping.runner import run_control

EXAMPLE = Path(__file__).resolve().parents[1] / 'examples/surface-hopping/fragment-sites.toml'


class FragmentControlTests(unittest.TestCase):
    def test_reference_masses_site_and_geometry(self):
        config, _, source, _ = load_control(EXAMPLE)
        params, digest = read_fragment_model(EXAMPLE.with_suffix('.json'))
        self.assertEqual(config.masses, params.masses)
        self.assertEqual(config.coordinates, params.coordinates)
        self.assertEqual(config.fragment_site, 'A')
        self.assertEqual(config.masses[3:6], (310000.,) * 3)
        self.assertIn(digest, source['sources'].values())

    def test_invalid_parameter_combinations(self):
        model = json.loads(EXAMPLE.with_suffix('.json').read_text())
        bad_pairs = [
            {**model['pairs'][0], 'sites': ['A', 'A']}, {**model['pairs'][0], 'sites': ['A', 'X']},
            {**model['pairs'][0], 'sites': ['B', 'A']},
        ]
        with tempfile.TemporaryDirectory() as root:
            path = Path(root) / 'invalid.json'
            for extra_pair in bad_pairs:
                with self.subTest(pair=extra_pair), self.assertRaises(ValueError):
                    path.write_text(json.dumps({**model, 'pairs': [*model['pairs'], extra_pair]}))
                    read_fragment_model(path)
            for replacement in [{'coupling_convention': 'absolute'}, {'basis': 'moving'},
                                {'pairs': [{**model['pairs'][0], 'decay_per_bohr': -1}]}]:
                with self.subTest(replacement=replacement), self.assertRaises(ValueError):
                    path.write_text(json.dumps({**model, **replacement}))
                    read_fragment_model(path)
            overlap = copy.deepcopy(model)
            overlap['sites'][1]['center_bohr'] = overlap['sites'][0]['center_bohr']
            path.write_text(json.dumps(overlap))
            with self.assertRaisesRegex(ValueError, 'nonzero reference separation'):
                read_fragment_model(path)


@unittest.skipIf(np is None, 'install dynamics extra for fragment-site dynamics')
class FragmentPhysicsTests(unittest.TestCase):
    def setUp(self):
        params, _ = read_fragment_model(EXAMPLE.with_suffix('.json'))
        self.model = FragmentSitesModel(params)
        self.q = np.asarray(params.coordinates) + np.array([.3, -.1, .2, -.2, -.3, .1, .1, .1, -.3])

    def test_reference_hamiltonian_and_nonlinear_radial_law(self):
        p = self.model.parameters
        h0, _ = self.model.evaluate(p.coordinates)
        np.testing.assert_array_equal(np.diag(h0), p.energies)
        self.assertEqual(h0[0, 1], -.0036)
        self.assertEqual(h0[0, 2], 0.)
        h, _ = self.model.evaluate(self.q)
        r = np.linalg.norm(self.q[:3] - self.q[3:6])
        self.assertAlmostEqual(h[0, 1], -.0036 * np.exp(-1.2 * (r - p.contacts[0].reference)))
        collision = self.q.copy()
        collision[:3] = collision[3:6]
        with self.assertRaisesRegex(ValueError, 'undefined radial derivative'):
            self.model.evaluate(collision)

    def test_analytic_cartesian_hamiltonian_adiabatic_forces_and_nac(self):
        _, dh = self.model.evaluate(self.q)
        frame = adiabatic_frame(self.model, self.q)
        delta = 2.e-6
        for k, displacement in enumerate(np.eye(len(self.q)) * delta):
            plus, minus = (self.model.evaluate(self.q + dq)[0] for dq in (displacement, -displacement))
            np.testing.assert_allclose((plus - minus)/(2*delta), dh[..., k], atol=4e-11)
            right, left = (adiabatic_frame(self.model, self.q + dq, frame.vectors)
                           for dq in (displacement, -displacement))
            np.testing.assert_allclose((right.energies - left.energies)/(2*delta), frame.gradients[:, k], atol=5e-11)
            numerical_nac = frame.vectors.T @ ((right.vectors - left.vectors)/(2*delta))
            np.testing.assert_allclose(numerical_nac, frame.nac[..., k], atol=5e-9)

    def test_rigid_translation_and_rotation(self):
        h, dh = self.model.evaluate(self.q)
        points = self.q.reshape(-1, 3)
        rotation = np.array([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]])
        moved = (points @ rotation.T + [5., -4., 2.]).ravel()
        other_h, other_dh = self.model.evaluate(moved)
        np.testing.assert_allclose(other_h, h, atol=1e-15)
        grad = dh.reshape(3, 3, 3, 3)
        np.testing.assert_allclose(other_dh.reshape(grad.shape), grad @ rotation.T, atol=1e-15)
        np.testing.assert_allclose(grad.sum(axis=2), 0., atol=1e-15)

    def test_localized_ensemble_preserves_momentum_and_resume_checks(self):
        control = EXAMPLE.read_text().replace('"fragment-sites.json"', json.dumps(str(EXAMPLE.with_suffix('.json'))))
        control = control.replace('max-steps = 200', 'max-steps = 8').replace('trajectories = 16', 'trajectories = 3')
        with tempfile.TemporaryDirectory() as root, redirect_stdout(io.StringIO()):
            path = Path(root) / 'run.toml'
            path.write_text(control)
            directory = Path(root) / 'out'
            result = run_control(path, output=directory)
            self.assertEqual(result['status'], 'completed')
            self.assertEqual(result['site_ids'], ['A', 'B', 'C'])
            self.assertAlmostEqual(sum(result['initial_active_probabilities']), 1.)
            self.assertLess(result['max_norm_error'], 1e-11)
            self.assertLess(result['max_energy_drift_hartree'], 1e-8)
            rows = [json.loads(line) for line in (directory/'trajectories.jsonl').read_text().splitlines()]
            initial_p = np.asarray(result['config']['momenta']).reshape(-1, 3).sum(axis=0)
            for row in rows:
                if row['step'] == 0:
                    np.testing.assert_allclose(row['diabatic_populations'], [1., 0., 0.], atol=2e-14)
                np.testing.assert_allclose(np.array(row['p_au']).reshape(-1, 3).sum(axis=0), initial_p, atol=1e-10)
            self.assertTrue(run_control(path, output=directory, resume=True)['reused'])
