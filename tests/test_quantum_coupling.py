"""Optional quantum algebra tests; no SCF is needed for these invariants."""
import os
from pathlib import Path
import tempfile
import unittest

try:
    import numpy as np
except ImportError:
    np = None


@unittest.skipIf(np is None, 'install the quantum extra for NumPy algebra tests')
class CouplingAlgebraTests(unittest.TestCase):
    def setUp(self):
        from universal_kmc.quantum import coupling
        self.c = coupling

    def test_two_state_overlap_formula_and_energy_origin(self):
        for s in (-.6, -.1, 0, .3, .8):
            metric = np.array([[1., s], [s, 1.]])
            h = np.array([[-.4, .023], [.023, -.31]])
            effective, x, w = self.c.lowdin(h, metric)
            expected = (h[0, 1]-s*np.trace(h)/2)/(1-s*s)
            self.assertAlmostEqual(effective[0, 1], expected, places=13)
            np.testing.assert_allclose(x@metric@x, np.eye(2), atol=2e-14)
            shifted, _, _ = self.c.lowdin(h+4.73*metric, metric)
            np.testing.assert_allclose(shifted, effective+4.73*np.eye(2), atol=2e-14)

    def test_multistate_phase_permutation_and_generalized_spectrum(self):
        rng = np.random.default_rng(1729)
        raw = rng.normal(size=(5, 5))
        metric = raw.T@raw + 2*np.eye(5)
        raw = rng.normal(size=(5, 5))
        h = raw+raw.T
        eff, _, _ = self.c.lowdin(h, metric)
        np.testing.assert_allclose(np.sort(np.linalg.eigvals(np.linalg.solve(metric, h))),
                                   np.linalg.eigvalsh(eff), atol=1e-13)
        transform = np.eye(5)[:, [2, 4, 1, 3, 0]]@np.diag([-1, 1, -1, 1, 1])
        changed, _, _ = self.c.lowdin(transform.T@h@transform, transform.T@metric@transform)
        np.testing.assert_allclose(changed, transform.T@eff@transform, atol=1e-13)

    def test_singular_invalid_and_nonsymmetric_metrics_rejected(self):
        for metric in (np.ones((2, 2)), np.diag([1., -1.]), np.diag([1., 1e-12])):
            with self.assertRaises(self.c.ApplicabilityError):
                self.c.lowdin(np.eye(2), metric)
        for h, s in (([[1, 2], [3, 4]], np.eye(2)),
                     (np.eye(2), [[1, 2], [0, 1]]),
                     ([[float('nan')]], [[1.]]), ([[1j]], [[1.]])):
            with self.assertRaises(ValueError):
                self.c.lowdin(h, s)

    def test_gmh_recovers_known_diabatic_hamiltonian(self):
        for delta in (0., .03, .7):
            h = np.array([[delta, .017], [.017, -.13]])
            energies, u = np.linalg.eigh(h)
            dipole_diabatic = np.zeros((3, 2, 2))
            dipole_diabatic[:, 0, 0] = [-1., 2., .3]
            dipole_diabatic[:, 1, 1] = [3., -1., 2.]
            mu = np.einsum('pi,xpq,qj->xij', u, dipole_diabatic, u)
            self.assertAlmostEqual(self.c.gmh(energies, mu)['magnitude'], .017, places=14)
            mu += np.array([14., -6., 21.])[:, None, None]*np.eye(2)
            rotated = np.einsum('xy,yij->xij', np.array([[0, -1, 0], [0, 0, 1], [-1, 0, 0]]), mu)
            self.assertAlmostEqual(self.c.gmh(energies+70., rotated)['magnitude'], .017, places=13)
            phase = np.diag([-1., 1.])
            mu = np.einsum('pi,xpq,qj->xij', phase, mu, phase)
            self.assertAlmostEqual(self.c.gmh(energies, mu)['magnitude'], .017, places=13)

    def test_gmh_symmetric_limit_and_undefined_dipoles(self):
        mu = np.zeros((3, 2, 2))
        mu[0, 0, 1] = mu[0, 1, 0] = 2.
        self.assertAlmostEqual(self.c.gmh([-.1, .1], mu)['magnitude'], .1)
        with self.assertRaises(self.c.ApplicabilityError):
            self.c.gmh([-.1, .1], np.zeros((3, 2, 2)))
        with self.assertRaises(self.c.ApplicabilityError):
            self.c.gmh([-.1, .1], mu, axis=[0., 1., 0.])

    def test_esid_requires_equivalence_and_site_check(self):
        self.assertAlmostEqual(self.c.energy_splitting([-.1, .1], [-.3, -.3], equivalent=True), .1)
        for kwargs, sites in (({}, [-.3, -.3]), ({'equivalent': True}, [-.3, -.2])):
            with self.assertRaises(self.c.ApplicabilityError):
                self.c.energy_splitting([-.1, .1], sites, **kwargs)

    def test_input_and_config_validation(self):
        from universal_kmc.quantum.transfer import SCFConfig, read_xyz
        for kwargs in ({'conv_tol': float('nan')}, {'max_cycle': True}, {'backend': 'silent-fallback'}, {'max_pair_leakage': 1}):
            with self.assertRaises(ValueError):
                SCFConfig(**kwargs).validate()
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'geometry.xyz'
            path.write_text('1\ncomment\nH 0 0 0\n1\nextra frame\nH 1 0 0\n')
            with self.assertRaises(ValueError):
                read_xyz(path)

    @unittest.skipUnless(os.environ.get('CUBE_GPU_TESTS') == '1', 'explicit GPU opt-in required')
    def test_gpu_lowdin_matches_numpy(self):
        import cupy as cp
        h, s = np.array([[-.3, .01], [.01, -.2]]), np.array([[1, .2], [.2, 1]])
        gpu, _, _ = self.c.lowdin(cp.asarray(h), cp.asarray(s))
        self.assertIsInstance(gpu, cp.ndarray)
        cpu, _, _ = self.c.lowdin(h, s)
        np.testing.assert_allclose(gpu.get(), cpu, atol=1e-13)


if __name__ == '__main__':
    unittest.main()
