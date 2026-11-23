"""DOS sum rules, analytic Bloch spectra and finite-sample disorder identities."""
import unittest
try:
    import numpy as np
    from wren.dos import broaden, bloch_spectra, disorder_statistics, kmesh
except ImportError:
    np = None

@unittest.skipIf(np is None, 'DOS NumPy dependency required')
class DosCoreTests(unittest.TestCase):
    def test_integral_and_moments_include_spin_and_fractional_weights(self):
        x = np.linspace(-3, 4, 14001)
        weights = np.array([.5, 1.5])
        dos = broaden([-1., 2.], x, .1, weights=weights, chunk=67)
        self.assertAlmostEqual(float(np.trapezoid(dos, x)), 2., places=11)
        self.assertAlmostEqual(float(np.trapezoid(dos*x, x))/2, 1.25, places=10)
        self.assertAlmostEqual(float(np.trapezoid(dos*x*x, x))/2, 3.25+.01, places=10)

    def test_periodic_chain_has_cosine_dispersion_and_correct_moments(self):
        points, weights = kmesh([512, 1, 1])
        e = bloch_spectra([1.], [(0, 0, [1, 0, 0], -.2)], points)[:, 0]
        np.testing.assert_allclose(e, 1.-.4*np.cos(2*np.pi*points[:, 0]), atol=1e-14)
        self.assertAlmostEqual(float(weights@e), 1.)
        self.assertAlmostEqual(float(weights@(e-1.)**2), .08)

    def test_dimer_and_interference_need_signed_couplings(self):
        points, _ = kmesh([2, 1, 1])
        positive = bloch_spectra([0., 0.], [(0, 1, [0, 0, 0], .1), (0, 1, [1, 0, 0], .1)], points)
        negative = bloch_spectra([0., 0.], [(0, 1, [0, 0, 0], .1), (0, 1, [1, 0, 0], -.1)], points)
        np.testing.assert_allclose(positive, [[-.2, .2], [0., 0.]], atol=1e-16)
        np.testing.assert_allclose(negative, positive[::-1], atol=1e-16)
        with self.assertRaisesRegex(ValueError, 'double count'):
            bloch_spectra([0., 0.], [(0, 1, [1, 0, 0], .1), (1, 0, [-1, 0, 0], .1)], points)

    def test_static_dynamic_total_variance_identity(self):
        levels = np.array([[1, 3, 5], [2, 2, 7], [-1, 4, 6]], dtype=float)
        result = disorder_statistics(levels)
        self.assertAlmostEqual(result['variance_total_ev2'], result['variance_static_ev2']+result['variance_dynamic_ev2'])
        self.assertFalse(result['statistical_uncertainty_estimated'])

    def test_invalid_inputs_and_no_finite_window_rescaling(self):
        self.assertLess(float(np.trapezoid(broaden([0.], [-.01, .01], .1), [-.01, .01])), .1)
        for args in (([float('nan')], [-1., 1.], .1), ([0.], [1., -1.], .1), ([0.], [-1., 1.], 0.)):
            with self.assertRaises(ValueError):
                broaden(*args)


if __name__ == '__main__':
    unittest.main()
