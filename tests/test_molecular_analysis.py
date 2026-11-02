"""Optional analytical checks for the real-system validation script (needs numpy)."""
import importlib.util
from pathlib import Path
import unittest

HAS_NUMPY = importlib.util.find_spec('numpy') is not None


@unittest.skipUnless(HAS_NUMPY, 'optional molecular analysis requires numpy')
class MolecularBenchmarkTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import numpy as np
        spec = importlib.util.spec_from_file_location('molecular_validation', Path(__file__).resolve().parents[1]/'scripts/run_molecular_kmc.py')
        cls.module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.module)
        cls.np = np

    def test_closed_dimer_has_zero_diffusion_and_known_transient(self):
        m, np = self.module, self.np
        edges = [m.Edge(0, 1, 3, (1e-9, 0, 0)), m.Edge(1, 0, 3, (-1e-9, 0, 0))]
        grid = np.array([0, .1, .5, 1.])
        tensor, msd, residual = m.exact_benchmark(2, edges, grid)
        np.testing.assert_allclose(tensor, 0, atol=1e-28)
        np.testing.assert_allclose(msd, .5*(1-np.exp(-6*grid)), atol=1e-14)
        self.assertLess(residual, 1e-14)

    def test_periodic_walk_retains_winding_displacements(self):
        m, np = self.module, self.np
        # Two internal sites; each has independent left/right channels to the other.
        edges = [m.Edge(i, 1-i, 2, (dx*1e-9, 0, 0)) for i in range(2) for dx in [-1, 1]]
        grid = np.array([0, .1, .5, 1.])
        tensor, msd, _ = m.exact_benchmark(2, edges, grid)
        np.testing.assert_allclose(tensor, np.diag([2., 0, 0]), atol=1e-14)
        np.testing.assert_allclose(msd, 4*grid, atol=1e-14)


if __name__ == '__main__':
    unittest.main()
