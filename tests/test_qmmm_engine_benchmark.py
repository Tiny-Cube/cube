"""Energy convention and derivative sign at the Gaussian/PySCF interface."""
import importlib.util
from pathlib import Path
import unittest

try:
    import numpy as np
    spec=importlib.util.spec_from_file_location('qmmm_engine_benchmark',Path(__file__).resolve().parents[1]/'scripts/benchmark_qmmm_engines.py')
    benchmark=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(benchmark)
except ImportError:
    np=None


@unittest.skipIf(np is None,'quantum workflow dependencies required')
class GaussianConventionTests(unittest.TestCase):
    # Numeric fields from the real H2/STO-3G + q=+0.1/-0.1 at y=+3/-3 A probe.
    log='''Self energy of the charges = -0.0008819620 a.u.
SCF Done: E(RHF) = -1.11764126953 A.U. after 1 cycles
Normal termination of Gaussian 16
'''
    fchk='''Number of basis functions                  I                2
Total Energy                               R     -1.117641269528463E+00
Cartesian Gradient                         R   N=           6
 -2.7679599400D-02 -4.4305844000D-05 0.0  2.7679599400D-02 -4.4305844000D-05 0.0
'''

    def test_mm_self_energy_removed_and_gradient_not_negated(self):
        result=benchmark.gaussian_result(self.log,self.fchk,2,gradient=True)
        self.assertAlmostEqual(result['electronic_energy_hartree'],-1.1167593073964255,places=8)
        self.assertLess(result['gradient_hartree_per_bohr'][0][0],0)
        g=np.array(result['gradient_hartree_per_bohr'])
        self.assertEqual(g[:,0].sum(),0)
        # A frozen external field may exert a net force on QM; do not remove it.
        self.assertAlmostEqual(g[:,1].sum(),-8.8611688e-5)

    def test_unconverged_or_incomplete_quantities_are_not_accepted(self):
        with self.assertRaises(RuntimeError):
            benchmark.gaussian_result(self.log+'Error termination',self.fchk,2,gradient=True)
        with self.assertRaises(RuntimeError):
            benchmark.gaussian_result(self.log,self.fchk,3,gradient=True)
        with self.assertRaises(ValueError):
            benchmark.fchk_field('Cartesian Gradient R N= 6\n .1 .2','Cartesian Gradient')


if __name__=='__main__':
    unittest.main()
