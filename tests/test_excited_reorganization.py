from importlib.util import find_spec
import json
from pathlib import Path
import tempfile
import unittest


def singlet_settings(**extra):
    settings = {'qc-engine': 'pyscf', 'environment': 'vacuum', 're-method': '4-points',
                'method': 'hf', 'basis-set': '6-31g', 'auxbasis': 'none', 'threads': 1,
                'grid-level': 3, 'conv-tol': 1e-10, 'max-cycle': 100, 'max-memory-mb': 800,
                'opt-maxsteps': 30}
    return {**settings, **extra}


@unittest.skipUnless(find_spec('numpy'), 'requires numpy')
class SingletSafeguards(unittest.TestCase):
    def test_collapsed_gap_and_root_character_swaps_stop_calculation(self):
        from wren.workflow.excited_reorganization import _screen_roots
        for energies, overlaps in (([2., 2.01], None), ([-.02, 2.], None),
                                   ([.01, 2.], None), ([2., 4.], [.02, .99]),
                                   ([2., 4.], [.70, .71])):
            with self.subTest(energies=energies, overlaps=overlaps), self.assertRaises(RuntimeError):
                _screen_roots(energies, overlaps, .05, .7)


@unittest.skipUnless(find_spec('geometric') and find_spec('pyscf'), 'requires PySCF and geomeTRIC')
class PhysicalSingletFourPoint(unittest.TestCase):
    def test_h2_two_optimizations_against_explicit_tda_matrix(self):
        import numpy as np
        from pyscf import gto, scf
        from pyscf.data.nist import HARTREE2EV
        from wren.workflow.excited_reorganization import calculate_excited_reorganization

        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)
            xyz = path / 'h2.xyz'
            xyz.write_text('2\nH2 starting structure, Angstrom\nH 0 0 -.39\nH 0 0 .39\n')
            output = path / 'result'
            result = calculate_excited_reorganization(xyz, singlet_settings(), output)
            self.assertEqual(result['status'], 'completed')
            energies = result['analysis']['energies']
            for label in ('RA', 'RB'):
                mol = gto.M(atom=str(output / f'{label}.xyz'), basis='6-31g', unit='Angstrom', verbose=0)
                mf = scf.RHF(mol).run(conv_tol=1e-12)
                response_a, _ = mf.TDA().get_ab()
                n = response_a.shape[0] * response_a.shape[1]
                omega = np.linalg.eigvalsh(response_a.reshape(n, n))[0]
                self.assertAlmostEqual(energies[f'E_A_{label}'], mf.e_tot, places=8)
                # The SCF energy threshold does not set the same tolerance
                # for orbital eigenvalues entering two independent TDA solves.
                self.assertAlmostEqual(energies[f'E_B_{label}'], mf.e_tot + omega, delta=1e-6)
            for key, value in result['analysis']['relaxation'].items():
                self.assertGreater(value, 0)
                self.assertAlmostEqual(result['relaxation_ev'][key], value * HARTREE2EV)
            trace = json.loads((output / 'RB.trace.json').read_text())
            self.assertGreater(len(trace), 1)
            self.assertTrue(all(step['transition_amplitude_overlaps'][0] >= .7 for step in trace))
            self.assertFalse(result['analysis']['diagnostics']['minima_verified'])
            self.assertTrue(calculate_excited_reorganization(xyz, singlet_settings(), output,
                                                              resume=True)['reused'])
            (output / 'RA.xyz').write_text('modified')
            with self.assertRaisesRegex(ValueError, 'changed or missing'):
                calculate_excited_reorganization(xyz, singlet_settings(), output, resume=True)

    def test_unsupported_combinations_rejected_before_output(self):
        from wren.workflow.excited_reorganization import calculate_excited_reorganization
        for overrides in ({'qc-engine': 'gpu-pyscf'}, {'environment': 'qmmm'}, {'auxbasis': 'weigend'}):
            with self.subTest(overrides=overrides), tempfile.TemporaryDirectory() as tmp:
                out = Path(tmp) / 'rejected'
                with self.assertRaises(ValueError):
                    calculate_excited_reorganization([], singlet_settings(**overrides), out)
                self.assertFalse(out.exists())


if __name__ == '__main__':
    unittest.main()
