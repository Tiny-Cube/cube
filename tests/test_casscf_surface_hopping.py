from importlib.util import find_spec
import json
from pathlib import Path
import tempfile
import unittest

from wren.surface_hopping.control import load_control


ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = ROOT/'examples/surface-hopping/lih-casscf.toml'


@unittest.skipUnless(find_spec('pyscf') and find_spec('numpy'), 'requires PySCF and NumPy')
class CASSCFDynamics(unittest.TestCase):
    def test_raw_signed_coupling_and_state_forces_against_wavefunction_differences(self):
        import numpy as np
        from pyscf import gto
        from wren.surface_hopping.casscf_input import read_casscf_model
        from wren.surface_hopping.pyscf_casscf import SACASSCF, wavefunction_overlaps

        settings, q, _, _ = read_casscf_model(EXAMPLE.with_suffix('.json'))
        model = SACASSCF(dict(settings, use_etfs=False))
        q = np.array(q)
        center = model.evaluate_adiabatic(q)
        displacement = np.zeros_like(q)
        displacement[-1] = .002
        plus = model.evaluate_adiabatic(q+displacement, center)
        minus = model.evaluate_adiabatic(q-displacement, center)
        fd_grad = (plus.energies - minus.energies)/.004
        np.testing.assert_allclose(fd_grad, center.gradients[:, -1], atol=1e-5, rtol=0)

        def explicit_overlap(left, right):
            # Independently enumerate full Slater determinants: one closed-shell
            # Li 1s core plus one active electron of each spin in two orbitals.
            cross = left.mo.T @ gto.intor_cross('int1e_ovlp', left.mol, right.mol) @ right.mo
            det = np.array([[np.linalg.det(cross[np.ix_([0,i+1],[0,j+1])]) for j in range(2)]
                            for i in range(2)])
            return np.array([[np.einsum('ab,ac,bd,cd', a, det, det, b) for b in right.ci] for a in left.ci])

        left, right = center.wavefunctions, plus.wavefunctions
        ov_plus = explicit_overlap(left, right)
        np.testing.assert_allclose(wavefunction_overlaps(left, right), ov_plus, atol=1e-12, rtol=0)
        ov_minus = explicit_overlap(left, minus.wavefunctions)
        fd_nac = (ov_plus[0,1] - ov_minus[0,1])/.004
        self.assertGreater(abs(fd_nac), .05)
        self.assertAlmostEqual(fd_nac, center.nac[0,1,-1], delta=3e-4)
        np.testing.assert_allclose(center.nac, -center.nac.swapaxes(0,1))
        self.assertIsNone(center.vectors)

    def test_native_trajectory_has_qc_evidence_and_no_diabatic_population(self):
        import numpy as np
        from wren.surface_hopping.runner import run_control
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)
            for ext in ('.toml','.json'):
                (path/('lih-casscf'+ext)).write_text(EXAMPLE.with_suffix(ext).read_text().replace(
                    'max-steps = 16','max-steps = 2'))
            (path/'lih.xyz').write_bytes((EXAMPLE.parent/'lih.xyz').read_bytes())
            result = run_control(path/EXAMPLE.name, output=path/'out')
            self.assertEqual(result['status'], 'completed')
            self.assertEqual(result['scope'], 'on_the_fly_molecular_sa_casscf')
            self.assertNotIn('mean_final_diabatic_populations', result)
            sample = result['final_states'][0]
            self.assertNotIn('diabatic_populations', sample)
            self.assertGreater(sample['adiabatic_populations'][0], 0)
            self.assertGreater(np.min(np.diag(sample['quantum']['adjacent_wavefunction_overlaps'])), .99)
            d = np.asarray(sample['derivative_couplings_per_bohr'])[0,1].reshape(2,3)
            np.testing.assert_allclose(d.sum(axis=0), 0, atol=1e-7)
            self.assertLess(result['max_energy_drift_hartree'], 1e-5)
            self.assertTrue(run_control(path/EXAMPLE.name, output=path/'out', resume=True)['reused'])
            (path/'lih.xyz').write_text((path/'lih.xyz').read_text().replace('1.852', '1.851'))
            with self.assertRaisesRegex(ValueError, 'identical'):
                run_control(path/EXAMPLE.name, output=path/'out', resume=True)

    def test_overlapping_roots_fail_with_wrong_tracking_and_with_small_gap(self):
        import numpy as np
        from unittest.mock import patch
        from wren.surface_hopping.pyscf_casscf import align_roots, SACASSCF, Wavefunctions
        from wren.surface_hopping.casscf_input import read_casscf_model
        mock = Wavefunctions(None, np.eye(2), [np.array([1]),np.array([1])],0,2,(1,1))
        with patch('wren.surface_hopping.pyscf_casscf.wavefunction_overlaps', return_value=np.array([[.7,.8],[.8,.7]])):
            with self.assertRaisesRegex(RuntimeError, 'identity'):
                align_roots(mock,mock,.6)
        config,q,_,_=read_casscf_model(EXAMPLE.with_suffix('.json'))
        with self.assertRaisesRegex(RuntimeError, 'near-degenerate'):
            SACASSCF(dict(config,min_gap_hartree=2)).evaluate_adiabatic(q)

    def test_nonzero_qc_flux_and_downward_hop_conserve_energy_and_momentum(self):
        from unittest.mock import Mock
        import numpy as np
        from wren.surface_hopping.casscf_input import read_casscf_model
        from wren.surface_hopping.pyscf_casscf import SACASSCF
        from wren.surface_hopping.engine import FSSH
        config, q, masses, _ = read_casscf_model(EXAMPLE.with_suffix('.json'))
        engine = FSSH(SACASSCF(config), masses, dt=.5)
        start = engine.initialize(q, [0,0,-20,0,0,20], active=1)
        # Controlled uniform draw selects the nonzero rare-hop branch. This
        # checks momentum adjustment with an actual QC NAC; it is not a rate.
        rng = Mock()
        rng.random.return_value = 0.
        end, event = engine.step(start, rng)
        self.assertEqual((event['event'], end.active), ('accepted', 0))
        self.assertGreater(event['probability_sum'], 0.)
        self.assertLess(abs(event['jump_energy_error']), 1e-12)
        np.testing.assert_allclose(end.p.reshape(-1,3).sum(axis=0), 0, atol=1e-8)


class CASSCFControl(unittest.TestCase):
    def test_check_hashes_xyz_and_versions_without_qc_run(self):
        config, _, provenance, _ = load_control(EXAMPLE)
        self.assertEqual(config.model, 'pyscf-sa-casscf')
        self.assertAlmostEqual(config.coordinates[-1], 3.5)
        self.assertIn(str(EXAMPLE.parent/'lih.xyz'), provenance['sources'])
        self.assertIn('pyscf', provenance['versions'])
        self.assertIn('pyscf_casscf.py', provenance['implementation_sha256'])


if __name__ == '__main__':
    unittest.main()
