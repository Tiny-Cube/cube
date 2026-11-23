from importlib.util import find_spec
import json
from pathlib import Path
import tempfile
import unittest

from wren.surface_hopping.casscf_input import read_casscf_model
from wren.surface_hopping.control import load_control

EXAMPLES = Path(__file__).resolve().parents[1]/'examples/surface-hopping'


class FragmentInput(unittest.TestCase):
    def test_partition_and_carrier_sign(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)
            cfg = json.loads((EXAMPLES/'ethene-hole.json').read_text())
            cfg['xyz'] = str(EXAMPLES/'ethene-dimer.xyz')
            cfg['fragments']['B'][0] = 1
            file = path/'ethene-hole.json'
            file.write_text(json.dumps(cfg))
            with self.assertRaisesRegex(ValueError, 'partition'):
                read_casscf_model(file)
            cfg['fragments']['B'][0] = 7
            file.write_text(json.dumps(cfg))
            request = (EXAMPLES/'ethene-hole.toml').read_text()
            control = path/'input.toml'
            control.write_text(request.replace('hole-fragment', 'electron-fragment'))
            with self.assertRaisesRegex(ValueError, 'charge'):
                load_control(control)
            control.write_text(request.replace('hole-fragment = "A"', 'hole-fragment = "A"\nactive-state=1'))
            with self.assertRaisesRegex(ValueError,'excluding active-state'):
                load_control(control)


@unittest.skipUnless(find_spec('pyscf') and find_spec('numpy'), 'needs PySCF and NumPy')
class DoubletFragmentCharges(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from wren.surface_hopping.pyscf_casscf import SACASSCF
        cfg, cls.q, _, _ = read_casscf_model(EXAMPLES/'lih-casscf.json')
        cls.provider = SACASSCF(dict(cfg, charge=1, multiplicity=2, active_electrons=1, use_etfs=False,
                                     fragments={'Li': [1], 'H': [2]}, orbital_gradient_tolerance=1e-6))
        cls.frame = cls.provider.evaluate_adiabatic(cls.q)

    def test_doublet_response_against_independent_determinant_difference(self):
        import numpy as np
        from pyscf import gto
        frame = self.frame
        dq = np.zeros_like(self.q)
        dq[-1] = .002
        plus = self.provider.evaluate_adiabatic(self.q + dq, frame)
        minus = self.provider.evaluate_adiabatic(self.q - dq, frame)
        np.testing.assert_allclose(frame.diagnostics['spin_squared'], [.75,.75], atol=1e-12)
        np.testing.assert_allclose((plus.energies-minus.energies)/.004,frame.gradients[:,-1],rtol=0,atol=1e-5)

        def overlap(other):
            left = frame.wavefunctions
            right = other.wavefunctions
            s = left.mo.T @ gto.intor_cross('int1e_ovlp', left.mol, right.mol) @ right.mo
            # LiH+: one core beta, one core and one active alpha; enumerate full determinants.
            a = np.array([[np.linalg.det(s[np.ix_([0,i+1], [0,j+1])]) for j in range(2)] for i in range(2)])
            return s[0,0] * left.ci[0].ravel() @ a @ right.ci[1].ravel()
        self.assertAlmostEqual((overlap(plus)-overlap(minus))/.004, frame.nac[0,1,-1],delta=2e-4)

    def test_coherence_against_independent_ao_density_and_localized_extremum(self):
        import numpy as np
        from pyscf import fci
        from dataclasses import replace
        from wren.surface_hopping.fragment_charges import charge_operators, fragment_observables, localize_charge
        frame = self.frame
        ops = frame.diagnostics['fragment_charge_operators_e']
        nuclei = frame.diagnostics['fragment_nuclear_charges_e']
        c = np.array([.5, .5+np.sqrt(.5)*1j])
        charge = fragment_observables(c, ops, nuclei)['fragment_net_charges_e']
        self.assertAlmostEqual(sum(charge.values()), 1.)
        w = frame.wavefunctions
        combined = sum(ci*amp for ci, amp in zip(w.ci, c))
        density = sum(fci.direct_spin1.make_rdm1(z.copy(), w.ncas, w.nelecas)
                      for z in (combined.real, combined.imag))
        core = w.mo[:,:w.ncore]
        active = w.mo[:,w.ncore:w.ncore+w.ncas]
        dao = 2 * core @ core.T + active @ density @ active.T
        eigen, vectors = np.linalg.eigh(w.mol.intor_symmetric('int1e_ovlp'))
        sqrt_s = (vectors*np.sqrt(eigen)) @ vectors.T
        ao = (sqrt_s @ dao @ sqrt_s).diagonal()
        last = w.mol.aoslice_by_atom()[0,3]
        self.assertAlmostEqual(3.-sum(ao[:last]),charge['Li'],places=11)
        mixture = np.dot(abs(c)**2, np.diag(ops['Li']))
        self.assertGreater(abs(charge['Li']-mixture),1e-3)
        signs = np.array([-1., 1.])
        flipped = replace(w, ci=[z*ci for z,ci in zip(signs, w.ci)])
        operators_flipped, _ = charge_operators(flipped, {'Li': [1], 'H': [2]})
        in_new_gauge = fragment_observables(c*signs, operators_flipped, nuclei)
        self.assertAlmostEqual(in_new_gauge['fragment_net_charges_e']['Li'], charge['Li'])
        localized, optimum = localize_charge(ops, 'H', 'hole')
        self.assertAlmostEqual(optimum, fragment_observables(localized,ops,nuclei)['fragment_net_charges_e']['H'])
        self.assertGreater(optimum,max(np.diag(ops['H'])))


if __name__ == '__main__':
    unittest.main()
