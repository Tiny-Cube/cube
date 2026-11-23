import copy
from importlib.util import find_spec
import json
import math
from pathlib import Path
import tempfile
import unittest

from wren.workflow.vertical_gradient import read_vertical_gradient


EXAMPLE = Path(__file__).resolve().parents[1] / 'examples/reorganization/vg-modes.json'


class VerticalGradientTest(unittest.TestCase):
    def setUp(self):
        self.work = tempfile.TemporaryDirectory()
        self.addCleanup(self.work.cleanup)
        self.data = json.loads(EXAMPLE.read_text())

    def analyze(self, data):
        path = Path(self.work.name) / 'vg.json'
        path.write_text(json.dumps(data))
        result = read_vertical_gradient(path)
        json.dumps(result, allow_nan=False)
        return result

    def projected(self):
        data = copy.deepcopy(self.data)
        del data['cartesian']
        # k = 4 eV/(Å² u); minimum at Q = +0.25 Å sqrt(u).
        # E(Q) = E_R + gQ + k Q²/2 has g = -1 and lower minimum by 0.125 eV.
        # Standard conversion for 1/(2πc)*sqrt(eV/(Å² u)) in cm^-1.
        data['modes'] = [{'id': 'first', 'wavenumber_cm1': 2 * 521.470899,
                          'gradient_ev_per_angstrom_sqrt_u': -1.0}]
        return data

    def test_known_harmonic_curvature_and_surface(self):
        data = self.projected()
        data['reference'] = 'b'
        result = self.analyze(data)
        mode = result['modes'][0]
        self.assertEqual((result['reference'], result['target']), ('b', 'a'))
        self.assertAlmostEqual(mode['displacement_angstrom_sqrt_u'], .25, places=6)
        self.assertAlmostEqual(mode['curvature_ev_per_angstrom2_u'], 4., places=5)
        self.assertAlmostEqual(result['relaxation_ev'], .125, places=6)
        self.assertAlmostEqual(mode['huang_rhys'] * mode['quantum_ev'], result['relaxation_ev'])
        self.assertIsNone(result['projection']['unrepresented_gradient_norm_ev_per_angstrom_sqrt_u'])

    @unittest.skipUnless(find_spec('numpy'), 'Cartesian projection requires numpy')
    def test_mass_weighting_and_remaining_gradient(self):
        self.data['cartesian'] = {'atom_ids': ['anchor'], 'masses_u': [4.0],
                                  'gradient_ev_per_angstrom': [[2.0, 1.0, 0.0]]}
        self.data['modes'][0]['mass_weighted_direction'] = [[-1.0, 0.0, 0.0]]
        result = self.analyze(self.data)
        self.assertAlmostEqual(result['modes'][0]['gradient_ev_per_angstrom_sqrt_u'], -1.0)
        self.assertAlmostEqual(result['projection']['unrepresented_gradient_norm_ev_per_angstrom_sqrt_u'], .5)
        self.assertEqual(result['projection']['cartesian_dimension'], 3)
        self.data['modes'].append({**self.data['modes'][0], 'id': 'duplicate-vector'})
        with self.assertRaisesRegex(ValueError, 'orthonormal'):
            self.analyze(self.data)

    def test_degenerate_orthogonal_subspace_invariant(self):
        data = self.projected()
        data['modes'].append({**data['modes'][0], 'id': 'second', 'gradient_ev_per_angstrom_sqrt_u': 0.})
        baseline = self.analyze(data)['relaxation_ev']
        for mode in data['modes']:
            mode['gradient_ev_per_angstrom_sqrt_u'] = 1 / math.sqrt(2)
        self.assertAlmostEqual(self.analyze(data)['relaxation_ev'], baseline)

    def test_zero_imaginary_nonfinite_wrong_units_and_unconfirmed_stationarity(self):
        for key, value in (('wavenumber_cm1', 0), ('wavenumber_cm1', -100),
                           ('wavenumber_cm1', math.nan), ('gradient_ev_per_angstrom_sqrt_u', True),
                           ('gradient_ev_per_angstrom', -1)):
            data = self.projected()
            data['modes'][0][key] = value
            with self.subTest(key=key, value=value), self.assertRaises(ValueError):
                self.analyze(data)
        data = self.projected()
        data['reference_stationary'] = False
        with self.assertRaisesRegex(ValueError, 'stationary'):
            self.analyze(data)


if __name__ == '__main__':
    unittest.main()
