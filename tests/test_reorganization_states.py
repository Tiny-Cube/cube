"""Two-surface physical invariants, rather than state-specific formula copies."""
from dataclasses import replace
import json
import math
import subprocess
import sys
import unittest

from wren.reorganization import (
    ElectronicState, FourPointEnergies, analyze_four_point_pair, four_point_relaxations,
)


class StatePairReorganizationTests(unittest.TestCase):
    def setUp(self):
        self.a = ElectronicState('S0', 0, 1, solver='scf')
        self.b = ElectronicState('S1', 0, 1, root=1, solver='tddft')
        # A(q)=q^2/2; B(q)=3(q-3)^2/2+2. Different curvatures mean
        # different relaxation energies; no equal-curvature assumption is valid.
        self.energies = FourPointEnergies(0., 15.5, 2., 4.5, unit='eV')

    def test_asymmetric_harmonic_surfaces_and_gap_identity(self):
        result = analyze_four_point_pair(self.a, self.b, self.energies)
        self.assertEqual(result['relaxation'], {'on_a': 4.5, 'on_b': 13.5, 'sum': 18.})
        self.assertEqual(result['gaps']['adiabatic'], 2.)
        self.assertEqual(result['gaps']['vertical_at_ra']-result['gaps']['vertical_at_rb'],
                         result['relaxation']['sum'])
        self.assertFalse(result['diagnostics']['minima_verified'])
        self.assertEqual(result['states']['b']['root'], 1)
        self.assertEqual(json.loads(json.dumps(result, allow_nan=False)), result)

    def test_reversing_state_pair_swaps_terms_and_reverses_driving_energy(self):
        e = self.energies
        forward = analyze_four_point_pair(self.a, self.b, e)
        reverse = analyze_four_point_pair(self.b, self.a, FourPointEnergies(
            e.e_b_rb, e.e_a_rb, e.e_a_ra, e.e_b_ra, e.unit))
        self.assertEqual(reverse['relaxation']['on_a'], forward['relaxation']['on_b'])
        self.assertEqual(reverse['relaxation']['on_b'], forward['relaxation']['on_a'])
        self.assertEqual(reverse['relaxation']['sum'], forward['relaxation']['sum'])
        self.assertEqual(reverse['gaps']['adiabatic'], -forward['gaps']['adiabatic'])

    def test_energy_origin_and_state_offset_invariants(self):
        e = self.energies
        result = analyze_four_point_pair(self.a, self.b, e)
        shifted = analyze_four_point_pair(self.a, self.b, FourPointEnergies(
            e.e_a_ra+100., e.e_b_ra+100., e.e_b_rb+100., e.e_a_rb+100., e.unit))
        self.assertEqual(result['relaxation'], shifted['relaxation'])
        self.assertEqual(result['gaps'], shifted['gaps'])
        # Offsetting only B changes the driving energy, not its deformation cost.
        offset = analyze_four_point_pair(self.a, self.b, replace(
            e, e_b_ra=e.e_b_ra+7., e_b_rb=e.e_b_rb+7.))
        self.assertEqual(result['relaxation'], offset['relaxation'])
        self.assertEqual(offset['gaps']['adiabatic'], result['gaps']['adiabatic']+7.)

    def test_charged_triplet_and_localized_states_use_same_analysis(self):
        for state in (ElectronicState('cation', 1, 2, solver='scf'),
                      ElectronicState('anion', -1, 2, solver='scf'),
                      ElectronicState('T1', 0, 3, root=1, solver='tddft'),
                      ElectronicState('localized_B', 0, 1)):
            result = analyze_four_point_pair(self.a, state, self.energies)
            self.assertEqual(result['relaxation']['sum'], 18.)
        # This is external-energy analysis, not proof any of these solvers ran.
        self.assertFalse(result['diagnostics']['state_identity_verified'])

    def test_negative_terms_are_reported_without_abs_or_clipping(self):
        result = analyze_four_point_pair(self.a, self.b, replace(self.energies, e_a_rb=-1.))
        self.assertEqual(result['relaxation']['on_a'], -1.)
        self.assertTrue(result['diagnostics']['negative_relaxation_a'])

    def test_invalid_units_states_and_nonfinite_inputs(self):
        for value in (True, '1', complex(1, 0), math.nan, math.inf):
            with self.assertRaises(ValueError):
                FourPointEnergies(value, 1., 0., 1., 'eV')
        for unit in ('', 'kJ/mol', 'mixed'):
            with self.assertRaises(ValueError):
                replace(self.energies, unit=unit)
        for kw in ({'name': ''}, {'charge': True}, {'multiplicity': 0}, {'root': -1}, {'solver': ''}):
            with self.assertRaises(ValueError):
                replace(self.a, **kw)
        with self.assertRaises(ValueError):
            analyze_four_point_pair(self.a, self.a, self.energies)
        with self.assertRaisesRegex(ValueError, 'overflowed'):
            four_point_relaxations(-1e308, 1e308, 0., 1e308)

    def test_algebra_import_has_no_quantum_or_array_dependencies(self):
        script = '''import sys
from wren.reorganization import ElectronicState, FourPointEnergies, analyze_four_point_pair
assert all(x not in sys.modules for x in ('numpy','pyscf','cupy','MDAnalysis'))
r=analyze_four_point_pair(ElectronicState('A',0,1),ElectronicState('B',0,1),FourPointEnergies(0,2,1,.5,'eV'))
assert r['relaxation']['sum']==1.5
'''
        result = subprocess.run([sys.executable, '-c', script], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)


if __name__ == '__main__':
    unittest.main()
