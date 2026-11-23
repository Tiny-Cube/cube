"""Import foreign four-point energies in their original units and state convention."""
import copy
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from wren.workflow.external_reorganization import read_external_four_points


ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = ROOT / 'examples/reorganization/external-points.json'


class ExternalFourPointTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.source = Path(temp.name) / 'dataset.json'
        self.pair = {
            'id': 'mol-1',
            'states': {
                'a': {'name': 'S0', 'charge': 0, 'multiplicity': 1, 'solver': 'KS', 'root': 0},
                'b': {'name': 'S1', 'charge': 0, 'multiplicity': 1, 'solver': 'LR-TDDFT', 'root': 1}},
            'energies': {'E_A_RA': -50.0, 'E_B_RA': -46.5, 'E_B_RB': -48.0, 'E_A_RB': -49.5}}
        self.dataset = {'schema': 'wren-four-point-input-1', 'energy_kind': 'total',
                        'energy_unit': 'eV', 'method': 'analytic test surfaces',
                        'environment': 'frozen point charges', 'pairs': [self.pair]}

    def read(self, data):
        self.source.write_text(json.dumps(data))
        return read_external_four_points(self.source)

    def test_same_charge_different_root_and_frozen_environment(self):
        result = self.read(self.dataset)
        row = result['pairs'][0]
        self.assertEqual(row['relaxation_ev'], {'on_a': 0.5, 'on_b': 1.5, 'sum': 2.0})
        self.assertEqual(row['gaps_ev'], {'vertical_at_ra': 3.5, 'vertical_at_rb': 1.5, 'adiabatic': 2.0})
        self.assertEqual(row['energies'], self.pair['energies'])
        self.assertEqual(row['states'], self.pair['states'])
        self.assertEqual(result['environment'], self.dataset['environment'])
        self.assertFalse(row['diagnostics']['state_identity_verified'])
        self.assertFalse(row['diagnostics']['minima_verified'])
        self.assertFalse(row['diagnostics']['environment_consistency_verified'])
        self.assertEqual(row['declared_convergence'], {})
        self.assertEqual(json.loads(json.dumps(result, allow_nan=False)), result)

    def test_archived_carrier_results_reproduced_in_ev(self):
        result = read_external_four_points(EXAMPLE)
        self.assertEqual(len(result['pairs']), 2)
        archive = ROOT / 'runs/cube_input_2026-11-07/cpu-F2-reorganization.json'
        archived = json.loads(archive.read_text())['results']
        self.assertEqual(result['source_sha256'][str(archive)], hashlib.sha256(archive.read_bytes()).hexdigest())
        for imported, carrier in zip(result['pairs'], ('hole', 'electron')):
            reference = archived[carrier]['relaxation_ev']
            self.assertEqual(imported['energy_unit'], 'hartree')
            self.assertEqual(imported['relaxation_ev']['on_a'], reference['neutral_relaxation'])
            self.assertEqual(imported['relaxation_ev']['on_b'], reference['charged_relaxation'])
            self.assertEqual(imported['relaxation_ev']['sum'], reference['total'])

    def test_sources_hashed_and_convergence_declared(self):
        log = self.source.with_name('external.out')
        log.write_text('external total energy\n')
        self.pair.update(sources={'E_A_RA': 'external.out'}, converged={'E_A_RA': True})
        result = self.read(self.dataset)
        row = result['pairs'][0]
        self.assertEqual(row['sources'], {'E_A_RA': str(log)})
        self.assertTrue(row['declared_convergence']['E_A_RA'])
        self.assertEqual(result['source_sha256'][str(log)], hashlib.sha256(log.read_bytes()).hexdigest())
        log.write_text('changed output\n')
        self.assertNotEqual(read_external_four_points(self.source)['source_sha256'], result['source_sha256'])
        log.unlink()
        with self.assertRaisesRegex(ValueError, 'source does not exist'):
            read_external_four_points(self.source)

    def test_signed_relaxations_and_state_gap(self):
        self.pair['energies']['E_A_RB'] = -50.25
        row = self.read(self.dataset)['pairs'][0]
        self.assertEqual(row['relaxation_ev']['on_a'], -0.25)
        self.assertTrue(row['diagnostics']['negative_relaxation_a'])
        self.assertEqual(row['gaps_ev']['vertical_at_ra'] - row['gaps_ev']['vertical_at_rb'],
                         row['relaxation_ev']['sum'])

    def test_invalid_input_is_not_silently_used(self):
        for replacement in [
            {'energy_unit': 'mixed'}, {'energy_kind': 'excitation'}, {'pairs': []},
            {'pairs': [self.pair, copy.deepcopy(self.pair)]}, {'energy_unit': None},
            {'method': ''}, {'environment': None}, {'unexpected': 5},
            {'pairs': [{**self.pair, 'converged': {'E_B_RB': False}}]},
            {'pairs': [{**self.pair, 'converged': {'E_B_RB': 1}}]},
            {'pairs': [{**self.pair, 'energies': {**self.pair['energies'], 'E_B_RA': True}}]},
            {'pairs': [{**self.pair, 'energies': {**self.pair['energies'], 'E_B_RA': float('nan')}}]},
        ]:
            with self.subTest(replacement=replacement), self.assertRaises(ValueError):
                self.read({**self.dataset, **replacement})
        self.source.write_text('{"energy_unit":"hartree", "energy_unit":"eV"}')
        with self.assertRaisesRegex(ValueError, 'duplicate JSON key'):
            read_external_four_points(self.source)

    def test_no_optional_science_packages_required(self):
        script = ('import sys; from wren.workflow.external_reorganization import read_external_four_points; '
                  'read_external_four_points(sys.argv[1]); '
                  'assert not set(sys.modules).intersection(("numpy", "pyscf", "MDAnalysis", "cupy"))')
        subprocess.run([sys.executable, '-c', script, str(EXAMPLE)], check=True)
