import csv
from contextlib import redirect_stdout
from copy import deepcopy
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from universal_kmc.workflow.reorganization_csv import COLUMNS, write_reorganization_csv


def completed_result(shift=0., carriers=('hole', 'electron')):
    # Nonzero absolute energies make an accidental Hartree/eV conversion visible.
    return {'schema': 'cube-reorganization-1', 'status': 'completed', 'reused': False,
            'results': {carrier: {
                'energies_hartree': {'E0_R0': -10.+shift, 'Eq_R0': -9.+shift,
                                     'Eq_Rq': -9.2+shift, 'E0_Rq': -9.8+shift},
                'relaxation_ev': {'charged_relaxation': 5.4422772491976,
                                  'neutral_relaxation': 5.4422772491976,
                                  'total': 10.8845544983952},
                'reorganization_energy_ev': 10.8845544983952} for carrier in carriers}}


class ReorganizationCsvTests(unittest.TestCase):
    def test_exact_header_units_ids_and_replacement(self):
        with tempfile.TemporaryDirectory() as directory:
            results = {'12': completed_result(), '2': completed_result(.5)}
            paths = write_reorganization_csv(results, directory, ('hole', 'electron'))
            for carrier, name in paths.items():
                with Path(name).open(newline='') as file:
                    reader = csv.DictReader(file)
                    self.assertEqual(reader.fieldnames, list(COLUMNS))
                    rows = list(reader)
                self.assertEqual([r['mol'] for r in rows], ['12', '2'])
                for row in rows:
                    block = results[row['mol']]['results'][carrier]
                    for key in COLUMNS[1:5]:
                        self.assertEqual(float(row[key]), block['energies_hartree'][key])
                    for key in COLUMNS[5:]:
                        self.assertEqual(float(row[key]), block['relaxation_ev'][key])
            write_reorganization_csv({'2': results['2']}, directory, ('hole', 'electron'))
            with Path(paths['hole']).open(newline='') as file:
                self.assertEqual([r['mol'] for r in csv.DictReader(file)], ['2'])

    def test_empty_and_single_channel(self):
        with tempfile.TemporaryDirectory() as directory:
            paths = write_reorganization_csv({}, directory, ('hole',))
            self.assertEqual(set(paths), {'hole'})
            self.assertFalse((Path(directory)/'electron.csv').exists())
            self.assertEqual(Path(paths['hole']).read_text().strip(), ','.join(COLUMNS))
            result = completed_result(carriers=('hole',))
            write_reorganization_csv({'molecule_A': result}, directory, ('hole',))
            with Path(paths['hole']).open(newline='') as file:
                self.assertEqual(next(csv.DictReader(file))['mol'], 'molecule_A')

    def test_failed_missing_or_nonfinite_data_never_replace_valid_tables(self):
        with tempfile.TemporaryDirectory() as directory:
            good = completed_result()
            paths = write_reorganization_csv({'1': good}, directory, ('hole', 'electron'))
            before = {c: Path(p).read_bytes() for c, p in paths.items()}
            bad_cases = []
            failed = deepcopy(good)
            failed['status'] = 'failed'
            bad_cases.append(failed)
            missing = deepcopy(good)
            del missing['results']['electron']['energies_hartree']['Eq_Rq']
            bad_cases.append(missing)
            invalid = deepcopy(good)
            invalid['results']['hole']['relaxation_ev']['total'] = float('nan')
            bad_cases.append(invalid)
            for result in bad_cases:
                with self.assertRaises(ValueError):
                    write_reorganization_csv({'1': result}, directory, ('hole', 'electron'))
                self.assertEqual({c: Path(p).read_bytes() for c, p in paths.items()}, before)


try:
    import numpy
except ImportError:
    numpy = None


@unittest.skipIf(numpy is None, 'workflow numpy dependency required')
class ReorganizationCsvWorkflowTests(unittest.TestCase):
    def test_multiframe_selected_molecules_and_cached_results_export(self):
        from universal_kmc.workflow.control import load_input
        from universal_kmc.workflow.dimers import file_hash
        from universal_kmc.workflow.runner import run_workflow
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            entries = []
            for frame in (0, 20):
                here = root/f'frame{frame:04d}'
                here.mkdir()
                manifest = here/'manifest.json'
                manifest.write_text(json.dumps({
                    'schema': 'cube-dimers-1', 'coordinate_unit': 'Angstrom', 'frame': frame,
                    'monomers': [{'id': str(i), 'elements': ['He'],
                                  'coordinates_angstrom': [[float(i), 0., 0.]]} for i in (1, 2, 3)],
                    'dimers': []}))
                entries.append({'frame': frame, 'manifest': str(manifest.relative_to(root)),
                                'sha256': file_hash(manifest)})
            (root/'frames.json').write_text(json.dumps({
                'schema': 'cube-frames-1', 'frame_count': 2, 'frames': entries}))
            control = root/'input'
            control.write_text('$reorganization-energy\ndo = 1\ninput = frames.json\n'
                               'environment = vacuum\ncarriers = hole,electron\n'
                               'molecules = 2,1\noutput = results/reorganization\n$end\n')

            def calculate(atoms, settings, output, *, resume, embedding):
                result = completed_result(atoms[0][1][0])
                result['reused'] = resume
                return result

            with patch('universal_kmc.workflow.reorganization.calculate_reorganization',
                       side_effect=calculate) as mocked, redirect_stdout(io.StringIO()):
                self.assertEqual(run_workflow(load_input(control)), 0)
                self.assertEqual(mocked.call_count, 4)
                for label in ('frame0000', 'frame0020'):
                    here = root/'results/reorganization'/label
                    for carrier in ('hole', 'electron'):
                        with (here/f'{carrier}.csv').open(newline='') as file:
                            rows = list(csv.DictReader(file))
                        self.assertEqual([r['mol'] for r in rows], ['2', '1'])
                        self.assertEqual(float(rows[0]['E0_R0']), -8.)
                        (here/f'{carrier}.csv').unlink()
                self.assertEqual(run_workflow(load_input(control), resume=True), 0)
            report = json.loads((root/'cube-results/workflow.json').read_text())
            for frame in report['tasks']['reorganization-energy']['frames'].values():
                self.assertEqual(frame['count'], 2)
                self.assertEqual(frame['csv_units']['E0_R0'], 'Hartree')
                self.assertEqual(frame['csv_units']['total'], 'eV')
                self.assertTrue(all(p['reused'] for p in frame['outputs'].values()))
                self.assertTrue(all(Path(p).is_file() for p in frame['csv'].values()))


if __name__ == '__main__':
    unittest.main()
