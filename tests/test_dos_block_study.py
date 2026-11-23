from contextlib import redirect_stdout
import csv
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


@unittest.skipUnless(importlib.util.find_spec('numpy'), 'DOS NumPy dependency required')
class OfflineDosBlocksTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.source = self.root/'results'
        self.source.mkdir()
        self.control = self.root/'dos.toml'
        self.control.write_text('[dos]\ninput="removed-md-files/frames.json"\noutput="results"\n')

    def archive(self, *, mode='monomer-orbitals', indices=range(10, 24, 2), grid=(-1., 0., 1.)):
        import numpy as np
        # Standard Wren output columns frame_0..frame_6 preserve the order in
        # spectra.json, rather than encoding the actual MD frame numbers.
        curves = np.array([[0, 2, 5], [2, 2, 7], [3, 2, 0], [5, 2, 0],
                           [6, 2, 0], [8, 2, 0], [100, 2, 100]])
        with (self.source/'dos.csv').open('w', newline='') as stream:
            writer = csv.writer(stream)
            writer.writerow(['energy_ev', 'dos', *[f'frame_{i}' for i in range(7)]])
            writer.writerows(zip(grid, curves.mean(axis=0), *curves))
        datasets = [{'frame': i, 'time_ps': i*.25} for i in indices]
        frames = [dict(d) for d in datasets]
        settings = {'mode': mode}
        if mode == 'periodic':
            settings['periodic-reference'] = 'external'
            for frame in frames:
                frame['energy_alignment'] = {'mode': 'external', 'reference_ev': .5}
        request = {'settings': settings, 'datasets': datasets}
        (self.source/'spectra.json').write_text(json.dumps({'schema': 'wren-dos-spectra-1', 'frames': frames}))
        self.result = {'schema': 'wren-dos-1', 'status': 'completed', 'mode': mode, 'request': request,
            'fingerprint': hashlib.sha256(json.dumps(request, sort_keys=True).encode()).hexdigest(),
            'frames': 'spectra.json', 'frame_count': 7, 'energy_unit': 'eV', 'energy_reference': 'test_ref',
            'normalization': 'states_per_model_per_eV', 'grid_spacing_ev': 1., 'broadening': {'sigma_ev': .1}}
        self.seal()

    def seal(self):
        self.result['output_sha256'] = {name: hashlib.sha256((self.source/name).read_bytes()).hexdigest()
                                        for name in ('spectra.json', 'dos.csv')}
        (self.source/'result.json').write_text(json.dumps(self.result))

    def test_normal_cli_reuses_only_the_archived_spectra_and_extends_sizes(self):
        import numpy as np
        from wren.dos.block_study import run_block_study
        from wren.dos.cli import main
        self.archive()
        output = self.root/'results-blocks'
        script = ('import sys; from wren.dos.cli import main; '
                  'assert main(sys.argv[1:])==0; assert "wren.dos.runner" not in sys.modules; '
                  'assert not any(k.startswith(("pyscf","gpu4pyscf")) for k in sys.modules)')
        process = subprocess.run([sys.executable, '-c', script, str(self.control), '--block-study', '1', '2', '--check'],
                                 cwd=self.root, text=True, capture_output=True)
        self.assertEqual(process.returncode, 0, process.stderr)
        self.assertFalse(output.exists())
        with redirect_stdout(io.StringIO()):
            self.assertEqual(main([str(self.control), '--block-study', '1', '2']), 0)
        report = json.loads((output/'block-study.json').read_text())
        self.assertEqual(report['frame_indices'], list(range(10, 24, 2)))
        self.assertEqual(report['studies'][1]['discarded_tail_frames'], 1)
        csv2 = np.genfromtxt(output/'block-2.csv', names=True, delimiter=',')
        np.testing.assert_allclose(csv2['retained_mean'], [4, 2, 2])
        np.testing.assert_allclose(csv2['standard_error'], [np.sqrt(3), 0, 2])
        self.assertTrue(run_block_study(self.source, [1, 2], resume=True)['reused'])
        self.assertEqual(run_block_study(self.source, [2, 3], resume=True)['block_sizes'], [2, 3])
        self.assertTrue((output/'block-1.csv').is_file())
        (output/'block-3.csv').write_text('edited')
        with self.assertRaisesRegex(ValueError, 'changed'):
            run_block_study(self.source, [2, 3], resume=True)

    def test_source_integrity_grid_and_frame_order_checked_before_outputs(self):
        from wren.dos.block_study import run_block_study
        for problem, pattern in (('request', 'fingerprint'), ('status', 'completed'), ('hash', 'hash'),
                                 ('grid', 'energy axis'), ('mean', 'mean DOS'), ('labels', 'ordered frame'),
                                 ('frame-order', 'source request order'), ('gapped', 'equally spaced')):
            with self.subTest(problem=problem):
                self.archive(indices=[10, 12, 14, 18, 20, 22, 24] if problem == 'gapped' else range(10, 24, 2),
                             grid=[-1, .2, 1.] if problem == 'grid' else [-1., 0., 1.])
                if problem in ('hash', 'labels', 'mean'):
                    csv_path = self.source/'dos.csv'
                    lines = csv_path.read_text().splitlines()
                    if problem == 'labels':
                        lines[0] = lines[0].replace('frame_1', 'frame_0')
                    else:
                        row = lines[1].split(',')
                        row[1] = '123.0'
                        lines[1] = ','.join(row)
                    csv_path.write_text('\n'.join(lines)+'\n')
                if problem == 'request':
                    self.result['request']['settings']['sigma'] = .5
                if problem == 'status':
                    self.result['status'] = 'failed'
                if problem == 'frame-order':
                    p = self.source/'spectra.json'
                    spectra = json.loads(p.read_text())
                    spectra['frames'].reverse()
                    p.write_text(json.dumps(spectra))
                if problem != 'hash':
                    self.seal()
                with self.assertRaisesRegex(ValueError, pattern):
                    run_block_study(self.source, [1, 2])
                self.assertFalse((self.root/'results-blocks').exists())

    def test_aligned_periodic_and_preflight_block_sizes(self):
        from wren.dos.block_study import run_block_study
        self.archive(mode='periodic')
        result = run_block_study(self.source, [3], expected_mode='periodic', check_only=True)
        self.assertEqual(result['status'], 'validated')
        for sizes in ([1, 1], [0], [4], [True], []):
            with self.subTest(sizes=sizes), self.assertRaises(ValueError):
                run_block_study(self.source, sizes)
        with self.assertRaisesRegex(ValueError, 'matching'):
            run_block_study(self.source, [2], expected_mode='monomer-orbitals')
        p = self.source/'spectra.json'
        spectra = json.loads(p.read_text())
        del spectra['frames'][0]['energy_alignment']
        p.write_text(json.dumps(spectra))
        self.seal()
        with self.assertRaisesRegex(ValueError, 'aligned'):
            run_block_study(self.source, [2])
        self.assertFalse((self.root/'results-blocks').exists())
