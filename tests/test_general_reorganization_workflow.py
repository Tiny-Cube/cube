"""Workflow-level checks for state-pair selection and two-ensemble semantics."""
import json
import math
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
import unittest

from wren.workflow.control import load_input
from wren.workflow.gap_reorganization import KB_EV_K, analyze_gap_ensembles
from wren.workflow.reorganization_task import _environment_cutoff
from wren.workflow.runner import run_workflow

ROOT = Path(__file__).resolve().parents[1]
EXAMPLES = ROOT/'examples/reorganization'


class GeneralReorganizationTests(unittest.TestCase):
    def test_linear_response_sign_equal_width_and_block_uncertainty(self):
        with TemporaryDirectory() as tmp:
            p = Path(tmp)
            (p/'input').write_text((EXAMPLES/'linear-response.input').read_text().replace(
                '../../build/reorganization/gap-ensembles', 'results'))
            def ensemble(key, values):
                (p/f'ensemble-{key}.csv').write_text('frame,E_A_eV,E_B_eV\n' +
                    ''.join(f'{i},0,{gap}\n' for i, gap in enumerate(values)))
            ensemble('a', [.5, .5, 1., 1.])
            ensemble('b', [-.5, -.5, 0., 0.])
            control = load_input(p/'input')
            self.assertEqual(control.plan()['tasks'], ['reorganization-energy'])
            run_workflow(control)
            path = p/'results/reorganization.json'
            result = json.loads(path.read_text())
            self.assertAlmostEqual(result['lambda_stokes_ev'], .5)
            self.assertAlmostEqual(result['delta_g_lr_ev'], .25)
            self.assertAlmostEqual(result['lambda_stokes_block_se_ev'], .25 / math.sqrt(2))
            self.assertAlmostEqual(result['ensembles']['a']['lambda_variance_ev'], 1/12/(2*KB_EV_K*300))
            self.assertFalse(result['diagnostics']['linear_response_verified'])
            run_workflow(control, resume=True)
            saved = path.read_bytes()
            ensemble('b', [1., 1., 1.5, 1.5])
            with self.assertRaisesRegex(ValueError, 'existing analysis'):
                run_workflow(control, resume=True)
            self.assertEqual(path.read_bytes(), saved)
            negative = analyze_gap_ensembles(control.sections['reorganization-energy'])
            self.assertTrue(negative['diagnostics']['negative_lambda_stokes'])

    def test_external_analysis_runs_without_scientific_packages(self):
        with TemporaryDirectory() as tmp:
            p = Path(tmp)
            (p/'input').write_text((EXAMPLES/'external.input').read_text().replace(
                'external-points.json', str(EXAMPLES/'external-points.json')).replace(
                '../../build/reorganization/external', 'results'))
            code = ("import sys; sys.path.insert(0, sys.argv[1]); from wren.workflow.cli import main; "
                    "assert main([sys.argv[2]]) == 0; assert not {'numpy','pyscf','cupy'} & sys.modules.keys()")
            r = subprocess.run([sys.executable, '-S', '-c', code, str(ROOT/'src'), str(p/'input')],
                               capture_output=True, text=True)
            self.assertEqual(r.returncode, 0, r.stderr)
            pairs = json.loads((p/'results/reorganization.json').read_text())['pairs']
            self.assertEqual(len(pairs), 2)

    def test_automatic_capabilities_and_conflicts_before_qc(self):
        with TemporaryDirectory() as tmp:
            p = Path(tmp)/'input'
            valid = (EXAMPLES/'singlet.input').read_text().replace('../dos/h2.xyz', str(ROOT/'examples/dos/h2.xyz'))
            for extra, message in [('qc-engine=gpu-pyscf', 'TDA'), ('environment=qmmm','TDA'),
                                   ('carriers=hole', 'either'), ('root-min-overlap=1.1', 'overlap')]:
                k = extra.split('=')[0]
                lines = [line for line in valid.splitlines() if not line.startswith(k + ' =')]
                p.write_text('\n'.join(lines).replace('$end', extra + '\n$end'))
                with self.assertRaisesRegex(ValueError, message):
                    load_input(p)
            p.write_text(valid.replace('state-b = S1','state-b = S0'))
            with self.assertRaisesRegex(ValueError, 'differ'):
                load_input(p)
            p.write_text(valid)
            self.assertEqual(load_input(p).sections['reorganization-energy']['adapter'], 'tda')

    def test_frozen_cluster_uses_extraction_manifest_distance(self):
        data = {'provenance': {'settings': {'monomers-qmmm-cutoff': 6.}}}
        self.assertEqual(_environment_cutoff({}, data, {}), 6.)
        self.assertEqual(_environment_cutoff({}, {}, {'cluster': {'cutoff_angstrom': 7.}}), 7.)
        self.assertEqual(_environment_cutoff({'mm-cutoff': 5}, data, {}), 5)
        self.assertEqual(_environment_cutoff({}, {}, {}), 12.)
