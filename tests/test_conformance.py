"""Native binaries must be built first: uv run python scripts/build.py."""
import csv
import io
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from universal_kmc import Engine
from universal_kmc.rates import marcus, arrhenius, miller_abrahams

ROOT = Path(__file__).resolve().parents[1]
COMMANDS = {
    'python': [sys.executable, '-m', 'universal_kmc'],
    'cpp': [str(ROOT/'build/ukmc-cpp')],
    'rust': [str(ROOT/'rust/target/release/ukmc-rust')],
    'fortran': [str(ROOT/'build/ukmc-fortran')],
}


def run(command, *args):
    return subprocess.run([*command, *map(str, args)], cwd=ROOT, text=True, capture_output=True, timeout=30)


class ConformanceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        missing = [command[0] for name, command in COMMANDS.items() if name != 'python' and not Path(command[0]).is_file()]
        if missing:
            raise RuntimeError('Build native binaries first with scripts/build.py: ' + ', '.join(missing))

    def compare(self, model, steps, seed, until=None):
        outputs = {}
        args = [model, steps, seed] + ([] if until is None else [until])
        for name, command in COMMANDS.items():
            result = run(command, *args)
            self.assertEqual(result.returncode, 0, f'{name}: {result.stderr}')
            outputs[name] = list(csv.DictReader(io.StringIO(result.stdout)))
            self.assertIn('status=', result.stderr)
        reference = outputs['python']
        for name, rows in outputs.items():
            self.assertEqual(len(rows), len(reference), name)
            for actual, expected in zip(rows, reference):
                self.assertEqual(set(actual), set(expected))
                for key in expected:
                    if key in ('time', 'x', 'y', 'z'):
                        self.assertTrue(math.isclose(float(actual[key]), float(expected[key]), rel_tol=3e-13, abs_tol=1e-25), (name, key, actual, expected))
                    else:
                        self.assertEqual(int(actual[key]), int(expected[key]), (name, key))
        return reference

    def test_graph_and_reaction_trajectories(self):
        for filename in ['two_state', 'periodic_walk', 'birth_death', 'dimerization', 'absorption']:
            for seed in [0, 42, 987654321]:
                with self.subTest(model=filename, seed=seed):
                    self.compare(ROOT/f'examples/{filename}.kmc', 300, seed)
                    rows = self.compare(ROOT/f'examples/{filename}.kmc', 10000, seed, 10.0)
                    self.assertEqual(float(rows[-1]['time']), 10.0)

    def test_zero_events_zero_steps_zero_horizon(self):
        with tempfile.TemporaryDirectory() as tmp:
            model=Path(tmp)/'model.kmc'
            model.write_text('KMC1 graph 1 0 0')
            self.assertEqual(len(self.compare(model, 5, 42)),1)
            self.assertEqual(len(self.compare(model, 5, 42, 10)),2)
            self.assertEqual(len(self.compare(model, 0, 42, 10)),1)
            self.assertEqual(len(self.compare(model, 5, 42, 0)),1)

    def test_invalid_input_rejected_in_all_languages(self):
        invalid = [
            'KMC1 graph 1 1 0 0 0 -1 0 0 0',
            'KMC1 graph 1 1 0 0 0 NaN 0 0 0',
            'KMC1 graph 1 1 0 0 0 1e309 0 0 0',
            'KMC1 graph 1 1 0 0 1 1 0 0 0',
            'KMC1 graph 1 0 0 extra',
            'KMC1 reaction 0 0', 'KMC1 reaction 1 1 1 1 2',
            'KMC1 graph 1.0 0 0', 'KMC1 graph +1 0 0',
            'KMC1 graph 1 1 0 0 0 1d0 0 0 0',
            'KMC1 graph 1 2 0 0 0 1e308 0 0 0 0 0 1e308 0 0 0',
            'KMC1 reaction 1 1 9007199254740991 1 0 1',
        ]
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'invalid.kmc'
            for text in invalid:
                path.write_text(text)
                for name, command in COMMANDS.items():
                    with self.subTest(text=text, language=name):
                        result=run(command,path,2,1)
                        self.assertNotEqual(result.returncode,0,result.stdout)

    def test_generic_queue_all_languages(self):
        import importlib.util
        spec = importlib.util.spec_from_file_location('generic_example', ROOT/'examples/generic.py')
        module = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = module
        spec.loader.exec_module(module)
        model=module.queue_model()
        engine=Engine(model,42)
        list(engine.run(10000,100))
        expected=[len(model.data.jobs),model.data.completed,engine.steps]
        for binary in ['build/generic-cpp','build/generic-fortran','rust/target/release/examples/generic']:
            result=run([str(ROOT/binary)])
            self.assertEqual(result.returncode,0,result.stderr)
            self.assertEqual([int(n) for n in result.stdout.split()],expected,binary)

    def test_native_rate_laws_match_python(self):
        expected=[marcus(.01,.2,.05),arrhenius(1e12,.1),miller_abrahams(1e12,1,2,.1)]
        for binary in ['build/test-cpp','build/test-fortran']:
            result=run([str(ROOT/binary)])
            self.assertEqual(result.returncode,0,result.stderr)
            actual=[float(x) for x in result.stdout.split()]
            for x,y in zip(actual,expected,strict=True):
                self.assertTrue(math.isclose(x,y,rel_tol=1e-13),(binary,x,y))


if __name__=='__main__':
    unittest.main()
