"""Control-file tests run without NumPy, PySCF or a GPU."""
from contextlib import redirect_stderr, redirect_stdout
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from universal_kmc.quantum.cli import main
from universal_kmc.quantum.control import load_control, validate_output


CONTROL = '''schema_version = 1
[input]
fragment_a = "a.chk"
fragment_b = "b.chk"
dimer = "ab.chk"
[calculation]
backend = "cpu"
methods = ["dipro"]
carriers = ["hole", "electron"]
[fragments]
a_atoms = "1-2, 5"
b_atoms = [3, 4, 6]
[orbitals.hole]
a = [7, 8]
b = "7-8"
[output]
file = "nested/transfer.json"
'''


class TransferControlTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.directory = Path(self.temp.name)
        for name in ('a.chk', 'b.chk', 'ab.chk'):
            (self.directory/name).write_bytes(b'placeholder; --check must not read matrices')
        self.path = self.directory/'transfer.toml'
        self.path.write_text(CONTROL)

    def test_paths_indices_and_dependency_free_check_from_other_cwd(self):
        job = load_control(self.path)
        self.assertEqual(job.output, self.directory/'nested/transfer.json')
        self.assertEqual(job.atom_indices, ([0, 1, 4], [2, 3, 5]))
        self.assertEqual(job.orbital_indices, {'hole': ([6, 7], [6, 7])})
        self.assertFalse(job.assume_converged)
        code = '''import sys
from universal_kmc.quantum.cli import main
assert main([sys.argv[1], '--check']) == 0
assert 'numpy' not in sys.modules and 'pyscf' not in sys.modules and 'cupy' not in sys.modules
'''
        env = dict(os.environ, PYTHONPATH=str(Path(__file__).resolve().parents[1]/'src'))
        result = subprocess.run([sys.executable, '-c', code, str(self.path)], cwd='/tmp', env=env,
                                capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        summary = json.loads(result.stdout)
        self.assertEqual(summary['validation'], 'paths_and_settings_only_no_checkpoint_matrices_read')
        self.assertFalse(job.output.exists())

    def test_bad_settings_rejected_before_import_or_output(self):
        variants = [
            CONTROL.replace('schema_version = 1', 'schema_version = true'),
            CONTROL.replace('backend = "cpu"', 'backend = "typo"'),
            CONTROL.replace('backend = "cpu"', 'backend = "cpu"\nbasis = "def2-svp"'),
            CONTROL.replace('methods = ["dipro"]', 'methods = ["frozen_density"]'),
            CONTROL.replace('methods = ["dipro"]', 'methods = ["dipro", "dipro"]'),
            CONTROL.replace('backend = "cpu"', 'backend = "cpu"\nassume_converged = 1'),
            CONTROL.replace('a_atoms = "1-2, 5"', 'a_atoms = "0, 1"'),
            CONTROL.replace('a_atoms = "1-2, 5"', 'a_atoms = "2-1"'),
            CONTROL.replace('a_atoms = "1-2, 5"', 'a_atoms = "1, 1"'),
            CONTROL.replace('a_atoms = "1-2, 5"', 'a_atoms = "1-3"'),
            CONTROL.replace('a = [7, 8]', 'a = [true]'),
            CONTROL.replace('b = "7-8"', 'b = "7-1000000"'),
            CONTROL.replace('fragment_b = "b.chk"', 'fragment_b = "a.chk"'),
            CONTROL.replace('fragment_b = "b.chk"', 'fragment_b = "missing.chk"'),
            CONTROL.replace('schema_version = 1', 'schema_version = 1\ntypo = 1'),
        ]
        for content in variants:
            with self.subTest(content=content):
                self.path.write_text(content)
                with self.assertRaises(ValueError):
                    load_control(self.path)
        self.assertFalse((self.directory/'nested').exists())

    def test_output_cannot_overwrite_inputs_control_or_hardlink(self):
        for output in ('a.chk', 'transfer.toml'):
            self.path.write_text(CONTROL.replace('nested/transfer.json', output))
            with redirect_stderr(io.StringIO()):
                self.assertEqual(main([str(self.path)]), 1)
            self.assertEqual((self.directory/'a.chk').read_bytes(), b'placeholder; --check must not read matrices')
            self.assertTrue(self.path.read_text().startswith('schema_version'))
        alias = self.directory/'alias.chk'
        os.link(self.directory/'a.chk', alias)
        with self.assertRaises(ValueError):
            validate_output(alias, [self.directory/'a.chk'])

    def test_no_ambiguous_flags_and_check_does_not_write_output(self):
        with redirect_stderr(io.StringIO()), self.assertRaises(SystemExit) as context:
            main([str(self.path), '--backend', 'cpu'])
        self.assertEqual(context.exception.code, 2)
        with redirect_stdout(io.StringIO()):
            self.assertEqual(main([str(self.path), '--check']), 0)
        self.assertFalse((self.directory/'nested').exists())

    def test_plain_logs_are_not_accepted_as_complete_checkpoint(self):
        (self.directory/'a.log').write_text('converged SCF energy = -1')
        self.path.write_text(CONTROL.replace('a.chk', 'a.log'))
        with self.assertRaisesRegex(ValueError, 'text SCF logs'):
            load_control(self.path)


if __name__ == '__main__':
    unittest.main()
