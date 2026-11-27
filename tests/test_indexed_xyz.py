"""Index addressing and conversion must preserve the numeric geometry and atom order."""
import io
from contextlib import redirect_stdout, redirect_stderr
from pathlib import Path
import tempfile
import unittest
from zipfile import ZipFile

from wren.quantum.xyz2gjf import main as convert
from wren.workflow.frame_geometry import main as retrieve


class IndexedXYZTests(unittest.TestCase):
    def test_lookup_and_convert_plain_or_archived_xyz(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            frame = root/'frame0020'
            frame.mkdir()
            xyz = '2\nframe20.mol19; Angstrom\nH 1.000000000000001 0 0\nF 1.91 -1 0\n'
            with ZipFile(frame/'monomers.zip', 'w') as archive:
                archive.writestr('frame20.mol19.xyz', xyz)
            text = io.StringIO()
            with redirect_stdout(text):
                retrieve(['frame20.mol19', '--root', str(root)])
            self.assertEqual(text.getvalue(), xyz)
            loose, gjf = root/'sample.xyz', root/'sample.gjf'
            retrieve(['frame20.mol19', '--root', str(root), '-o', str(loose)])
            flags = ['--method','HF','--basis','STO-3G','--charge','-1','--multiplicity','2']
            convert([str(loose), str(gjf), *flags])
            expected = gjf.read_text()
            self.assertIn('HF/STO-3G', expected)
            coordinates = expected.split('-1 2\n')[1].splitlines()
            self.assertEqual([line.split()[0] for line in coordinates if line.strip()], ['H', 'F'])
            self.assertEqual(float(coordinates[0].split()[1]), 1.000000000000001)
            convert(['frame0020.mol19', str(gjf), '--root', str(root), *flags])
            self.assertEqual(gjf.read_text(), expected)
            for key in ['frame20.mol1', 'frame21.mol19', '../../frame20.mol19']:
                with self.subTest(key=key), self.assertRaises(SystemExit), redirect_stderr(io.StringIO()):
                    retrieve([key, '--root', str(root), '-o', str(root/'missing.xyz')])
                self.assertFalse((root/'missing.xyz').exists())
