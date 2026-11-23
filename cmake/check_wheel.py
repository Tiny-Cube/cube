"""Run the wheel's actual console entry point without importing the checkout."""
from importlib.metadata import distribution
from pathlib import Path
import subprocess
import sys
import tempfile
import zipfile


def launch(extracted):
    sys.path.insert(0, str(extracted))
    dist = distribution('wren')
    if Path(dist.locate_file('')).resolve() != extracted.resolve():
        raise RuntimeError('imported Wren package metadata outside the staged wheel')
    sys.argv = ['wren', '--help']
    entry = dist.entry_points.select(group='console_scripts', name='wren')
    return next(iter(entry)).load()()


if __name__ == '__main__':
    if sys.argv[1] == '--extracted':
        raise SystemExit(launch(Path(sys.argv[2])))
    with tempfile.TemporaryDirectory() as root:
        with zipfile.ZipFile(sys.argv[1]) as wheel:
            wheel.extractall(root)
        subprocess.run([sys.executable, '-I', '-S', str(Path(__file__).resolve()),
                        '--extracted', root], cwd=root, check=True)
