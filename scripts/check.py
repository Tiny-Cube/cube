#!/usr/bin/env python3
"""Full build, native tests, statistical tests, conformance and uv packaging."""
import os
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault('UV_CACHE_DIR', str(Path(tempfile.gettempdir()) / 'universal-kmc-uv-cache'))
os.environ.setdefault('UV_LINK_MODE', 'copy')


def run(*command):
    print('+', ' '.join(command), flush=True)
    subprocess.run(command, cwd=ROOT, check=True)


if __name__ == '__main__':
    run(sys.executable, 'scripts/build.py')
    run('cargo', 'test', '--manifest-path', 'rust/Cargo.toml', '--offline')
    run(str(ROOT / 'build/test-cpp'))
    run(str(ROOT / 'build/test-fortran'))
    run('uv', 'run', '--locked', 'python', '-m', 'unittest', 'discover', '-s', 'tests', '-v')
    run('uv', 'build')
