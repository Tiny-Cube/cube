#!/usr/bin/env python3
"""Compatibility entry point for the project's CMake build."""
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
BUILD = Path(os.environ.get('WREN_BUILD_DIR', ROOT/'build/cmake')).resolve()


def main():
    subprocess.run(['cmake', '-S', str(ROOT), '-B', str(BUILD),
                    '-DPython3_EXECUTABLE=' + sys.executable], check=True)
    subprocess.run(['cmake', '--build', str(BUILD), '--config', 'Release', '--parallel'], check=True)


if __name__ == '__main__':
    main()
