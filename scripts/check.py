#!/usr/bin/env python3
"""Build Wren with CMake and run its registered CTest suite."""
import subprocess

from build import BUILD, main as build


if __name__ == '__main__':
    build()
    subprocess.run(['cmake', '--build', str(BUILD), '--config', 'Release', '--target', 'wren_package'], check=True)
    subprocess.run(['ctest', '--test-dir', str(BUILD), '-C', 'Release', '--output-on-failure'], check=True)
