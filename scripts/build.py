#!/usr/bin/env python3
"""Build dependency-free native libraries, CLIs, tests and generic examples."""
import os
from pathlib import Path
import shlex
import subprocess

ROOT = Path(__file__).resolve().parents[1]
BUILD = ROOT / 'build'


def run(command):
    print('+', shlex.join(map(str, command)), flush=True)
    subprocess.run(list(map(str, command)), cwd=ROOT, check=True)


def main():
    (BUILD / 'mod').mkdir(parents=True, exist_ok=True)
    cxx = shlex.split(os.environ.get('CXX', 'c++'))
    fc = shlex.split(os.environ.get('FC', 'gfortran'))
    cflags = ['-std=c++17', '-O2', '-Wall', '-Wextra', '-Wpedantic', '-Werror', '-Icpp/include']
    for source, output in [('cpp/src/main.cpp', 'ukmc-cpp'), ('cpp/tests/test.cpp', 'test-cpp'), ('examples/generic.cpp', 'generic-cpp')]:
        run([*cxx, *cflags, source, '-o', BUILD / output])
    fflags = ['-std=f2008', '-O2', '-g', '-Wall', '-Wextra', '-Wno-compare-reals', '-Wno-unused-dummy-argument', '-fcheck=all', '-fbacktrace', '-Jbuild/mod', '-Ibuild/mod']
    objects = []
    for module in ['ukmc', 'ukmc_io', 'ukmc_events']:
        obj = BUILD / f'{module}.o'
        run([*fc, *fflags, '-c', f'fortran/{module}.f90', '-o', obj])
        objects.append(obj)
    run(['ar', 'rcs', BUILD / 'libukmc_fortran.a', *objects])
    for source, output in [('fortran/main.f90', 'ukmc-fortran'), ('fortran/test.f90', 'test-fortran'), ('examples/generic.f90', 'generic-fortran')]:
        run([*fc, *fflags, source, BUILD / 'libukmc_fortran.a', '-o', BUILD / output])
    run(['cargo', 'build', '--manifest-path', 'rust/Cargo.toml', '--release', '--examples', '--offline'])
    run(['cargo', 'build', '--manifest-path', 'rust/Cargo.toml', '--release', '--offline'])


if __name__ == '__main__':
    main()
