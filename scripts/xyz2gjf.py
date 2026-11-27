#!/usr/bin/env python3
"""Convert a single XYZ or indexed monomer to GJF from a source checkout."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from wren.quantum.xyz2gjf import main

if __name__ == '__main__':
    main()
