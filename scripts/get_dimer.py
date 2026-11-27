#!/usr/bin/env python3
"""Write indexed monomer XYZs and dimer/cluster membership tables."""
import sys
from wren.workflow.cli import main

if __name__ == '__main__':
    raise SystemExit(main([*sys.argv[1:], '--only', 'get_dimers_monomers']))
