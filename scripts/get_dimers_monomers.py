#!/usr/bin/env python3
"""Extract monomers, frozen-MM shells and dimers selected by $get_dimers_monomers in the supplied input control."""
import sys
from universal_kmc.workflow.cli import main

if __name__ == '__main__':
    raise SystemExit(main([*sys.argv[1:], '--only', 'get_dimers_monomers']))
