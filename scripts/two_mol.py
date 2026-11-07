#!/usr/bin/env python3
"""Run $two-mol, its single-point prerequisites and transfer-integral analysis."""
import sys
from universal_kmc.workflow.cli import main

if __name__ == '__main__':
    raise SystemExit(main([*sys.argv[1:], '--only', 'two-mol']))
