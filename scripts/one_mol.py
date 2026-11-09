#!/usr/bin/env python3
"""Run $one-mol from input; extract prerequisites if requested by the control."""
import sys
from universal_kmc.workflow.cli import main

if __name__ == '__main__':
    raise SystemExit(main([*sys.argv[1:], '--only', 'one-mol']))
