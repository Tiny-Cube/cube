#!/usr/bin/env python3
"""Run vacuum four-point $reorganization-energy from an input control."""
import sys
from wren.workflow.cli import main

if __name__ == '__main__':
    raise SystemExit(main([*sys.argv[1:], '--only', 'reorganization-energy']))
