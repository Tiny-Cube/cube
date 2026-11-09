#!/usr/bin/env python3
"""Stage 2: transfer.toml referencing A/B/AB .chk -> coupling JSON. No SCF.

uv run --extra quantum python scripts/calculate_transfer.py transfer.toml
Equivalent installed command: cube-transfer transfer.toml
"""
from universal_kmc.quantum.cli import main


if __name__ == '__main__':
    raise SystemExit(main())
