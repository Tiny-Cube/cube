"""Run DOS from a [dos] TOML control, independently of transport workflows."""
import argparse
import json
from pathlib import Path
import sys
import tomllib

from .control import normalize


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('control', type=Path)
    parser.add_argument('--check', action='store_true')
    parser.add_argument('--resume', action='store_true')
    args = parser.parse_args(argv)
    try:
        source = args.control.resolve()
        raw = tomllib.loads(source.read_text())
        if set(raw) != {'dos'}:
            raise ValueError('DOS control requires exactly one [dos] table')
        settings = normalize(raw['dos'], source.parent)
        if args.check:
            print(json.dumps(settings, indent=2))
            return 0
        from .runner import run
        result = run(settings, resume=args.resume)
        print(f'wren: {result["status"]}; DOS result: {result["output"]}')
        return 0 if result['status'] == 'completed' else 2
    except (ValueError, TypeError, OSError, ImportError, KeyError, RuntimeError) as exc:
        print(f'wren-dos: {type(exc).__name__}: {exc}', file=sys.stderr)
        return 1


if __name__ == '__main__':
    raise SystemExit(main())
