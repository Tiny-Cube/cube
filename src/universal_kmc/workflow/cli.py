"""cube input: interpret section controls and dispatch the requested tasks."""
import argparse
import json
import sys

from .control import load_input
from .calculations import section_names


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', nargs='?', default='input', help='section control file; default: input')
    parser.add_argument('--plan', '--check', action='store_true', dest='plan', help='validate controls and show task order without MD/QC/GPU execution')
    parser.add_argument('--only', action='append', choices=(*section_names(), 'get_dimer'), help='select main workflows or submodule subsets (repeatable); submodule do=1 runs independently; respects explicit do=0')
    parser.add_argument('--resume', action='store_true', help='reuse only identical, successful, hash-verified tasks')
    args = parser.parse_args(argv)
    try:
        control = load_input(args.input)
        plan = control.plan(args.only)
        if args.plan:
            print(json.dumps(plan, indent=2, allow_nan=False))
            return 0
        control.execution_plan(args.only)
        from .runner import run_workflow
        return run_workflow(control, only=args.only, resume=args.resume)
    except (ImportError, ValueError, RuntimeError, OSError, KeyError, TypeError) as exc:
        print(f'cube: {type(exc).__name__}: {exc}', file=sys.stderr)
        if isinstance(exc, ImportError):
            print('Install optional dependencies: uv sync --extra workflow [--extra quantum-cuda13]', file=sys.stderr)
        return 1


if __name__ == '__main__':
    raise SystemExit(main())
