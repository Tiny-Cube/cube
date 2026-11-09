"""Run model-Hamiltonian FSSH ensembles using a TOML control file."""
import argparse
import json
import sys
from .control import load_control


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('control')
    parser.add_argument('--output', help='override output directory')
    parser.add_argument('--check', action='store_true')
    parser.add_argument('--resume', action='store_true')
    args = parser.parse_args(argv)
    try:
        config, directory, provenance, fingerprint = load_control(args.control, args.output)
        if args.check:
            print(json.dumps({'output': str(directory), 'provenance': provenance, 'fingerprint': fingerprint}, indent=2))
            return 0
        from .runner import run_control
        result = run_control(args.control, output=directory, resume=args.resume)
        print(f'cube: {result["status"]}; FSSH result: {result["output"]}')
        return 0 if result['status'] == 'completed' else 2
    except (ImportError, ValueError, RuntimeError, OSError, KeyError, TypeError) as exc:
        print(f'cube-surface-hopping: {type(exc).__name__}: {exc}', file=sys.stderr)
        return 1


if __name__ == '__main__':
    raise SystemExit(main())
