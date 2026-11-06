"""Command-line entry point for reproducible fragment-coupling calculations."""
import argparse
import json
import os
from pathlib import Path
import sys
import tempfile


def _write_result(path, result):
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(mode='w', dir=path.parent, prefix=path.name+'.', delete=False) as file:
        temporary = Path(file.name)
        try:
            json.dump(result, file, indent=2, allow_nan=False)
            file.write('\n')
        except BaseException:
            temporary.unlink(missing_ok=True)
            raise
    try:
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fragment-a', required=True, type=Path)
    parser.add_argument('--fragment-b', required=True, type=Path)
    parser.add_argument('--backend', choices=('cpu', 'gpu'), default='cpu')
    parser.add_argument('--basis', required=True)
    parser.add_argument('--xc', required=True, help='DFT functional or hf')
    parser.add_argument('--method', action='append', choices=('dipro', 'frozen_density', 'esid', 'orbital_gmh'))
    parser.add_argument('--carrier', action='append', choices=('hole', 'electron'))
    parser.add_argument('--assume-equivalent', action='store_true', help='explicit ESID symmetry assumption, still subject to numerical checks')
    parser.add_argument('--unit', choices=('Angstrom', 'Bohr'), default='Angstrom')
    parser.add_argument('--threads', type=int, default=4)
    parser.add_argument('--grid-level', type=int, default=3)
    parser.add_argument('--max-cycle', type=int, default=100)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args(argv)
    try:
        from .transfer import SCFConfig, calculate_transfer
        config = SCFConfig(backend=args.backend, basis=args.basis, xc=args.xc,
                           unit=args.unit, threads=args.threads, grid_level=args.grid_level,
                           max_cycle=args.max_cycle)
        result = calculate_transfer(args.fragment_a, args.fragment_b, config=config,
                                    methods=args.method or ('dipro',), carriers=args.carrier or ('hole', 'electron'),
                                    assume_equivalent=args.assume_equivalent)
        _write_result(args.output, result)
        rejected = any(r['status'] != 'ok' for c in result['results'].values() for r in c['methods'].values())
        if rejected:
            print('One or more methods were not applicable; inspect output JSON.', file=sys.stderr)
            return 2
        return 0
    except (ImportError, ValueError, RuntimeError, OSError) as exc:
        # Replace stale successes with an explicit failure if this run fails.
        try:
            _write_result(args.output, {'schema': 'cube-transfer-1', 'status': 'failed',
                                        'error': {'type': type(exc).__name__, 'message': str(exc)}})
        except OSError:
            pass
        print(f'Transfer calculation failed: {exc}', file=sys.stderr)
        return 1


if __name__ == '__main__':
    raise SystemExit(main())
