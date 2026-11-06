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
    parser.add_argument('control', nargs='?', type=Path, help='TOML control referencing three existing PySCF checkpoints')
    parser.add_argument('--check', action='store_true', help='validate control paths/settings only; do not read matrices or run SCF')
    parser.add_argument('--fragment-a', type=Path)
    parser.add_argument('--fragment-b', type=Path)
    parser.add_argument('--backend', choices=('cpu', 'gpu'))
    parser.add_argument('--basis')
    parser.add_argument('--xc', help='DFT functional or hf')
    parser.add_argument('--method', action='append', choices=('dipro', 'frozen_density', 'esid', 'orbital_gmh'))
    parser.add_argument('--carrier', action='append', choices=('hole', 'electron'))
    parser.add_argument('--assume-equivalent', action='store_true', default=None, help='explicit ESID symmetry assumption, still subject to numerical checks')
    parser.add_argument('--unit', choices=('Angstrom', 'Bohr'))
    parser.add_argument('--threads', type=int)
    parser.add_argument('--grid-level', type=int)
    parser.add_argument('--max-cycle', type=int)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args(argv)
    overrides = {k: v for k, v in vars(args).items() if k not in ('control', 'check') and v is not None}
    if args.control is not None and overrides:
        parser.error('a control file cannot be combined with calculation flags; edit the control file')
    if args.control is None:
        if args.check:
            parser.error('--check requires a control file')
        for key in ('fragment_a', 'fragment_b', 'basis', 'xc', 'output'):
            if getattr(args, key) is None:
                parser.error(f'--{key.replace("_", "-")} is required without a control file')
    output, provenance = None, None
    try:
        from .control import load_control, validate_output
        if args.control is not None:
            job = load_control(args.control)
            if args.check:
                print(json.dumps(job.summary(), indent=2, allow_nan=False))
                return 0
            output, provenance = job.output, job.provenance
            from .checkpoints import calculate_from_checkpoints
            result = calculate_from_checkpoints(job)
        else:
            output = validate_output(args.output, [args.fragment_a, args.fragment_b])
            from .transfer import SCFConfig, calculate_transfer
            settings = {k: getattr(args, k) for k in ('backend', 'basis', 'xc', 'unit', 'threads', 'grid_level', 'max_cycle') if getattr(args, k) is not None}
            config = SCFConfig(**settings)
            result = calculate_transfer(args.fragment_a, args.fragment_b, config=config,
                                        methods=args.method or ('dipro',), carriers=args.carrier or ('hole', 'electron'),
                                        assume_equivalent=bool(args.assume_equivalent))
        _write_result(output, result)
        rejected = any(r['status'] != 'ok' for c in result['results'].values() for r in c['methods'].values())
        if rejected:
            print('One or more methods were not applicable; inspect output JSON.', file=sys.stderr)
            return 2
        return 0
    except (ImportError, ValueError, RuntimeError, OSError, KeyError, TypeError) as exc:
        # Only a validated destination may replace a stale result on failure.
        if output is not None:
            try:
                failure = {'schema': 'cube-transfer-1', 'status': 'failed',
                           'error': {'type': type(exc).__name__, 'message': str(exc)}}
                if provenance is not None:
                    failure['control'] = provenance
                _write_result(output, failure)
            except OSError:
                pass
        print(f'Transfer calculation failed: {exc}', file=sys.stderr)
        return 1


if __name__ == '__main__':
    raise SystemExit(main())
