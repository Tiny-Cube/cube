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
    study = parser.add_mutually_exclusive_group()
    study.add_argument('--kmesh-study', nargs='+', metavar='NxNxN',
                        help='compare successively refined periodic meshes, e.g. 1x1x1 2x2x2 3x3x3')
    study.add_argument('--block-study', nargs='+', type=int, metavar='FRAMES',
                       help='read the archived DOS from output and reblock its frame curves without rerunning QC')
    args = parser.parse_args(argv)
    try:
        source = args.control.resolve()
        raw = tomllib.loads(source.read_text())
        if set(raw) != {'dos'} or not isinstance(raw['dos'], dict):
            raise ValueError('DOS control requires exactly one [dos] table')
        if args.block_study is not None:
            from .block_study import run_block_study
            cfg = raw['dos']
            if 'kmesh-study' in cfg:
                raise ValueError('--block-study requires a single completed DOS output, not a kmesh study')
            directory = source.parent/Path(cfg.get('output', 'dos-results')).expanduser()
            mode = str(cfg['mode']).lower() if 'mode' in cfg else None
            result = run_block_study(directory, args.block_study, check_only=args.check,
                                     expected_mode=mode, resume=args.resume)
            print(f'wren: {result["status"]}; block study: {result["output"]}')
            return 0
        if args.kmesh_study:
            raw['dos']['kmesh-study'] = args.kmesh_study
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
