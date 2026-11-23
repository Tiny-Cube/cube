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
    parser.add_argument('--kmesh-study', nargs='+', metavar='NxNxN',
                        help='compare successively refined periodic meshes, e.g. 1x1x1 2x2x2 3x3x3')
    args = parser.parse_args(argv)
    try:
        source = args.control.resolve()
        raw = tomllib.loads(source.read_text())
        if set(raw) != {'dos'}:
            raise ValueError('DOS control requires exactly one [dos] table')
        settings = normalize(raw['dos'], source.parent)
        meshes = []
        for token in args.kmesh_study or ():
            mesh = normalize(dict(raw['dos'], kmesh=token.lower().split('x')), source.parent)['kmesh']
            if meshes and (mesh == meshes[-1] or any(a < b for a, b in zip(mesh, meshes[-1]))):
                raise ValueError('each study mesh must refine at least one axis and never coarsen another')
            meshes.append(mesh)
        if meshes and (settings['mode'] != 'periodic' or len(meshes) < 2):
            raise ValueError('--kmesh-study needs periodic mode and at least two distinct meshes')
        if args.check:
            print(json.dumps(dict(settings, study_meshes=meshes) if meshes else settings, indent=2))
            return 0
        if meshes:
            from .mesh_study import run_study
            report = run_study(settings, meshes, resume=args.resume)
            print(f'wren: completed mesh study: {report}')
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
