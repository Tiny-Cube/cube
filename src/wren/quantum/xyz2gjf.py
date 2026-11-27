"""Convert one ordinary XYZ or indexed monomer to a Gaussian single-point input."""
import argparse
from pathlib import Path
from zipfile import BadZipFile

from .geometry import parse_xyz, read_xyz
from wren.workflow.frame_geometry import indexed_xyz


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', help='single-frame .xyz path, or frameX.molY with --root')
    parser.add_argument('output', type=Path, help='destination .gjf')
    parser.add_argument('--root', type=Path, help='mol directory for indexed input')
    parser.add_argument('--method', default='B3LYP')
    parser.add_argument('--basis', default='6-31G(d,p)')
    parser.add_argument('--charge', type=int, default=0)
    parser.add_argument('--multiplicity', type=int, default=1)
    args = parser.parse_args(argv)
    if args.multiplicity < 1:
        parser.error('--multiplicity must be positive')
    try:
        if args.root is None:
            if args.output.resolve() == Path(args.source).resolve():
                raise ValueError('output must differ from input XYZ')
            atoms = read_xyz(args.source)
        else:
            atoms = parse_xyz(indexed_xyz(args.root, args.source))
        lines = [f'#p {args.method}/{args.basis} SP NoSymm Units=Angstrom', '',
                 args.output.stem, '', f'{args.charge} {args.multiplicity}']
        lines.extend(f'{e} {x:.17g} {y:.17g} {z:.17g}' for e, (x, y, z) in atoms)
        args.output.write_text('\n'.join(lines)+'\n\n', encoding='utf-8')
    except (KeyError, ValueError, OSError, BadZipFile) as exc:
        parser.exit(1, f'xyz2gjf: {exc}\n')


if __name__ == '__main__':
    main()
