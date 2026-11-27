"""Read a monomer by frame/molecule key or assemble a periodic QM pair in memory."""
import argparse
import json
from pathlib import Path
import re
import sys
from zipfile import BadZipFile, ZipFile

from wren.quantum.geometry import parse_xyz


def indexed_xyz(root, key):
    match = re.fullmatch(r'frame([0-9]+)\.mol([0-9]+)', key)
    if not match:
        raise ValueError('expected frameX.molY, e.g. frame2000.mol8')
    frame, molecule = map(int, match.groups())
    with ZipFile(Path(root)/f'frame{frame:04d}'/'monomers.zip') as archive:
        return archive.read(f'frame{frame}.mol{molecule}.xyz').decode('utf-8')


def load_frame(path, data=None):
    """Load centered molecular coordinates for existing QM, embedding and DOS clients."""
    path = Path(path)
    if data is None:
        data = json.loads(path.read_text())
    if data.get('coordinate_unit') != 'Angstrom' or data.get('schema') not in ('cube-dimers-1', 'wren-dimers-2'):
        raise ValueError('expected a single Angstrom extraction manifest')
    if data['schema'] == 'cube-dimers-1':
        return data
    import numpy as np
    from .dimers import file_hash
    archive_path = path.parent/'monomers.zip'
    if file_hash(archive_path) != data['output_sha256']['monomers.zip']:
        raise ValueError(f'monomer archive changed since extraction: {archive_path}')
    with ZipFile(archive_path) as archive:
        for record in data['monomers']:
            atoms = parse_xyz(archive.read(record['xyz_member']).decode('utf-8'))
            xyz = np.array([r for _, r in atoms])
            record['elements'] = [e for e, _ in atoms]
            record['coordinates_angstrom'] = (xyz-record['transport_position_angstrom']).tolist()
    return data


def pair_geometry(pair, monomers, lattice):
    """Return the QM atoms centered together, in monomer A then image B order."""
    if 'coordinates_angstrom' in pair:  # Inline coordinates in legacy manifests.
        return pair
    import numpy as np
    a, b = monomers[pair['mol_i']], monomers[pair['mol_j']]
    xyz_a = np.asarray(a['coordinates_angstrom']) + a['transport_position_angstrom']
    xyz_b = np.asarray(b['coordinates_angstrom']) + b['transport_position_angstrom']
    xyz = np.concatenate((xyz_a, xyz_b + np.asarray(pair['image']) @ lattice))
    return {**pair, 'elements': a['elements'] + b['elements'],
            'coordinates_angstrom': (xyz-xyz.mean(axis=0)).tolist()}


def main(argv=None):
    parser = argparse.ArgumentParser(description='Read one complete XYZ from a Wren extraction.')
    parser.add_argument('key', help='frameX.molY')
    parser.add_argument('--root', type=Path, default=Path('mol'))
    parser.add_argument('-o', '--output', type=Path, help='XYZ destination (default: stdout)')
    args = parser.parse_args(argv)
    try:
        text = indexed_xyz(args.root, args.key)
        if args.output is None:
            sys.stdout.write(text)
        else:
            args.output.write_text(text, encoding='utf-8')
    except (ValueError, KeyError, OSError, BadZipFile) as exc:
        parser.exit(1, f'wren-xyz: {exc}\n')


if __name__ == '__main__':
    main()
