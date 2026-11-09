#!/usr/bin/env python3
"""Export a verified GROMACS QM/MM+LJ cluster for another compute machine."""
import argparse
from dataclasses import replace
import json
from pathlib import Path

from universal_kmc.quantum.geometry import read_xyz
from universal_kmc.workflow.dimers import file_hash
from universal_kmc.workflow.embedding import load_gromacs_environment, cluster_from_frame, read_charges, snapshot
from universal_kmc.workflow.lj import read_gromacs_lj, bind_gromacs_lj, geometry_key, sites_key


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--molecule', required=True)
    parser.add_argument('--cutoff', type=float, default=12.)
    parser.add_argument('--gmx-executable', default='gmx')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text())
    if manifest.get('schema') != 'cube-dimers-1' or manifest.get('coordinate_unit') != 'Angstrom':
        raise ValueError('expected cube-dimers-1 Angstrom manifest')
    records = [r for r in manifest['monomers'] if str(r['id']) == args.molecule]
    if len(records) != 1:
        raise ValueError('molecule must identify exactly one extracted full molecule')
    record = records[0]
    if args.output.exists():
        raise ValueError('choose a new output directory')
    frame = load_gromacs_environment(manifest)
    environment = cluster_from_frame(frame, record, args.cutoff)
    parameters = read_gromacs_lj(manifest['provenance']['settings']['tpr'], args.gmx_executable)
    model = bind_gromacs_lj(parameters, record, environment)
    args.output.mkdir(parents=True)
    xyz, pc, lj_file = (args.output/n for n in ('molecule.xyz', 'environment.pc', 'lj.json'))
    xyz.write_text(f'{len(record["elements"])}\nQM molecule {args.molecule}; Angstrom\n'+''.join(
        f'{e} {r[0]:.17g} {r[1]:.17g} {r[2]:.17g}\n' for e, r in zip(record['elements'], record['coordinates_angstrom'])))
    pc.write_text('# x y z q; Angstrom/e; frozen MM\n'+''.join(
        f'{r[0]:.17g} {r[1]:.17g} {r[2]:.17g} {q:.17g}\n' for r, q in zip(environment.coordinates_angstrom, environment.charges_e)))
    model = replace(model, provenance={**model.provenance,
        'qm_initial_geometry_fingerprint': geometry_key(read_xyz(xyz)), 'mm_sites_fingerprint': sites_key(read_charges(pc))})
    lj_file.write_text(json.dumps(model.payload(), indent=2)+'\n')
    snapshot(replace(environment, vdw=model), args.output/'source-mm-environment.json')
    metadata = {'qm_atoms': len(record['elements']), 'mm_atoms': len(environment.charges_e), 'frame': manifest['frame'],
                'molecule': args.molecule, 'cutoff_angstrom': args.cutoff,
                'file_sha256': {p.name: file_hash(p) for p in (xyz, pc, lj_file)}}
    (args.output/'export.json').write_text(json.dumps(metadata, indent=2)+'\n')
    print(json.dumps(metadata, indent=2))


if __name__ == '__main__':
    main()
