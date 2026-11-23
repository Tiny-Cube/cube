#!/usr/bin/env python3
"""Export one whole-molecule QM/MM cluster from a single-frame GROMACS manifest."""
import argparse
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import shutil
import tempfile

from wren.workflow.dimers import file_hash


def _write_json(path, data):
    path.write_text(json.dumps(data, indent=2, allow_nan=False) + '\n', encoding='utf-8')


def _read_topology(tpr):
    import MDAnalysis as mda
    # Read explicit TPR chemical bonds; never guess from interatomic distances.
    top = mda.Universe(str(tpr), to_guess=())
    neighbors = [[] for _ in top.atoms]
    for a, b in top.bonds.indices:
        neighbors[a].append(int(b))
        neighbors[b].append(int(a))
    return top.atoms.elements.tolist(), top.atoms.types.tolist(), top.atoms.charges.tolist(), neighbors


def _cluster_files(directory, manifest, record, environment, parameters, topology):
    elements, ff_names, charges, neighbors = topology
    lookup = {tuple(m['global_atom_indices_zero_based']): str(m['id']) for m in manifest['monomers']}
    central_ids = record['global_atom_indices_zero_based']
    groups = [(central_ids, [0, 0, 0], record['coordinates_angstrom'], 'QM')]
    for mapping in environment.provenance['molecular_images']:
        start, count = mapping['site_start'], mapping['site_count']
        ids = mapping['global_atom_indices_zero_based']
        if len(ids) != count or [charges[idx] for idx in ids] != list(environment.charges_e[start:start+count]):
            raise ValueError('frozen charge/order differs from original topology')
        groups.append((ids, mapping['image'], environment.coordinates_angstrom[start:start+count], 'MM'))
    metadata, bonds, lines = [], [], []
    for source_ids, image, coords, region in groups:
        if tuple(source_ids) not in lookup or len(source_ids) != len(coords):
            raise ValueError('incomplete molecule in the selected environment')
        mol_id = lookup[tuple(source_ids)]
        local = {source: len(metadata)+i for i, source in enumerate(source_ids, 1)}
        for source, xyz in zip(source_ids, coords):
            atom_id = local[source]
            ff = f"tpr-{parameters['types'][source]}" if parameters else ff_names[source]
            metadata.append({'atom_id': atom_id, 'molecule_id': mol_id, 'source_atom_id': source+1,
                             'image': list(image), 'element': elements[source], 'ff_type': ff, 'charge': charges[source]})
            row = ' '.join(format(float(v), '.17g') for v in xyz)
            lines.append(f"{elements[source]} {row} {region} {'T' if region == 'MM' else 'F'} {atom_id}")
            for neighbor in neighbors[source]:
                if neighbor not in local:
                    raise ValueError('selected molecule cuts a TPR bond')
                if neighbor > source:
                    bonds.append([atom_id, local[neighbor]])
    header = 'Properties=species:S:1:pos:R:3:region:S:1:frozen:L:1:atom_id:I:1 coordinate_unit=Angstrom'
    (directory/'coordinates.extxyz').write_text('\n'.join([str(len(lines)), header, *lines, '']), encoding='utf-8')
    _write_json(directory/'topology.json', {'schema': 'wren-qmmm-topology-1', 'charge_unit': 'e',
                                            'atoms': metadata, 'bonds': bonds})
    if parameters:
        lj = environment.vdw
        used = sorted({tuple(sorted((i, j))) for i in set(lj.qm_types) for j in set(lj.mm_types)})
        pairs = [{'types': [f'tpr-{i}', f'tpr-{j}'], 'c6': lj.c6_kj_mol_nm6[i][j],
                  'c12': lj.c12_kj_mol_nm12[i][j]} for i, j in used]
        _write_json(directory/'forcefield.json', {'schema': 'wren-qmmm-forcefield-1', 'model': 'LJ12-6',
                    'units': {'energy': 'kJ/mol', 'length': 'nm'},
                    'source': f"Resolved GROMACS TPR LJ_SR type pairs, SHA-256 {parameters['provenance']['tpr_sha256']}",
                    'pairs': pairs})


def _round_trip(directory, record, environment):
    import numpy as np
    from wren.workflow.separated_embedding import read_separated_embedding
    ff = directory/'forcefield.json' if environment.vdw else None
    qm, new, mol_id = read_separated_embedding(directory/'coordinates.extxyz', directory/'topology.json', ff)
    source_qm = list(zip(record['elements'], record['coordinates_angstrom']))
    if (mol_id != str(record['id']) or [e for e, _ in qm] != record['elements'] or
            not np.array_equal([r for _, r in qm], record['coordinates_angstrom']) or
            new.coordinates_angstrom != environment.coordinates_angstrom or new.charges_e != environment.charges_e):
        raise ValueError('QM atoms or frozen MM geometry/charges changed on export')
    check = {'qm_elements_and_order_equal': True, 'qm_coordinates_exact': True,
             'mm_coordinate_order_exact': True, 'mm_charge_order_exact': True}
    if environment.vdw:
        old_e, old_g = environment.classical(source_qm)
        new_e, new_g = new.classical(qm)
        grad_error = float(np.max(np.abs(new_g-old_g)))
        if abs(new_e-old_e) > 1e-12 or grad_error > 1e-12:
            raise ValueError('export changed QM-MM LJ energy or gradient')
        check.update(original_lj_hartree=old_e, imported_lj_hartree=new_e,
                     max_qm_lj_gradient_difference_hartree_per_bohr=grad_error,
                     max_abs_reference_qm_lj_gradient_hartree_per_bohr=float(np.max(np.abs(old_g))))
    return check


def export_cluster(manifest_file, molecule, directory, *, cutoff=None, lj=False, gmx='gmx'):
    from wren.workflow.embedding import cluster_from_frame, load_gromacs_environment
    from wren.workflow.lj import bind_gromacs_lj, read_gromacs_lj

    path = Path(manifest_file).resolve()
    directory = Path(directory).resolve()
    if directory.exists():
        raise FileExistsError(f'output must be a new directory: {directory}')
    raw = path.read_bytes()
    data = json.loads(raw)
    if data.get('schema') != 'cube-dimers-1' or data.get('coordinate_unit') != 'Angstrom':
        raise ValueError('select a single frame manifest.json with cube-dimers-1 and Angstrom coordinates')
    selected = [m for m in data['monomers'] if str(m['id']) == str(molecule)]
    if len(selected) != 1:
        raise ValueError('molecule must identify exactly one monomer in the frame')
    record = selected[0]
    if cutoff is None:
        cutoff = record.get('cluster', {}).get('cutoff_angstrom')
    if cutoff is None:
        raise ValueError('no cluster cutoff recorded for this molecule; supply --cutoff in Angstrom')
    frame = load_gromacs_environment(data)
    env = cluster_from_frame(frame, record, cutoff)
    tpr = Path(data['provenance']['settings']['tpr'])
    topology = _read_topology(tpr)
    parameters = read_gromacs_lj(tpr, gmx) if lj else None
    if file_hash(tpr) != frame['provenance']['sha256']['tpr']:
        raise ValueError('TPR changed since reading the MD frame')
    if lj:
        env = replace(env, vdw=bind_gromacs_lj(parameters, record, env))
    directory.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.wren-export-', dir=directory.parent) as temp:
        stage = Path(temp)
        _cluster_files(stage, data, record, env, parameters, topology)
        verification = _round_trip(stage, record, env)
        report = {'schema': 'wren-qmmm-cluster-export-1', 'frame': data['frame'], 'molecule': str(molecule),
                  'cutoff_angstrom': cutoff, 'qm_atom_count': len(record['elements']), 'mm_atom_count': len(env.charges_e),
                  'mm_molecular_image_count': len(env.provenance['molecular_images']), 'mm_charge_e': sum(env.charges_e),
                  'bare_qm_mm_lj': lj, 'manifest_sha256': hashlib.sha256(raw).hexdigest(),
                  'source_sha256': frame['provenance']['sha256'],
                  'file_sha256': {p.name: file_hash(p) for p in sorted(stage.iterdir())}, 'verification': verification}
        _write_json(stage/'export.json', report)
        # Claim a new directory after scientific validation. Exclusive writes
        # also work on removable filesystems without hard-link support.
        directory.mkdir()
        published = []
        try:
            for p in stage.iterdir():
                target = directory/p.name
                with target.open('xb') as dest, p.open('rb') as source:
                    published.append(target)
                    shutil.copyfileobj(source, dest)
        except OSError:
            for p in published:
                p.unlink()
            if not any(directory.iterdir()):
                directory.rmdir()
            raise
    return report


def _legacy_export(args):
    """Preserve the original XYZ, charge-list and bound LJ export."""
    from wren.quantum.geometry import read_xyz
    from wren.workflow.embedding import load_gromacs_environment, cluster_from_frame, read_charges, snapshot
    from wren.workflow.lj import read_gromacs_lj, bind_gromacs_lj, geometry_key, sites_key
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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('frame_manifest', nargs='?', help='single-frame manifest.json with accessible original TPR/GRO/XTC')
    parser.add_argument('--manifest', type=Path, help='legacy invocation; defaults to XYZ + charge list + bound LJ')
    parser.add_argument('--format', choices=('legacy', 'separated'), help='override default export format')
    parser.add_argument('--molecule', required=True, help='manifest molecule ID')
    parser.add_argument('--cutoff', type=float, help='strict shortest atom distance, Angstrom; default from stored cluster')
    parser.add_argument('--lj', action='store_true', help='export resolved QM-MM C6/C12 pairs via gmx dump')
    parser.add_argument('--gmx', '--gmx-executable', dest='gmx_executable', default='gmx', help='GROMACS executable for --lj')
    parser.add_argument('--output', type=Path, required=True, help='new directory')
    args = parser.parse_args()
    if bool(args.frame_manifest) == bool(args.manifest):
        parser.error('provide exactly one positional manifest or --manifest')
    kind = args.format or ('legacy' if args.manifest else 'separated')
    args.manifest = args.manifest or Path(args.frame_manifest)
    try:
        if kind == 'legacy':
            if args.cutoff is None:
                args.cutoff = 12.
            _legacy_export(args)
            return
        result = export_cluster(args.manifest, args.molecule, args.output, cutoff=args.cutoff, lj=args.lj,
                                gmx=args.gmx_executable)
    except (ValueError, OSError, KeyError, RuntimeError) as exc:
        parser.exit(1, f'wren: {exc}\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
