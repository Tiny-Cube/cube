"""Extract indexed monomer geometries and lists of periodic QM/MM neighbours."""
from contextlib import contextmanager
import csv
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import re
import tempfile
from zipfile import ZIP_DEFLATED, ZipFile, ZipInfo

import numpy as np


def file_hash(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as file:
        for block in iter(lambda: file.read(1024*1024), b''):
            digest.update(block)
    return digest.hexdigest()


def fingerprint(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, allow_nan=False).encode()).hexdigest()


def read_core(path):
    text = ' '.join(line.split('#', 1)[0] for line in Path(path).read_text().splitlines())
    tokens = re.split(r'[\s,]+', text.strip())
    if not tokens or any(not re.fullmatch('[1-9][0-9]*', t) for t in tokens):
        raise ValueError('fragment-file must contain positive one-based local atom numbers')
    indices = [int(t)-1 for t in tokens]
    if len(indices) != len(set(indices)):
        raise ValueError('fragment-file contains duplicate atom numbers')
    return indices


def contact_images(cores, lattice, cutoff, *, include_self=False):
    """All i<j and integer images with min(core_i, core_j+image@cell)<cutoff.

    If contact is possible, |center_j+image@cell-center_i|<Ri+Rj+cutoff.
    Multiplying this bound by each reciprocal-vector norm bounds every image
    integer. Unlike per-atom minimum images, one image moves the WHOLE molecule.
    """
    from scipy.spatial.distance import cdist
    lattice = np.asarray(lattice, dtype=float)
    if lattice.shape != (3, 3) or not np.isfinite(lattice).all() or np.linalg.det(lattice) <= 1e-8:
        raise ValueError('a finite, right-handed nondegenerate 3x3 periodic box is required')
    if not math.isfinite(cutoff) or cutoff <= 0:
        raise ValueError('cutoff must be positive and finite')
    cores = [np.asarray(c, dtype=float) for c in cores]
    if any(c.ndim != 2 or c.shape[1] != 3 or not len(c) or not np.isfinite(c).all() for c in cores):
        raise ValueError('each core must contain finite 3D coordinates')
    centers = [c.mean(axis=0) for c in cores]
    radii = [np.linalg.norm(c-center, axis=1).max() for c, center in zip(cores, centers)]
    inverse = np.linalg.inv(lattice)
    reciprocal_norm = np.linalg.norm(inverse, axis=0)
    pairs = itertools.combinations_with_replacement if include_self else itertools.combinations
    for i, j in pairs(range(len(cores)), 2):
        fraction = (centers[j]-centers[i])@inverse
        bounds = (radii[i]+radii[j]+cutoff)*reciprocal_norm
        lower = np.ceil(-fraction-bounds-1e-12).astype(int)
        upper = np.floor(-fraction+bounds+1e-12).astype(int)
        counts = upper-lower+1
        if np.any(counts <= 0):
            continue
        if math.prod(map(int, counts)) > 100000:
            raise ValueError('periodic image search exceeds 100000 candidates per pair; inspect box/core/cutoff')
        for shift in itertools.product(*(range(int(lo), int(hi)+1) for lo, hi in zip(lower, upper))):
            if i == j and shift <= (0, 0, 0):
                continue  # Reverse images are supplied by the neighbour list.
            shift = np.asarray(shift, dtype=int)
            distance = float(cdist(cores[i], cores[j]+shift@lattice).min())
            if distance < cutoff:
                yield i, j, shift, distance


def xyz_text(elements, coordinates, title):
    return f'{len(elements)}\n{title}; Angstrom\n'+''.join(
        f'{e:2s} {x:.10f} {y:.10f} {z:.10f}\n' for e, (x, y, z) in zip(elements, coordinates))


def selected_frames(start, step, count):
    """-1 selects only the last frame; otherwise sample from start to the end."""
    if type(start) is not int or start < -1 or type(step) is not int or step <= 0:
        raise ValueError('frame must be >= -1 and frame_step must be a positive integer')
    if count <= 0 or start >= count:
        raise ValueError(f'frame {start} exceeds {count} available frames (zero-based)')
    return [count-1] if start == -1 else list(range(start, count, step))


def _mm_shells(positions, lattice, cutoff):
    """Both directions of every whole-molecule contact, including self images."""
    shells = [set() for _ in positions]
    for i, j, shift, distance in contact_images(positions, lattice, cutoff, include_self=True):
        if distance < 0.5:
            raise ValueError('QM atom is within 0.5 Angstrom of an MM site; check QM exclusion/coordinates')
        shells[i].add((j, tuple(map(int, shift))))
        shells[j].add((i, tuple(map(int, -shift))))
    return shells


def _write_members(writer, qm, shell, mol_ids):
    # CSV contains only molecule numbers. Offsets refer to positions in that row,
    # so even two images of the same molecule can be reconstructed uniquely.
    members = sorted(shell)
    writer.writerow([qm, '-'.join(str(mol_ids[j]) for j, _ in members)])
    return {k: offset for k, (_, offset) in enumerate(members) if any(offset)}


def _extract_frame(frame, frame_count, read_frame, labels, atom_ids, settings, provenance, resume):
    base = Path(settings['output'])
    out = base/f'frame{frame:04d}'
    key = fingerprint({'provenance': provenance, 'frame': frame, 'extraction_version': 4})
    if out.exists():
        manifest_path = out/'manifest.json'
        if not resume or not manifest_path.is_file():
            raise ValueError(f'output exists: {out}; use --resume for identical inputs or choose another output')
        manifest = json.loads(manifest_path.read_text())
        if manifest.get('fingerprint') != key:
            raise ValueError('existing extraction belongs to different inputs/settings')
        for relative, expected in manifest['output_sha256'].items():
            if not (out/relative).is_file() or file_hash(out/relative) != expected:
                raise ValueError(f'cannot resume: missing/changed extracted file {relative}')
        return manifest_path
    molecules, lattice, time_ps, details = read_frame(frame)
    mol_ids = labels if labels is not None else list(range(1, len(molecules)+1))
    local = read_core(settings['fragment-file'])
    if not molecules or max(local) >= min(map(len, molecules)):
        raise ValueError('fragment-file atom number exceeds a complete molecule')
    if any(not np.array_equal(m.elements, molecules[0].elements) or not np.array_equal(m.names, molecules[0].names) for m in molecules):
        raise ValueError('one fragment-file currently requires identical molecules and consistent local atom ordering')
    if any(abs(float(m.charges.sum())) >= .01 for m in molecules):
        raise ValueError('this neutral closed-shell extraction requires neutral molecules')
    positions = [m.positions.astype(float).copy() for m in molecules]
    cores = [p[local] for p in positions]
    cutoff = settings.get('dimers_cutoff', settings.get('cutoff', 6.))
    contacts = list(contact_images(cores, lattice, cutoff))
    manifest = {'schema': 'wren-dimers-2', 'fingerprint': key, 'provenance': provenance,
                'frame': frame, 'frame_count': frame_count, 'time_ps': float(time_ps), **details,
                'coordinate_unit': 'Angstrom', 'distance_criterion': 'strict_minimum_selected_core_atom_distance',
                'cutoff_angstrom': cutoff, 'core_atom_numbers_one_based': [i+1 for i in local],
                'lattice_angstrom': lattice.tolist(), 'molecule_count': len(molecules),
                'monomer_archive': 'monomers.zip', 'cluster_table': 'clusters.csv', 'dimer_table': 'dimers.csv',
                'dimer_count': len(contacts), 'monomers': [], 'dimers': [], 'output_sha256': {}}
    mm_cutoff = settings.get('monomers-qmmm-cutoff', 12.)
    shells = _mm_shells(positions, lattice, mm_cutoff)
    base.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.wren-extract-', dir=base) as temporary:
        stage = Path(temporary)/'frame'
        stage.mkdir()
        with (ZipFile(stage/'monomers.zip', 'w', ZIP_DEFLATED) as archive,
              (stage/'clusters.csv').open('w', newline='') as table):
            writer = csv.writer(table)
            writer.writerow(['cluster_qm', 'cluster_mm'])
            for i, (identity, mol, xyz) in enumerate(zip(mol_ids, molecules, positions, strict=True)):
                member = f'frame{frame}.mol{identity}.xyz'
                # Store unwrapped MD positions once, at round-trip precision.
                lines = [str(len(xyz)), f'frame{frame}.mol{identity}; Angstrom']
                lines.extend(f'{e} {x:.17g} {y:.17g} {z:.17g}'
                             for e, (x, y, z) in zip(mol.elements, xyz, strict=True))
                archive.writestr(ZipInfo(member), '\n'.join(lines)+'\n', compress_type=ZIP_DEFLATED)
                offsets = _write_members(writer, str(identity), shells[i], mol_ids)
                record = {'id': str(identity), 'xyz_member': member,
                          'transport_position_angstrom': xyz.mean(axis=0).tolist(),
                          'global_atom_indices_zero_based': mol.indices.tolist(),
                          'cluster': {'cutoff_angstrom': mm_cutoff, 'mm_image_offsets': offsets}}
                if atom_ids is not None:
                    record['lammps_atom_ids'] = [atom_ids[j] for j in mol.indices]
                    record['lammps_molecule_id'] = identity
                manifest['monomers'].append(record)
        with (stage/'dimers.csv').open('w', newline='') as table:
            writer = csv.writer(table)
            writer.writerow(['dimer_qm', 'dimer_mm'])
            for i, j, image, distance in contacts:
                identity = f'{mol_ids[i]}_{mol_ids[j]}'
                shift = tuple(map(int, image))
                if any(shift):
                    identity += '_img_'+'_'.join(map(str, shift))
                shell = shells[i].union((k, tuple(a+b for a, b in zip(offset, shift)))
                                        for k, offset in shells[j])
                shell.difference_update([(i, (0, 0, 0)), (j, shift)])
                offsets = _write_members(writer, f'{mol_ids[i]}-{mol_ids[j]}', shell, mol_ids)
                ma, mb = molecules[i].masses[local], molecules[j].masses[local]
                if not np.all(ma > 0) or not np.all(mb > 0):
                    raise ValueError('positive topology core masses are required for core COM reporting')
                ca = np.average(cores[i], axis=0, weights=ma)
                cb = np.average(cores[j], axis=0, weights=mb) + image@lattice
                manifest['dimers'].append({'id': identity, 'mol_i': str(mol_ids[i]), 'mol_j': str(mol_ids[j]),
                    'image': shift, 'core_min_angstrom': distance, 'core_com_angstrom': float(np.linalg.norm(ca-cb)),
                    'atom_count_a': len(positions[i]), 'mm_image_offsets': offsets})
        for file in stage.rglob('*'):
            if file.is_file():
                manifest['output_sha256'][str(file.relative_to(stage))] = file_hash(file)
        (stage/'manifest.json').write_text(json.dumps(manifest, separators=(',', ':'), allow_nan=False)+'\n')
        os.replace(stage, out)
    return out/'manifest.json'


def _environment_frame(frame, lattice, molecules, provenance, labels=None):
    fragments = [{'indices': m.indices.tolist(), 'elements': m.elements.tolist(),
                  'coordinates': m.positions.astype(float).tolist(), 'charges': m.charges.astype(float).tolist()}
                 for m in molecules]
    if labels is not None:
        for record, label in zip(fragments, labels):
            record['molecule_id'] = str(label)
    source = 'LAMMPS_data_charges_dump_coordinates' if labels is not None else 'GROMACS_TPR_charges_XTC_coordinates'
    return {'lattice': lattice.tolist(), 'fragments': fragments,
            'provenance': {'source': source, 'frame': frame,
                           'sha256': {k:v for k,v in provenance['sha256'].items() if k != 'fragment-file'}}}


@contextmanager
def molecular_trajectory(settings):
    if 'lammps-data' in settings:
        from .lammps import LammpsTrajectory
        reader = LammpsTrajectory(settings)
        try:
            yield len(reader.offsets), reader.read, reader.mol_ids, reader.ids
        finally:
            reader.stream.close()
        return

    import MDAnalysis as mda
    from MDAnalysis.exceptions import NoDataError
    universe = mda.Universe(settings['tpr'], settings['xtc'])
    try:
        gro = mda.Universe(settings['gro'])
        try:
            if (len(gro.atoms) != len(universe.atoms) or not np.array_equal(gro.atoms.names, universe.atoms.names)
                    or not np.array_equal(gro.atoms.resindices, universe.atoms.resindices)):
                raise ValueError('GRO/TPR atom count, names or residue ordering do not match')
        finally:
            gro.trajectory.close()
        def read_frame(frame):
            ts = universe.trajectory[frame]
            universe.atoms.unwrap(compound='fragments', reference=None)
            return universe.atoms.fragments, np.asarray(ts.triclinic_dimensions, dtype=float), ts.time, {}
        yield len(universe.trajectory), read_frame, None, None
    except NoDataError as exc:
        raise ValueError(f'TPR must provide bonds, elements, masses and charges: {exc}') from exc
    finally:
        universe.trajectory.close()


def extract_dimers_monomers(settings, *, one_settings=None, two_settings=None, resume=False):
    """Store monomer XYZs and membership tables; legacy GJF settings are ignored."""
    if settings.get('format', ['xyz']) != ['xyz']:
        raise ValueError('extraction requires format = xyz; use xyz2gjf for Gaussian files')
    source_files = ('lammps-data', 'lammps-dump', 'element-map') if 'lammps-data' in settings else ('tpr','gro','xtc')
    provenance = {'settings': settings,
                  'sha256': {k: file_hash(settings[k]) for k in (*source_files, 'fragment-file')}}
    with molecular_trajectory(settings) as (count, read_frame, labels, atom_ids):
        frames = selected_frames(settings['frame'], settings.get('frame_step', 1), count)
        index_path = Path(settings['output'])/'frames.json'
        index_key = fingerprint({'provenance': provenance, 'frames': frames, 'extraction_version': 4})
        if len(frames) > 1 and index_path.exists():
            previous = json.loads(index_path.read_text())
            if not resume or previous.get('fingerprint') != index_key:
                raise ValueError('existing frame collection belongs to different inputs/settings; choose another output')
        entries = []
        for frame in frames:
            print(f'wren: extracting frame {frame}', flush=True)
            manifest_path = _extract_frame(frame, count, read_frame, labels, atom_ids, settings, provenance, resume)
            manifest = json.loads(manifest_path.read_text())
            entries.append({'frame': frame, 'time_ps': manifest['time_ps'],
                            'manifest': str(manifest_path.relative_to(index_path.parent)),
                            'sha256': file_hash(manifest_path), 'molecule_count': manifest['molecule_count'],
                            'dimer_count': manifest['dimer_count']})
        if len(entries) == 1:
            return manifest_path
        index = {'schema': 'cube-frames-1', 'coordinate_unit': 'Angstrom', 'fingerprint': index_key,
                 'provenance': provenance, 'frame_count': len(frames), 'trajectory_frame_count': count,
                 'frames': entries}
        with tempfile.NamedTemporaryFile(mode='w', prefix='.frames-', dir=index_path.parent, delete=False) as stream:
            temporary_index = Path(stream.name)
            json.dump(index, stream, indent=2, allow_nan=False)
            stream.write('\n')
        os.replace(temporary_index, index_path)
        return index_path


# Compatibility for callers of the original extraction helper.
extract_dimers = extract_dimers_monomers
