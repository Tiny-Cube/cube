"""Full-molecule dimer extraction with explicit whole-molecule PBC images.

The user's extract_dimers.py defines strict minimum selected-core atom
distance. This implementation retains that criterion, independently writes
the extraction/manifest layer, and needs no Gaussian/Slurm/MOMAP templates.
"""
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import re
import tempfile

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


def contact_images(cores, lattice, cutoff):
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
    for i, j in itertools.combinations(range(len(cores)), 2):
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
            shift = np.asarray(shift, dtype=int)
            distance = float(cdist(cores[i], cores[j]+shift@lattice).min())
            if distance < cutoff:
                yield i, j, shift, distance


def xyz_text(elements, coordinates, title):
    return f'{len(elements)}\n{title}; Angstrom\n'+''.join(
        f'{e:2s} {x:.10f} {y:.10f} {z:.10f}\n' for e, (x, y, z) in zip(elements, coordinates))


def gjf_text(elements, coordinates, title, settings):
    method, basis = settings.get('method', 'b3lyp'), settings.get('basis-set', '6-31G(d,p)')
    threads, memory = settings.get('threads', 4), settings.get('max-memory-mb', 4000)
    return (f'%chk={title}.chk\n%nprocshared={threads}\n%mem={memory}MB\n'
            f'#p {method}/{basis} SP NoSymm SCF=(XQC,MaxCycle=512)\n\n'
            f'{title}; full frozen MD molecules; Angstrom\n\n0 1\n'+''.join(
                f'{e:2s} {x:.10f} {y:.10f} {z:.10f}\n' for e, (x, y, z) in zip(elements, coordinates))+'\n')


def _write_formats(directory, stem, elements, coordinates, formats, settings):
    files = {}
    for extension in formats:
        relative = f'{stem}.{extension}'
        renderer = xyz_text if extension == 'xyz' else gjf_text
        args = (elements, coordinates, Path(stem).name)
        text = renderer(*args) if extension == 'xyz' else renderer(*args, settings)
        target = directory/relative
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(text)
        files[extension] = relative
    return files


def extract_dimers(settings, *, one_settings=None, two_settings=None, resume=False):
    """TPR topology, GRO order check, XTC coordinates; output only chosen formats."""
    import MDAnalysis as mda
    from MDAnalysis.exceptions import NoDataError
    header_keys = ('method', 'basis-set', 'threads', 'max-memory-mb')
    one_settings = {k: v for k, v in (one_settings or {}).items() if k in header_keys}
    two_settings = {k: v for k, v in (two_settings or {}).items() if k in header_keys}
    provenance = {'settings': settings, 'one_gjf': one_settings, 'two_gjf': two_settings,
                  'sha256': {k: file_hash(settings[k]) for k in ('tpr', 'gro', 'xtc', 'fragment-file')}}
    key = fingerprint(provenance)
    universe = mda.Universe(settings['tpr'], settings['xtc'])
    gro = mda.Universe(settings['gro'])
    try:
        if len(gro.atoms) != len(universe.atoms) or not np.array_equal(gro.atoms.names, universe.atoms.names) or not np.array_equal(gro.atoms.resindices, universe.atoms.resindices):
            raise ValueError('GRO/TPR atom count, names or residue ordering do not match')
        if settings['frame'] >= len(universe.trajectory):
            raise ValueError(f'frame {settings["frame"]} exceeds {len(universe.trajectory)} available frames (zero-based)')
        universe.trajectory[settings['frame']]
        frame = universe.trajectory.frame
        base = Path(settings['output'])
        out = base/f'frame{frame:04d}'
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
        universe.atoms.unwrap(compound='fragments', reference=None)
        molecules = universe.atoms.fragments
        local = read_core(settings['fragment-file'])
        if not molecules or max(local) >= min(map(len, molecules)):
            raise ValueError('fragment-file atom number exceeds a complete molecule')
        if any(not np.array_equal(m.elements, molecules[0].elements) or not np.array_equal(m.names, molecules[0].names) for m in molecules):
            raise ValueError('one fragment-file currently requires identical molecules and consistent local atom ordering')
        if any(abs(float(m.charges.sum())) >= .01 for m in molecules):
            raise ValueError('this neutral closed-shell extraction requires neutral molecules')
        lattice = np.asarray(universe.trajectory.ts.triclinic_dimensions, dtype=float)
        positions = [m.positions.astype(float).copy() for m in molecules]
        cores = [p[local] for p in positions]
        contacts = list(contact_images(cores, lattice, settings['cutoff']))
        manifest = {'schema': 'cube-dimers-1', 'fingerprint': key, 'provenance': provenance,
                    'frame': frame, 'frame_count': len(universe.trajectory), 'time_ps': float(universe.trajectory.time),
                    'coordinate_unit': 'Angstrom', 'distance_criterion': 'strict_minimum_selected_core_atom_distance',
                    'cutoff_angstrom': settings['cutoff'], 'core_atom_numbers_one_based': [i+1 for i in local],
                    'lattice_angstrom': lattice.tolist(), 'molecule_count': len(molecules),
                    'dimer_count': len(contacts), 'monomers': [], 'dimers': [], 'output_sha256': {}}
        base.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix='.cube-extract-', dir=base) as temporary:
            stage = Path(temporary)/'frame'
            stage.mkdir()
            for i, (mol, xyz) in enumerate(zip(molecules, positions), 1):
                xyz = xyz-xyz.mean(axis=0)
                elements = mol.elements.tolist()
                files = _write_formats(stage, f'monomers/{i}', elements, xyz, settings['format'], one_settings)
                manifest['monomers'].append({'id': str(i), 'elements': elements, 'coordinates_angstrom': xyz.tolist(),
                    'global_atom_indices_zero_based': mol.indices.tolist(), 'files': files})
            for i, j, image, distance in contacts:
                identity = f'{i+1}_{j+1}'
                if np.any(image):
                    identity += '_img_'+'_'.join(str(int(n)) for n in image)
                a, b = positions[i].copy(), positions[j]+image@lattice
                center = np.concatenate((a, b)).mean(axis=0)
                a, b = a-center, b-center
                elements = molecules[i].elements.tolist()+molecules[j].elements.tolist()
                xyz = np.concatenate((a, b))
                files = _write_formats(stage, f'dimers/{identity}', elements, xyz, settings['format'], two_settings)
                masses_a, masses_b = molecules[i].masses[local], molecules[j].masses[local]
                if not np.all(masses_a > 0) or not np.all(masses_b > 0):
                    raise ValueError('positive TPR core masses are required for core COM reporting')
                com = np.linalg.norm(np.average(a[local], axis=0, weights=masses_a)-np.average(b[local], axis=0, weights=masses_b))
                manifest['dimers'].append({'id': identity, 'mol_i': str(i+1), 'mol_j': str(j+1),
                    'image': image.tolist(), 'core_min_angstrom': distance, 'core_com_angstrom': float(com),
                    'elements': elements, 'coordinates_angstrom': xyz.tolist(), 'atom_count_a': len(a), 'files': files})
            import csv
            with (stage/'dimers.csv').open('w', newline='') as file:
                keys = ('id', 'mol_i', 'mol_j', 'image_a', 'image_b', 'image_c', 'core_min_angstrom', 'core_com_angstrom')
                writer = csv.DictWriter(file, fieldnames=keys)
                writer.writeheader()
                for pair in manifest['dimers']:
                    writer.writerow({**{k: pair[k] for k in ('id','mol_i','mol_j','core_min_angstrom','core_com_angstrom')},
                                     **dict(zip(('image_a','image_b','image_c'), pair['image']))})
            for file in stage.rglob('*'):
                if file.is_file():
                    manifest['output_sha256'][str(file.relative_to(stage))] = file_hash(file)
            (stage/'manifest.json').write_text(json.dumps(manifest, indent=2, allow_nan=False)+'\n')
            os.replace(stage, out)
        return out/'manifest.json'
    except NoDataError as exc:
        raise ValueError(f'TPR must provide bonds, elements, masses and charges: {exc}') from exc
    finally:
        universe.trajectory.close()
        gro.trajectory.close()
