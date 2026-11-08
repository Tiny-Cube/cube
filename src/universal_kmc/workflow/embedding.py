"""Frozen point-charge electrostatic embedding, without MM relaxation or Ewald.

PySCF model and energy terms: https://pyscf.org/user/qmmm.html
The full target molecule is QM; covalent QM/MM cuts/link atoms are unsupported.
"""
from dataclasses import dataclass
import itertools
import json
import math
from pathlib import Path

import numpy as np

from .dimers import file_hash, fingerprint


@dataclass(frozen=True)
class FrozenMM:
    coordinates_angstrom: tuple
    charges_e: tuple
    provenance: dict
    vdw: object = None

    def __post_init__(self):
        coordinates = np.asarray(self.coordinates_angstrom, dtype=float)
        charges = np.asarray(self.charges_e, dtype=float)
        if (coordinates.ndim != 2 or coordinates.shape[1:] != (3,) or not len(coordinates)
                or charges.shape != (len(coordinates),) or not np.isfinite(coordinates).all()
                or not np.isfinite(charges).all()):
            raise ValueError('MM environment needs nonempty finite Nx3 coordinates and N charges')
        # Own immutable coordinates/charges even when the caller supplies arrays.
        object.__setattr__(self, 'coordinates_angstrom', tuple(tuple(map(float, r)) for r in coordinates))
        object.__setattr__(self, 'charges_e', tuple(map(float, charges)))
        if self.vdw is not None and len(self.vdw.mm_types) != len(charges):
            raise ValueError('LJ MM type count differs from frozen charge/site count')

    def payload(self):
        result = {'schema': 'cube-frozen-mm-1', 'model': 'fixed_point_charge_electrostatic_embedding',
                'coordinate_unit': 'Angstrom', 'charge_unit': 'elementary_charge',
                'coordinates_angstrom': self.coordinates_angstrom, 'charges_e': self.charges_e,
                'atom_count': len(self.charges_e), 'total_charge_e': float(sum(self.charges_e)),
                'mm_frozen': True, 'charges_frozen': True, 'electrostatics': 'finite_cluster_no_Ewald',
                'qm_mm_vdw_included': self.vdw is not None, 'polarizable_mm': False, 'provenance': self.provenance}
        if self.vdw is not None:
            result.update(model='fixed_point_charge_electrostatic_embedding_plus_LJ12_6',
                          lennard_jones=self.vdw.payload())
        return result

    def classical(self, atoms):
        if self.vdw is None:
            return 0., np.zeros((len(atoms), 3))
        return self.vdw.evaluate([r for e, r in atoms], self.coordinates_angstrom)

    @property
    def fingerprint(self):
        return fingerprint(self.payload())

    def validate_qm(self, atoms):
        from scipy.spatial.distance import cdist
        qm = np.asarray([a[1] for a in atoms], dtype=float)
        if not len(qm) or not np.isfinite(qm).all():
            raise ValueError('QM coordinates must be nonempty and finite')
        if float(cdist(qm, self.coordinates_angstrom).min()) < 0.5:
            raise ValueError('QM atom is within 0.5 Angstrom of an MM site; check QM exclusion/coordinates')

    def assert_attached(self, mf):
        mm = getattr(mf, 'mm_mol', None)
        if (mm is None or mm.atom_coords().shape != (len(self.charges_e), 3)
                or mm.atom_charges().shape != (len(self.charges_e),)
                or not np.allclose(mm.atom_coords(unit='Angstrom'), self.coordinates_angstrom,
                                         rtol=0, atol=1e-10)
                or not np.allclose(mm.atom_charges(), self.charges_e, rtol=0, atol=1e-14)):
            raise RuntimeError('frozen MM coordinates/charges changed or embedding was lost')


def read_charges(filename):
    """Four whitespace columns x y z q; Angstrom/e; optional # comments."""
    path = Path(filename)
    raw = path.read_bytes()
    rows = []
    for line_number, line in enumerate(raw.decode('utf-8-sig').splitlines(), 1):
        line = line.split('#', 1)[0].strip()
        if not line:
            continue
        try:
            values = list(map(float, line.split()))
            if len(values) != 4:
                raise ValueError('expected x y z q')
            rows.append(values)
        except ValueError as exc:
            raise ValueError(f'{path}:{line_number}: expected four numeric columns x y z q (Angstrom/e)') from exc
    import hashlib
    return FrozenMM(tuple(tuple(r[:3]) for r in rows), tuple(r[3] for r in rows),
                    {'source': 'explicit_file', 'path': str(path.resolve()),
                     'sha256': hashlib.sha256(raw).hexdigest(), 'selection': 'all supplied sites'})


def load_gromacs_environment(manifest):
    """Recover the original frame with verified topology charges and atom IDs.

    Reading the original frame also supports existing cube-dimers-1 manifests,
    which intentionally only stored individually centered QM geometries.
    """
    import MDAnalysis as mda
    provenance = manifest.get('provenance', {})
    paths = provenance.get('settings', {})
    expected = provenance.get('sha256', {})
    for key in ('tpr', 'gro', 'xtc'):
        if key not in paths or key not in expected:
            raise ValueError('MM source needs a cube GROMACS manifest with original TPR/GRO/XTC hashes')
        if file_hash(paths[key]) != expected[key]:
            raise ValueError(f'original {key} changed since extraction; cannot construct MM environment')
    u = mda.Universe(paths['tpr'], paths['xtc'])
    gro = mda.Universe(paths['gro'])
    try:
        if (len(u.atoms) != len(gro.atoms) or not np.array_equal(u.atoms.names, gro.atoms.names)
                or not np.array_equal(u.atoms.resindices, gro.atoms.resindices)):
            raise ValueError('GRO/TPR atom ordering differs')
        frame = manifest['frame']
        if type(frame) is not int or frame < 0 or frame >= len(u.trajectory):
            raise ValueError('invalid original frame index in manifest')
        u.trajectory[frame]
        u.atoms.unwrap(compound='fragments', reference=None)
        lattice = np.asarray(u.trajectory.ts.triclinic_dimensions, dtype=float)
        if not np.allclose(lattice, manifest['lattice_angstrom'], rtol=0, atol=1e-6):
            raise ValueError('manifest box differs from original frame')
        fragments = [{'indices': m.indices.tolist(), 'elements': m.elements.tolist(),
                      'coordinates': m.positions.astype(float).tolist(), 'charges': m.charges.astype(float).tolist()}
                     for m in u.atoms.fragments]
        return {'lattice': lattice.tolist(), 'fragments': fragments,
                'provenance': {'source': 'GROMACS_TPR_charges_XTC_coordinates', 'frame': frame,
                               'sha256': {k: expected[k] for k in ('tpr', 'gro', 'xtc')}}}
    finally:
        u.trajectory.close()
        gro.trajectory.close()


def select_mm_sites(frame, record, cutoff):
    """Select whole molecular images; extraction also permits an empty shell."""
    from scipy.spatial.distance import cdist
    if not math.isfinite(cutoff) or cutoff <= 0:
        raise ValueError('MM cutoff must be positive and finite')
    lattice = np.asarray(frame['lattice'], dtype=float)
    if lattice.shape != (3, 3) or not np.isfinite(lattice).all() or np.linalg.det(lattice) <= 1e-8:
        raise ValueError('MM cluster requires a right-handed nondegenerate periodic box')
    ids = record.get('global_atom_indices_zero_based')
    if not ids or len(set(ids)) != len(ids):
        raise ValueError('QM record needs unique global atom indices from extraction')
    candidates = [m for m in frame['fragments'] if m['indices'] == ids]
    if len(candidates) != 1:
        raise ValueError('QM region must be one complete topology fragment; covalent QM/MM cuts are unsupported')
    target = candidates[0]
    qm = np.asarray(record['coordinates_angstrom'], dtype=float)
    origin_qm = np.asarray(target['coordinates'], dtype=float)
    if record['elements'] != target['elements'] or qm.shape != origin_qm.shape:
        raise ValueError('QM atom elements/order differ from topology')
    offset = qm-origin_qm
    if not np.allclose(offset, offset[0], rtol=0, atol=1e-6):
        raise ValueError('QM geometry must match original MD frame up to a common translation')
    center = origin_qm.mean(axis=0)
    radius = np.linalg.norm(origin_qm-center, axis=1).max()
    inverse = np.linalg.inv(lattice)
    reciprocal = np.linalg.norm(inverse, axis=0)
    coordinates, charges, elements, mapping = [], [], [], []
    for molecule_index, fragment in enumerate(frame['fragments'], 1):
        xyz = np.asarray(fragment['coordinates'], dtype=float)
        q = np.asarray(fragment['charges'], dtype=float)
        if xyz.shape != (len(fragment['indices']), 3) or q.shape != (len(xyz),) or not np.isfinite(xyz).all() or not np.isfinite(q).all():
            raise ValueError('invalid MM topology coordinates/charges')
        c = xyz.mean(axis=0)
        bound = (radius+np.linalg.norm(xyz-c, axis=1).max()+cutoff)*reciprocal
        fraction = (c-center)@inverse
        lower, upper = np.ceil(-fraction-bound-1e-12).astype(int), np.floor(-fraction+bound+1e-12).astype(int)
        if math.prod(max(0, int(hi-lo+1)) for lo, hi in zip(lower, upper)) > 100000:
            raise ValueError('MM periodic image search exceeds 100000 candidates per molecule')
        for image in itertools.product(*(range(int(lo), int(hi)+1) for lo, hi in zip(lower, upper))):
            if fragment is target and image == (0, 0, 0):
                continue
            moved = xyz+np.asarray(image)@lattice
            if float(cdist(origin_qm, moved).min()) >= cutoff:
                continue
            coordinates.extend((moved+offset[0]).tolist())
            charges.extend(q.tolist())
            elements.extend(fragment['elements'])
            mapping.append({'global_atom_indices_zero_based': fragment['indices'], 'image': image,
                            'site_start': len(charges)-len(q), 'site_count': len(q),
                            'molecule_id': str(molecule_index)})
    provenance = {**frame['provenance'], 'qm_global_atom_indices_zero_based': ids,
         'cutoff_angstrom': cutoff, 'selection': 'whole molecular images with strict min QM-MM distance < cutoff',
         'qm_frame_translation_angstrom': offset[0].tolist(), 'molecular_images': mapping,
         'selection_fixed_at_initial_geometry': True}
    if coordinates and float(cdist(qm, coordinates).min()) < 0.5:
        raise ValueError('QM atom is within 0.5 Angstrom of an MM site; check QM exclusion/coordinates')
    return {'coordinates_angstrom': coordinates, 'charges_e': charges, 'elements': elements, 'provenance': provenance}


def cluster_from_frame(frame, record, cutoff):
    """Construct a nonempty fixed MM environment for quantum calculations."""
    sites = select_mm_sites(frame, record, cutoff)
    # Preserve existing calculation fingerprints; the molecule label is export metadata.
    provenance = {**sites['provenance'], 'molecular_images': [
        {k:v for k,v in image.items() if k != 'molecule_id'}
        for image in sites['provenance']['molecular_images']]}
    return FrozenMM(sites['coordinates_angstrom'], sites['charges_e'], provenance)


def snapshot(environment, filename):
    payload = environment.payload()
    payload['fingerprint'] = environment.fingerprint
    Path(filename).write_text(json.dumps(payload, indent=2, allow_nan=False)+'\n')
