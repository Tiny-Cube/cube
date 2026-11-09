"""Bare finite-cluster QM-MM LJ12-6, using resolved TPR pair coefficients.

V(r)=C12/r**12-C6/r**6. No intramolecular QM force-field terms or MM motion.
https://manual.gromacs.org/current/reference-manual/functions/nonbonded-interactions.html
"""
from dataclasses import dataclass
import json
from pathlib import Path
import re
import subprocess

import numpy as np

from .dimers import file_hash, fingerprint


def geometry_key(atoms):
    return fingerprint([(str(e), tuple(map(float, r))) for e, r in atoms])


def sites_key(environment):
    return fingerprint({'coordinates_angstrom': environment.coordinates_angstrom, 'charges_e': environment.charges_e})


@dataclass(frozen=True)
class LennardJones:
    qm_types: tuple
    mm_types: tuple
    c6_kj_mol_nm6: tuple
    c12_kj_mol_nm12: tuple
    provenance: dict

    def __post_init__(self):
        c6, c12 = np.asarray(self.c6_kj_mol_nm6, dtype=float), np.asarray(self.c12_kj_mol_nm12, dtype=float)
        if (c6.ndim != 2 or c6.shape[0] != c6.shape[1] or not c6.size or c12.shape != c6.shape
                or not np.isfinite(c6).all() or not np.isfinite(c12).all() or np.any(c6 < 0) or np.any(c12 < 0)
                or not np.allclose(c6, c6.T, rtol=1e-10, atol=0) or not np.allclose(c12, c12.T, rtol=1e-10, atol=0)):
            raise ValueError('LJ matrices must be finite nonnegative symmetric square C6/C12 tables')
        for name in ('qm_types', 'mm_types'):
            values = np.asarray(getattr(self, name))
            if (values.ndim != 1 or not len(values) or values.dtype.kind not in 'iu'
                    or np.any(values < 0) or np.any(values >= len(c6))):
                raise ValueError('LJ atom type arrays must contain valid zero-based integer type indices')
            object.__setattr__(self, name, tuple(map(int, values)))
        object.__setattr__(self, 'c6_kj_mol_nm6', tuple(tuple(map(float, r)) for r in c6))
        object.__setattr__(self, 'c12_kj_mol_nm12', tuple(tuple(map(float, r)) for r in c12))

    def payload(self):
        return {'schema': 'cube-qmmm-lj-1', 'model': 'bare_LJ12_6_finite_cluster_no_switch_no_tail_correction',
                'qm_types': self.qm_types, 'mm_types': self.mm_types,
                'c6_kj_mol_nm6': self.c6_kj_mol_nm6, 'c12_kj_mol_nm12': self.c12_kj_mol_nm12,
                'provenance': self.provenance, 'implementation_sha256': file_hash(__file__)}

    def evaluate(self, qm_angstrom, mm_angstrom):
        """Return LJ energy (Hartree), dE/dR_QM (Hartree/Bohr)."""
        from scipy.constants import Avogadro, physical_constants
        from pyscf.data.nist import BOHR
        qm, mm = np.asarray(qm_angstrom, dtype=float), np.asarray(mm_angstrom, dtype=float)
        if (qm.shape != (len(self.qm_types), 3) or mm.shape != (len(self.mm_types), 3)
                or not np.isfinite(qm).all() or not np.isfinite(mm).all()):
            raise ValueError('LJ coordinates do not match the bound QM/MM atom types')
        # Distances in nm match TPR coefficients; convert gradient at the end.
        delta = (qm[:, None, :]-mm[None, :, :])/10
        r2 = np.einsum('ijk,ijk->ij', delta, delta)
        if np.any(r2 < 1e-12):
            raise ValueError('coincident QM/MM sites in Lennard-Jones calculation')
        types = np.ix_(self.qm_types, self.mm_types)
        c6, c12 = np.asarray(self.c6_kj_mol_nm6)[types], np.asarray(self.c12_kj_mol_nm12)[types]
        inv6 = 1/r2**3
        attraction, repulsion = c6*inv6, c12*inv6**2
        kj_per_hartree = physical_constants['Hartree energy'][0]*Avogadro/1000
        energy = float((repulsion-attraction).sum()/kj_per_hartree)
        gradient_nm = np.einsum('ij,ijk->ik', (6*attraction-12*repulsion)/r2, delta)
        gradient = gradient_nm*(BOHR/10)/kj_per_hartree
        if not np.isfinite(energy) or not np.isfinite(gradient).all():
            raise ValueError('nonfinite Lennard-Jones energy/gradient')
        return energy, gradient


def parse_gromacs_dumps(parameters, system):
    """Strictly parse resolved LJ matrix and expanded global atom type indices."""
    n = re.search(r'\batnr\s*=\s*(\d+)', parameters)
    power = re.search(r'\breppow\s*=\s*([\d.eE+-]+)', parameters)
    if n is None or power is None or float(power[1]) != 12:
        raise ValueError('TPR must contain a resolved LJ table with reppow=12')
    if re.search(r'bIntermolecularInteractions\s*=\s*true', parameters, re.I):
        raise ValueError('intermolecular bonded interactions in TPR are unsupported')
    n = int(n[1])
    entries = re.findall(r'functype\[\s*(\d+)\]\s*=\s*LJ_SR,\s*c6\s*=\s*([\d.eE+-]+),\s*c12\s*=\s*([\d.eE+-]+)', parameters)
    pairs = {int(index): (float(a), float(b)) for index, a, b in entries}
    if n < 1 or set(pairs) != set(range(n*n)) or len(entries) != n*n:
        raise ValueError('incomplete/unsupported TPR nonbonded matrix; all type pairs must be LJ_SR')
    c6 = [[pairs[i*n+j][0] for j in range(n)] for i in range(n)]
    c12 = [[pairs[i*n+j][1] for j in range(n)] for i in range(n)]
    count = re.search(r'\bnatoms\s*=\s*(\d+)', system)
    atoms = {}
    for index, body in re.findall(r'atom\[\s*(\d+)\]\s*=\s*\{([^}]+)\}', system):
        fields = dict(re.findall(r'(\w+)\s*=\s*([^,}]+)', body))
        if 'type' not in fields:
            continue
        if (fields.get('ptype', '').strip() != 'Atom' or fields['type'].strip() != fields.get('typeB', '').strip()
                or float(fields['q']) != float(fields['qB'])):
            raise ValueError('virtual sites and perturbed A/B atom parameters are unsupported')
        if int(index) in atoms:
            raise ValueError('duplicate global atoms in gmx dump -sys')
        atoms[int(index)] = (int(fields['type']), int(fields['atomnumber']))
    if count is None or set(atoms) != set(range(int(count[1]))):
        raise ValueError('gmx dump -sys atom mapping is incomplete')
    values = [atoms[i] for i in range(len(atoms))]
    # Validate coefficients/types even before binding to a particular QM molecule.
    LennardJones((0,), tuple(v[0] for v in values), c6, c12, {})
    return {'c6': c6, 'c12': c12, 'types': [v[0] for v in values], 'atomic_numbers': [v[1] for v in values]}


def read_gromacs_lj(tpr, gmx='gmx'):
    import hashlib
    before = file_hash(tpr)
    dumps = []
    for options in ([], ['-sys']):
        process = subprocess.run([str(gmx), 'dump', '-s', str(tpr), *options], text=True, capture_output=True, timeout=120)
        if process.returncode:
            raise ValueError(f'gmx dump failed: {process.stderr[-1500:]}')
        dumps.append(process.stdout)
    if file_hash(tpr) != before:
        raise ValueError('TPR changed during LJ parameter extraction')
    result = parse_gromacs_dumps(*dumps)
    result['provenance'] = {'source': 'GROMACS_resolved_TPR_LJ_SR_matrix', 'tpr_sha256': before,
                            'dump_sha256': [hashlib.sha256(text.encode()).hexdigest() for text in dumps],
                            'repulsive_power': 12, 'mixing_rule_inferred': False}
    return result


def bind_gromacs_lj(parameters, record, environment):
    from pyscf.data.elements import charge
    ids = record['global_atom_indices_zero_based']
    if any(parameters['atomic_numbers'][i] != charge(e) for i, e in zip(ids, record['elements'])):
        raise ValueError('TPR LJ global atom mapping does not match QM elements')
    mm_ids = []
    for image in environment.provenance['molecular_images']:
        if image['site_start'] != len(mm_ids):
            raise ValueError('MM image mapping is not contiguous')
        mm_ids.extend(image['global_atom_indices_zero_based'])
    if len(mm_ids) != len(environment.charges_e):
        raise ValueError('MM image/site count mismatch in LJ mapping')
    atoms = list(zip(record['elements'], record['coordinates_angstrom']))
    return LennardJones(tuple(parameters['types'][i] for i in ids), tuple(parameters['types'][i] for i in mm_ids),
                        parameters['c6'], parameters['c12'],
                        {**parameters['provenance'], 'qm_initial_geometry_fingerprint': geometry_key(atoms),
                         'mm_sites_fingerprint': sites_key(environment), 'qm_global_atom_indices_zero_based': ids})


def read_lj_file(filename, atoms, environment):
    data = json.loads(Path(filename).read_text())
    if data.get('schema') != 'cube-qmmm-lj-1':
        raise ValueError('mm-vdw-file must use cube-qmmm-lj-1 schema')
    provenance = data['provenance']
    if (provenance.get('qm_initial_geometry_fingerprint') != geometry_key(atoms)
            or provenance.get('mm_sites_fingerprint') != sites_key(environment)):
        raise ValueError('LJ parameters are bound to different QM geometry/order or MM coordinates/charges')
    return LennardJones(data['qm_types'], data['mm_types'], data['c6_kj_mol_nm6'], data['c12_kj_mol_nm12'],
                        {**provenance, 'parameter_file_sha256': file_hash(filename)})
