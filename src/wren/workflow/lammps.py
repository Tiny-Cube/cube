"""Molecular LAMMPS data and native text custom dumps (real or metal units)."""
from dataclasses import dataclass
from itertools import product
from pathlib import Path

import numpy as np


@dataclass
class Molecule:
    indices: np.ndarray
    elements: np.ndarray
    names: np.ndarray
    masses: np.ndarray
    charges: np.ndarray
    positions: np.ndarray

    def __len__(self):
        return len(self.indices)


def _element_map(path):
    from MDAnalysis.guesser.tables import SYMB2Z
    elements = {}
    for raw in Path(path).read_text().splitlines():
        row = raw.partition('#')[0].split()
        if not row:
            continue
        if len(row) != 2 or row[1] not in SYMB2Z:
            raise ValueError('element-map requires numeric atom_type and chemical symbol on each line')
        atom_type = int(row[0])
        if atom_type <= 0 or atom_type in elements:
            raise ValueError('element-map needs distinct positive atom types')
        elements[atom_type] = row[1]
    return elements


def _topology(settings):
    mapping = _element_map(settings['element-map'])
    atoms, masses, bonds, counts = {}, {}, {}, {}
    section = None
    with open(settings['lammps-data']) as data:
        next(data, None)  # mandatory title line
        for raw in data:
            line, _, comment = raw.partition('#')
            row = line.split()
            if not row:
                continue
            if row[0][0].isalpha():
                section = ' '.join(row)
                if section == 'Atoms' and comment.strip() and comment.split()[0] != 'full':
                    raise ValueError('LAMMPS Atoms requires atom_style full')
            elif section is None and len(row) == 2 and row[1] in ('atoms', 'bonds'):
                counts[row[1]] = int(row[0])
            elif section == 'Masses':
                if len(row) != 2 or int(row[0]) in masses:
                    raise ValueError('invalid or duplicate LAMMPS Masses row')
                masses[int(row[0])] = float(row[1])
            elif section == 'Atoms':
                if len(row) not in (7, 10):
                    raise ValueError('LAMMPS Atoms requires full atom-ID molecule-ID type q x y z [ix iy iz]')
                atom, molecule, kind = map(int, row[:3])
                if min(atom, molecule, kind) <= 0 or atom in atoms:
                    raise ValueError('LAMMPS atoms require unique positive atom IDs and positive molecule IDs/types')
                if not np.isfinite([float(v) for v in row[3:7]]).all():
                    raise ValueError('nonfinite LAMMPS atom charge/coordinates')
                atoms[atom] = (molecule, kind, float(row[3]))
            elif section == 'Bonds':
                if len(row) != 4:
                    raise ValueError('invalid LAMMPS Bonds row')
                bond, _, a, b = map(int, row)
                if bond in bonds or a == b:
                    raise ValueError('duplicate or self LAMMPS bond')
                bonds[bond] = (a, b)
    if not atoms or counts.get('atoms') != len(atoms) or counts.get('bonds', 0) != len(bonds):
        raise ValueError('LAMMPS data atom/bond counts disagree with records')
    ids = sorted(atoms)
    index = {atom: i for i, atom in enumerate(ids)}
    groups, adjacency = {}, [[] for _ in ids]
    for i, atom in enumerate(ids):
        mol, kind, _ = atoms[atom]
        if kind not in mapping or not np.isfinite(masses.get(kind, np.nan)) or masses[kind] <= 0:
            raise ValueError(f'LAMMPS type {kind} needs an element-map entry and positive Masses value')
        groups.setdefault(mol, []).append(i)
    for a, b in bonds.values():
        if a not in atoms or b not in atoms or atoms[a][0] != atoms[b][0]:
            raise ValueError('LAMMPS bond has unknown atoms or crosses molecule IDs')
        adjacency[index[a]].append(index[b])
        adjacency[index[b]].append(index[a])
    molecules, trees = [], []
    for mol in sorted(groups):
        indices = groups[mol]
        order, seen, tree = [indices[0]], {indices[0]}, []
        for a in order:
            for b in adjacency[a]:
                if b not in seen:
                    order.append(b)
                    seen.add(b)
                    tree.append((a, b))
        if len(seen) != len(indices):
            raise ValueError(f'LAMMPS molecule {mol} must be one bonded fragment')
        kinds = [atoms[ids[i]][1] for i in indices]
        molecules.append(Molecule(np.array(indices), np.array([mapping[t] for t in kinds]),
            np.array(kinds), np.array([masses[t] for t in kinds]),
            np.array([atoms[ids[i]][2] for i in indices]), np.empty((0, 3))))
        trees.extend(tree)
    return molecules, trees, ids, atoms, sorted(groups), [(index[a], index[b]) for a,b in bonds.values()]


def _box(stream, header):
    labels = header.split()[3:]
    tilted = labels[:3] == ['xy', 'xz', 'yz']
    if labels != (['xy', 'xz', 'yz'] if tilted else []) + ['pp']*3:
        raise ValueError('LAMMPS requires a fully periodic orthogonal or restricted triclinic BOX BOUNDS')
    rows = [list(map(float, stream.readline().split())) for _ in range(3)]
    if any(len(row) != 2 + tilted for row in rows) or not np.isfinite(rows).all():
        raise ValueError('invalid LAMMPS BOX BOUNDS rows')
    origin = np.array([row[0] for row in rows])
    high = np.array([row[1] for row in rows])
    xy, xz, yz = [row[2] for row in rows] if tilted else (0., 0., 0.)
    origin -= [min(0., xy, xz, xy+xz), min(0., yz), 0.]
    high -= [max(0., xy, xz, xy+xz), max(0., yz), 0.]
    lx, ly, lz = high-origin
    if min(lx, ly, lz) <= 0:
        raise ValueError('nonpositive LAMMPS cell length')
    return origin, np.array([[lx, 0., 0.], [xy, ly, 0.], [xz, yz, lz]])


def _bond_vectors(delta, lattice):
    """Find shortest periodic representations, including tilted boxes."""
    inv = np.linalg.inv(lattice)
    fractional = delta @ inv
    residual = fractional - np.rint(fractional)
    best = residual @ lattice
    distances = np.linalg.norm(best, axis=1)
    # For a better candidate the fractional difference on each axis cannot
    # exceed its reciprocal-vector norm times the initial separation.
    bounds = np.floor(.5 + distances.max()*np.linalg.norm(inv, axis=0)).astype(int)
    if np.prod(2*bounds.astype(object)+1) > 100000:
        raise ValueError('LAMMPS box too ill-conditioned for molecular unwrapping')
    ambiguous = np.zeros(len(delta), dtype=bool)
    for shift in product(*(range(-n, n+1) for n in bounds)):
        if shift == (0, 0, 0):
            continue
        candidate = (residual+shift) @ lattice
        norm = np.linalg.norm(candidate, axis=1)
        tied = np.isclose(norm, distances, rtol=1e-13, atol=1e-12)
        change = (norm < distances) & ~tied
        ambiguous = (ambiguous & ~change) | tied
        best[change], distances[change] = candidate[change], norm[change]
    if ambiguous.any():
        raise ValueError('LAMMPS bond has ambiguous equally short periodic images')
    return best


class LammpsTrajectory:
    def __init__(self, settings):
        self.settings = settings
        self.molecules, self.tree, self.ids, self.atoms, self.mol_ids, self.bonds = _topology(settings)
        self.index = {atom: i for i, atom in enumerate(self.ids)}
        self.offsets = []
        self.stream = open(settings['lammps-dump'])
        try:
            self._index()
        except BaseException:
            self.stream.close()
            raise

    def _index(self):
        time = None
        units = self.settings['lammps-units']
        while header := self.stream.readline():
            header = header.strip()
            if not header:
                continue
            if header == 'ITEM: UNITS':
                if self.stream.readline().strip() != units:
                    raise ValueError('LAMMPS dump unit style disagrees with lammps-units')
                continue
            if header == 'ITEM: TIME':
                if time is not None:
                    raise ValueError('duplicate LAMMPS TIME')
                time = float(self.stream.readline())
                continue
            if header != 'ITEM: TIMESTEP':
                raise ValueError(f'expected LAMMPS TIMESTEP, found {header}')
            step = int(self.stream.readline())
            if time is None:
                if 'lammps-timestep' not in self.settings:
                    raise ValueError('LAMMPS dump needs ITEM: TIME or an explicit fixed lammps-timestep')
                time = step*self.settings['lammps-timestep']
            time_ps = time*(.001 if units == 'real' else 1.)
            if not np.isfinite(time_ps):
                raise ValueError('nonfinite LAMMPS physical time')
            if self.stream.readline().strip() != 'ITEM: NUMBER OF ATOMS' or int(self.stream.readline()) != len(self.ids):
                raise ValueError('LAMMPS dump must contain all topology atoms')
            header = self.stream.readline().strip()
            if not header.startswith('ITEM: BOX BOUNDS '):
                raise ValueError('missing LAMMPS BOX BOUNDS')
            origin, lattice = _box(self.stream, header)
            header = self.stream.readline().split()
            columns = header[2:]
            if header[:2] != ['ITEM:', 'ATOMS'] or len(columns) != len(set(columns)) or not {'id','type'} <= set(columns):
                raise ValueError('LAMMPS dump requires unique columns including id and type')
            choices = [('x','y','z'), ('xs','ys','zs'), ('xu','yu','zu'), ('xsu','ysu','zsu')]
            coords = next((c for c in choices if set(c) <= set(columns)), None)
            if coords is None:
                raise ValueError('LAMMPS dump requires a complete Cartesian or scaled coordinate triplet')
            self.offsets.append((self.stream.tell(), step, time_ps, origin, lattice, columns, coords))
            for _ in self.ids:
                if not self.stream.readline():
                    raise ValueError('truncated LAMMPS dump atom rows')
            time = None
        if not self.offsets or time is not None:
            raise ValueError('empty or incomplete LAMMPS dump')

    def read(self, frame):
        offset, step, time, origin, cell, columns, triplet = self.offsets[frame]
        self.stream.seek(offset)
        col = {name: i for i, name in enumerate(columns)}
        xyz = np.empty((len(self.ids), 3))
        seen = set()
        for _ in self.ids:
            row = self.stream.readline().split()
            if len(row) != len(col):
                raise ValueError('LAMMPS dump atom row has the wrong number of fields')
            atom = int(row[col['id']])
            if atom not in self.index or atom in seen:
                raise ValueError('duplicate or unknown LAMMPS dump atom ID')
            mol, kind, charge = self.atoms[atom]
            if int(row[col['type']]) != kind or ('mol' in col and int(row[col['mol']]) != mol):
                raise ValueError('LAMMPS dump atom type or molecule ID differs from data topology')
            if 'q' in col and not np.isclose(float(row[col['q']]), charge, rtol=1e-5, atol=1e-8):
                raise ValueError('time-dependent LAMMPS charges are unsupported; dump q differs from data')
            seen.add(atom)
            xyz[self.index[atom]] = [float(row[col[key]]) for key in triplet]
        if not np.isfinite(xyz).all():
            raise ValueError('nonfinite LAMMPS dump coordinates')
        # Both supported styles use Angstrom. Keep the native origin in metadata
        # and express positions relative to it, as for a GROMACS cell.
        xyz = xyz @ cell if triplet[0] in ('xs', 'xsu') else xyz - origin
        if self.tree and triplet[0] in ('x', 'xs'):
            parent, child = np.asarray(self.tree).T
            vectors = _bond_vectors(xyz[child]-xyz[parent], cell)
            a, b = np.asarray(self.bonds).T
            wanted = _bond_vectors(xyz[b]-xyz[a], cell)
            for (p, c), v in zip(self.tree, vectors):
                xyz[c] = xyz[p]+v
            if not np.allclose(xyz[b]-xyz[a], wanted, rtol=0, atol=1e-6):
                raise ValueError('bond cycle crosses box inconsistently; a finite whole molecule is required')
        elif self.bonds:
            a, b = np.asarray(self.bonds).T
            delta = xyz[b]-xyz[a]
            if not np.allclose(delta, _bond_vectors(delta, cell), rtol=0, atol=1e-6):
                raise ValueError('unwrapped LAMMPS coordinates contain inconsistent molecular image flags')
        inv = np.linalg.inv(cell)
        for molecule in self.molecules:
            positions = xyz[molecule.indices]
            anchor = positions[0] @ inv
            nearest = np.rint(anchor)
            anchor = np.where(np.isclose(anchor, nearest, rtol=0, atol=1e-13), nearest, anchor)
            molecule.positions = positions - np.floor(anchor) @ cell
        details = {'lammps_timestep': step, 'lammps_box_origin_angstrom': origin.tolist()}
        return self.molecules, cell, time, details
