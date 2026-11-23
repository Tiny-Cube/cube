"""Join an explicitly selected finite ExtXYZ cluster, topology and optional LJ table."""
from dataclasses import dataclass, replace
import hashlib
import json
import math
from pathlib import Path
import re
import shlex

from .embedding import FrozenMM
from .lj import LennardJones, geometry_key, sites_key


ELEMENTS = frozenset((
    'H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr '
    'Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er '
    'Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm '
    'Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og').split())


@dataclass(frozen=True)
class Position:
    atom_id: int
    element: str
    xyz: tuple
    region: str


def _properties(text):
    fields = text.split(':')
    if len(fields) % 3:
        raise ValueError('ExtXYZ Properties must consist of name:type:count triples')
    columns, offset = {}, 0
    for name, kind, count in zip(fields[::3], fields[1::3], fields[2::3]):
        if not name or name in columns or kind not in ('R', 'I', 'S', 'L') or not re.fullmatch('[1-9][0-9]*', count):
            raise ValueError('invalid/duplicate ExtXYZ property')
        columns[name] = (kind, int(count), offset)
        offset += int(count)
    required = {'species': ('S', 1), 'pos': ('R', 3), 'atom_id': ('I', 1), 'region': ('S', 1), 'frozen': ('L', 1)}
    if any(columns.get(name, ())[:2] != signature for name, signature in required.items()):
        raise ValueError('ExtXYZ requires species:S:1, pos:R:3, atom_id:I:1, region:S:1 and frozen:L:1')
    return {name: columns[name][2] for name in required}, offset


def _positions(raw):
    lines = raw.decode('utf-8-sig').splitlines()
    if len(lines) < 3:
        raise ValueError('expected a single nonempty ExtXYZ frame')
    count = int(lines[0])
    if count <= 0 or len(lines) < count + 2 or any(s.strip() for s in lines[count+2:]):
        raise ValueError('expected exactly one full nonempty ExtXYZ frame')
    meta = {}
    for item in shlex.split(lines[1]):
        key, sep, value = item.partition('=')
        if sep:
            if key in meta:
                raise ValueError(f'duplicate ExtXYZ header property {key}')
            meta[key] = value
    if meta.get('coordinate_unit') != 'Angstrom':
        raise ValueError('ExtXYZ requires coordinate_unit=Angstrom')
    if ('Lattice' in meta or 'pbc' in meta) and meta.get('pbc', '').lower().split() not in (
            ['f', 'f', 'f'], ['false', 'false', 'false']):
        raise ValueError('ExtXYZ cluster must be expanded to finite coordinates with PBC disabled')
    col, width = _properties(meta.get('Properties', ''))
    positions, occupied = {}, set()
    for lineno, line in enumerate(lines[2:count+2], 3):
        values = line.split()
        if len(values) != width:
            raise ValueError(f'ExtXYZ line {lineno}: wrong column count')
        identifier = int(values[col['atom_id']])
        if identifier < 1 or identifier in positions:
            raise ValueError(f'ExtXYZ line {lineno}: atom_id must be positive and unique')
        element = values[col['species']]
        if element not in ELEMENTS:
            raise ValueError(f'ExtXYZ line {lineno}: invalid element symbol')
        xyz = tuple(map(float, values[col['pos']:col['pos']+3]))
        if not all(map(math.isfinite, xyz)):
            raise ValueError(f'ExtXYZ line {lineno}: coordinates must be finite')
        if xyz in occupied:
            raise ValueError(f'ExtXYZ line {lineno}: duplicate Cartesian site')
        occupied.add(xyz)
        region, frozen = values[col['region']].upper(), values[col['frozen']].lower()
        if region not in ('QM', 'MM') or frozen not in ('t', 'f', 'true', 'false'):
            raise ValueError(f'ExtXYZ line {lineno}: expected QM/MM and a logical frozen flag')
        if (frozen in ('t', 'true')) != (region == 'MM'):
            raise ValueError(f'ExtXYZ line {lineno}: current optimizer requires mobile QM and frozen MM atoms')
        positions[identifier] = Position(identifier, element, xyz, region)
    return positions


def _identifier(value):
    return type(value) is int and value > 0


def _real(value):
    return type(value) in (int, float) and math.isfinite(value)


def _read_json(raw):
    def unique(fields):
        record = dict(fields)
        if len(record) != len(fields):
            raise ValueError('duplicate JSON field in QM/MM input')
        return record
    return json.loads(raw, object_pairs_hook=unique)


def _topology(data, positions):
    fields = {'schema', 'charge_unit', 'atoms', 'bonds'}
    if not isinstance(data, dict) or data.keys() != fields or data['schema'] != 'wren-qmmm-topology-1' or data['charge_unit'] != 'e':
        raise ValueError('expected wren-qmmm-topology-1, charge_unit=e, atoms and bonds')
    if not isinstance(data['atoms'], list) or len(data['atoms']) != len(positions):
        raise ValueError('ExtXYZ/topology atom counts differ')
    atoms, source_images, source_parameters = {}, set(), {}
    required = {'atom_id', 'molecule_id', 'element', 'ff_type', 'charge'}
    for atom in data['atoms']:
        if (not isinstance(atom, dict) or required - atom.keys() or
                atom.keys() - (required | {'source_atom_id', 'image'})):
            raise ValueError('invalid topology atom fields')
        identifier, group = atom['atom_id'], atom['molecule_id']
        if not _identifier(identifier) or identifier in atoms or identifier not in positions:
            raise ValueError('topology atom_id must uniquely match an ExtXYZ atom_id')
        if (not isinstance(group, str) and not _identifier(group)) or not re.fullmatch('[A-Za-z0-9_-]+', str(group)):
            raise ValueError('molecule_id must be a positive integer or plain identifier')
        if atom['element'] != positions[identifier].element:
            raise ValueError(f'ExtXYZ/topology elements differ for atom_id {identifier}')
        atom_type = atom['ff_type']
        if not isinstance(atom_type, str) or not atom_type.strip() or not _real(atom['charge']):
            raise ValueError(f'invalid topology ff_type or charge for atom_id {identifier}')
        image = atom.get('image', [0, 0, 0])
        if not isinstance(image, list) or len(image) != 3 or any(type(x) is not int for x in image):
            raise ValueError('image must have three integer lattice offsets')
        source_id = atom.get('source_atom_id', identifier)
        if not _identifier(source_id) or (source_id, tuple(image)) in source_images:
            raise ValueError('source_atom_id and image must uniquely identify a site')
        source_images.add((source_id, tuple(image)))
        parameters = (str(group), atom['element'], atom_type, atom['charge'])
        if source_parameters.setdefault(source_id, parameters) != parameters:
            raise ValueError('inconsistent topology for periodic copies of a source_atom_id')
        atoms[identifier] = {**atom, 'molecule_id': str(group), 'image': tuple(image), 'source_atom_id': source_id}
    if not isinstance(data['bonds'], list):
        raise ValueError('bonds must be an explicit list of atom_id pairs')
    seen = set()
    for bond in data['bonds']:
        if not isinstance(bond, list) or len(bond) != 2 or any(not _identifier(i) or i not in atoms for i in bond):
            raise ValueError('invalid topology bond: expected two known atom IDs')
        i, j = bond
        pair = tuple(sorted(bond))
        if i == j or pair in seen:
            raise ValueError('self or duplicate topology bond')
        seen.add(pair)
        if positions[i].region != positions[j].region:
            raise ValueError('covalent bond crosses the QM/MM boundary')
    return atoms


def _lj(data, qm_types, mm_types, provenance):
    fields = {'schema', 'model', 'units', 'source', 'pairs'}
    if not isinstance(data, dict) or data.keys() != fields or data['schema'] != 'wren-qmmm-forcefield-1':
        raise ValueError('expected wren-qmmm-forcefield-1 with model, units, source and explicit pairs')
    if data['model'] != 'LJ12-6' or data['units'] != {'energy': 'kJ/mol', 'length': 'nm'}:
        raise ValueError('forcefield requires LJ12-6 and kJ/mol, nm units for C6/C12')
    if not isinstance(data['source'], str) or not data['source'].strip() or not isinstance(data['pairs'], list):
        raise ValueError('forcefield requires a parameter source and list of explicit type pairs')
    table = {}
    for entry in data['pairs']:
        if not isinstance(entry, dict) or entry.keys() != {'types', 'c6', 'c12'}:
            raise ValueError('forcefield pair requires types, c6 and c12')
        types = entry['types']
        if not isinstance(types, list) or len(types) != 2 or any(not isinstance(t, str) or not t.strip() for t in types):
            raise ValueError('invalid forcefield type pair')
        pair = tuple(sorted(types))
        numbers = (entry['c6'], entry['c12'])
        if pair in table or any(not _real(x) or x < 0 for x in numbers):
            raise ValueError(f'duplicate or invalid LJ pair {pair}')
        table[pair] = numbers
    used = {tuple(sorted((q, m))) for q in set(qm_types) for m in set(mm_types)}
    missing = used - table.keys()
    if missing:
        raise ValueError(f'missing explicit QM-MM LJ coefficients for {sorted(missing)}')
    labels = sorted(set(qm_types) | set(mm_types))
    index = {label: i for i, label in enumerate(labels)}
    c6, c12 = ([[0.] * len(labels) for _ in labels] for _ in range(2))
    for a, b in used:
        i, j = index[a], index[b]
        c6[i][j] = c6[j][i] = table[a, b][0]
        c12[i][j] = c12[j][i] = table[a, b][1]
    return LennardJones(tuple(index[t] for t in qm_types), tuple(index[t] for t in mm_types), c6, c12,
                        {**provenance, 'parameter_source': data['source'], 'type_names': labels,
                         'resolved_QM_MM_pairs': sorted(used), 'mixing_rule_inferred': False})


def read_separated_embedding(coordinates, topology, forcefield=None):
    """Return ordered QM atoms, frozen finite MM environment, central molecule ID.

    Atom identifiers bind records across files. QC and MM site ordering retain
    their relative order from the coordinate file. All MM charges are used.
    """
    paths = {'coordinates': Path(coordinates).resolve(), 'topology': Path(topology).resolve()}
    if forcefield is not None:
        paths['forcefield'] = Path(forcefield).resolve()
    raw = {k: path.read_bytes() for k, path in paths.items()}
    hashes = {k: hashlib.sha256(value).hexdigest() for k, value in raw.items()}
    positions = _positions(raw['coordinates'])
    records = _topology(_read_json(raw['topology']), positions)
    qm = [p for p in positions.values() if p.region == 'QM']
    mm = [p for p in positions.values() if p.region == 'MM']
    if not qm or not mm:
        raise ValueError('cluster must contain both QM and MM atoms')
    def group(p):
        atom = records[p.atom_id]
        return atom['molecule_id'], atom['image']
    target = group(qm[0])
    if any(group(p) != target for p in qm) or any(group(p) == target for p in mm):
        raise ValueError('QM region must be one complete molecule image')
    atoms = [(p.element, p.xyz) for p in qm]
    provenance = {'source': 'separated_extxyz_topology', 'files': {k: str(v) for k, v in paths.items()},
                  'sha256': hashes, 'reader_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  'qm_molecule_id': target[0], 'qm_image': target[1],
                  'qm_atom_ids': [p.atom_id for p in qm], 'mm_atom_ids': [p.atom_id for p in mm],
                  'selection': 'all explicitly labeled MM atoms in input coordinates'}
    environment = FrozenMM(tuple(p.xyz for p in mm), tuple(records[p.atom_id]['charge'] for p in mm), provenance)
    environment.validate_qm(atoms)
    if forcefield is not None:
        lj = _lj(_read_json(raw['forcefield']), [records[p.atom_id]['ff_type'] for p in qm],
                 [records[p.atom_id]['ff_type'] for p in mm],
                 {'forcefield_sha256': hashes['forcefield'], 'qm_initial_geometry_fingerprint': geometry_key(atoms),
                  'mm_sites_fingerprint': sites_key(environment)})
        environment = replace(environment, vdw=lj)
    return atoms, environment, target[0]
