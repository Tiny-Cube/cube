"""Write and read the Wren fixed-topology HDF5 extraction schema."""
import errno
import hashlib
import json
import os
from pathlib import Path
import tempfile

import h5py
import numpy as np

from .dimers import (contact_images, dimer_shell, fingerprint, read_core, selected_frames,
                     validate_fragments, _mm_shells)

SCHEMA = 'wren-geometry/1'
FILENAME = 'geometry.wren.h5'


def archive_path(output):
    path = Path(output)
    return path if path.suffix.lower() in ('.h5', '.hdf5') else path/FILENAME


def _topology(h, molecules, ids, atom_ids, bonds):
    top = h.create_group('topology')
    indices = np.concatenate([m.indices for m in molecules])
    offsets = np.cumsum([0, *map(len, molecules)])
    top['molecule_id'] = np.asarray(ids, dtype='i8')
    top['molecule_offset'] = offsets
    top['global_atom_index'] = indices
    top['atom_id'] = (np.asarray(atom_ids)[indices] if atom_ids is not None else
                      np.concatenate([getattr(m, 'ids', m.indices) for m in molecules]))
    for name, values in [('element', [m.elements for m in molecules]), ('atom_name', [m.names for m in molecules])]:
        top[name] = np.asarray([str(v).encode() for v in np.concatenate(values)])
    top['charge_e'] = np.concatenate([m.charges for m in molecules])
    original = np.asarray(bonds, dtype='i8').reshape(-1, 2)
    top['bonds'] = np.argsort(indices)[original]
    return offsets[-1]


def _array(group, name, shape, dtype, chunks=True):
    return group.create_dataset(name, (0, *shape), dtype=dtype, maxshape=(None, *shape),
                                chunks=chunks, compression='gzip', compression_opts=4, shuffle=True,
                                fletcher32=True)


def _append(ds, rows):
    array = np.asarray(rows, dtype=ds.dtype).reshape((-1, *ds.shape[1:]))
    if len(array):
        start = len(ds)
        ds.resize(start+len(array), axis=0)
        ds[start:] = array


def _neighbors(h, local, cutoff, mm_cutoff):
    group = h.create_group('neighbors')
    group['frame_id'] = h['frames/id']  # The same dataset: every saved frame has both tables.
    group['core_local_atom_numbers_one_based'] = np.asarray(local)+1
    group.attrs.update(comparison='strictly_less_than', dimer_core_cutoff_angstrom=cutoff,
                       mm_all_atom_cutoff_angstrom=mm_cutoff)
    for name in ('monomer', 'dimer'):
        table = group.create_group(name)
        for key in ('frame_ptr', 'mm_ptr'):
            _append(_array(table, key, (), 'i8'), [0])
        _array(table, 'qm_mol_id', (2,) if name=='dimer' else (), 'i8')
        _array(table, 'mm_mol_id', (), 'i8')
        _array(table, 'mm_image', (3,), 'i4')
        _array(table, 'mm_min_angstrom', (), 'f8')
        if name=='dimer':
            _array(table, 'qm_b_image', (3,), 'i4')
            _array(table, 'core_min_angstrom', (), 'f8')
            _array(table, 'core_com_angstrom', (), 'f8')


def _write_rows(table, qm_ids, shells, molecule_ids):
    """Append CSR rows for one frame; offsets are absolute across the archive."""
    mm_ids, images, distances, ends = [], [], [], []
    offset = len(table['mm_mol_id'])
    for shell in shells:
        for (j, image), distance in sorted(shell.items()):
            mm_ids.append(molecule_ids[j])
            images.append(image)
            distances.append(distance)
        ends.append(offset+len(mm_ids))
    _append(table['qm_mol_id'], qm_ids)
    _append(table['mm_mol_id'], mm_ids)
    _append(table['mm_image'], images)
    _append(table['mm_min_angstrom'], distances)
    _append(table['mm_ptr'], ends)
    _append(table['frame_ptr'], [len(table['qm_mol_id'])])


def _write_frame(h, molecules, lattice, local, ids, cutoff, mm_cutoff):
    xyz = np.concatenate([m.positions for m in molecules])
    _append(h['geometry/positions'], [xyz])
    positions = [np.asarray(m.positions, dtype=float) for m in molecules]
    cores = [p[local] for p in positions]
    shells = _mm_shells(positions, lattice, mm_cutoff)
    _write_rows(h['neighbors/monomer'], ids, shells, ids)
    pairs = list(contact_images(cores, lattice, cutoff))
    groups = (dimer_shell(shells, i, j, tuple(map(int, image))) for i,j,image,_ in pairs)
    table = h['neighbors/dimer']
    _write_rows(table, [(ids[i],ids[j]) for i,j,_,_ in pairs], groups, ids)
    _append(table['qm_b_image'], [image for _,_,image,_ in pairs])
    _append(table['core_min_angstrom'], [d for _,_,_,d in pairs])
    masses = [m.masses[local] for m in molecules]
    if any(np.any(w <= 0) for w in masses):
        raise ValueError('positive core masses are required for core COM reporting')
    centers = [np.average(p, axis=0, weights=w) for p,w in zip(cores,masses)]
    _append(table['core_com_angstrom'], [np.linalg.norm(centers[j]+image@lattice-centers[i]) for i,j,image,_ in pairs])


def _digest(h):
    """Hash logical content, including attributes, in bounded batches."""
    sha = hashlib.sha256()
    def visit(name, obj):
        sha.update(name.encode()+b'\0')
        for key, value in sorted(obj.attrs.items()):
            if key != 'content_sha256':
                sha.update(key.encode()+b'\0'+str(value).encode()+b'\0')
        if isinstance(obj, h5py.Dataset):
            sha.update(str((obj.shape,obj.dtype.str)).encode())
            batch = max(1, (4<<20)//max(1, obj.dtype.itemsize*int(np.prod(obj.shape[1:]))))
            for i in range(0,len(obj),batch):
                sha.update(obj[i:i+batch].tobytes())
    visit('/', h)
    h.visititems(visit)
    return sha.hexdigest()


def _publish(staged, dest):
    """Publish a closed archive without replacing another writer's result."""
    try:
        os.link(staged, dest)
    except OSError as exc:
        if exc.errno not in (errno.EPERM, errno.EOPNOTSUPP, errno.ENOSYS):
            raise
        if os.name == 'nt':
            os.rename(staged, dest)  # On Windows this fails if dest exists.
            return
        # Linux removable filesystems can support exclusive rename without hard links.
        import ctypes
        rename = getattr(ctypes.CDLL(None, use_errno=True), 'renameat2', None)
        if rename is None:
            raise
        at_fdcwd, rename_noreplace = -100, 1
        if rename(at_fdcwd, os.fsencode(staged), at_fdcwd, os.fsencode(dest), rename_noreplace):
            error = ctypes.get_errno()
            raise OSError(error, os.strerror(error), dest) from exc


def write_archive(settings, provenance, trajectory, resume):
    count, read, labels, atom_ids, bonds = trajectory
    frames = selected_frames(settings['frame'], settings.get('frame_step',1), count)
    dest = archive_path(settings['output'])
    key = fingerprint({'schema':SCHEMA, 'provenance':provenance, 'frames':frames})
    if dest.exists():
        if not resume:
            raise ValueError(f'output exists: {dest}; use --resume or choose another output')
        with h5py.File(dest) as h:
            if (h.attrs.get('fingerprint') != key or h.attrs.get('schema') != SCHEMA
                    or 'neighbors/frame_id' not in h or 'frames/id' not in h
                    or h['neighbors/frame_id'].id != h['frames/id'].id
                    or h.attrs.get('content_sha256') != _digest(h)):
                raise ValueError('cannot resume: archive or extraction inputs changed')
        return dest
    cutoff = settings.get('dimers_cutoff', settings.get('cutoff',6.))
    mm_cutoff = settings.get('monomers-qmmm-cutoff',12.)
    local = read_core(settings['fragment-file'])
    dest.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.wren-hdf5-', dir=dest.parent) as tmp:
        staged = Path(tmp)/dest.name
        with h5py.File(staged, 'w') as h:
            h.attrs.update(schema=SCHEMA, coordinate_unit='angstrom', time_unit='ps', cell_vectors='rows',
                           image_rule='xyz + image @ cell', fingerprint=key, trajectory_frame_count=count,
                           provenance=json.dumps(provenance, sort_keys=True))
            f = h.create_group('frames')
            for name,shape,dtype in [('id',(),'i8'),('time_ps',(),'f8'),('cell',(3,3),'f8'),
                                     ('origin',(3,),'f8'),('pbc',(3,),'?')]:
                _array(f, name, shape, dtype)
            if labels is not None:
                _array(f, 'lammps_timestep', (), 'i8')
            for row, frame in enumerate(frames):
                print(f'wren: extracting HDF5 frame {frame}', flush=True)
                molecules, cell, time, details = read(frame)
                if row == 0:
                    validate_fragments(molecules, local)
                    ids = list(range(1,len(molecules)+1)) if labels is None else labels
                    n = _topology(h, molecules, ids, atom_ids, bonds())
                    _array(h.create_group('geometry'), 'positions', (n,3), molecules[0].positions.dtype,
                           chunks=(1,min(n,8*len(molecules[0])),3))
                    _neighbors(h, local, cutoff, mm_cutoff)
                _append(f['id'], [frame])
                _append(f['time_ps'], [time])
                _append(f['cell'], [cell])
                _append(f['origin'], [details.get('lammps_box_origin_angstrom',(0.,0.,0.))])
                _append(f['pbc'], [[True]*3])
                if labels is not None:
                    _append(f['lammps_timestep'], [details['lammps_timestep']])
                _write_frame(h, molecules, cell, local, ids, cutoff, mm_cutoff)
            h.attrs['content_sha256'] = _digest(h)
        # Publish only after the complete archive and its checksum are closed.
        _publish(staged, dest)
    return dest


def _open(path):
    h = h5py.File(path, 'r')
    if h.attrs.get('schema') != SCHEMA or h.attrs.get('coordinate_unit') != 'angstrom':
        h.close()
        raise ValueError(f'expected a {SCHEMA} Angstrom archive: {path}')
    return h


def frame_summary(path):
    with _open(path) as h:
        ids = h['frames/id'][:]
        counts = np.diff(h['neighbors/dimer/frame_ptr'][:])
        return [{'frame':int(f), 'time_ps':float(t), 'molecule_count':len(h['topology/molecule_id']),
                 'dimer_count':int(d)} for f,t,d in zip(ids,h['frames/time_ps'][:],counts,strict=True)]


def _row(values, wanted, label):
    found = np.flatnonzero(values[:] == wanted)
    if len(found) != 1:
        raise ValueError(f'unknown or duplicated {label} {wanted} in archive')
    return int(found[0])


def monomer_xyz(path, frame, molecule):
    with _open(path) as h:
        row = _row(h['frames/id'],frame,'frame')
        idx = _row(h['topology/molecule_id'],molecule,'molecule')
        a,b = h['topology/molecule_offset'][idx:idx+2]
        xyz = h['geometry/positions'][row,a:b]
        elements = h['topology/element'].asstr()[a:b]
    lines = [str(len(xyz)), f'frame{frame}.mol{molecule}; Angstrom']
    lines.extend(f'{e} {float(x):.17g} {float(y):.17g} {float(z):.17g}' for e,(x,y,z) in zip(elements,xyz))
    return '\n'.join(lines)+'\n'


def read_frame(path, frame=None):
    """Expose one archived geometry to the existing QM, QM/MM, DOS and KMC interfaces."""
    with _open(path) as h:
        f, top = h['frames'], h['topology']
        if frame is None:
            if len(f['id']) != 1:
                raise ValueError('choose a frame from the multi-frame HDF5 archive')
            frame = int(f['id'][0])
        row = _row(f['id'],frame,'frame')
        xyz = h['geometry/positions'][row].astype(float)
        ids, offsets = top['molecule_id'][:], top['molecule_offset'][:]
        elements, indices = top['element'].asstr()[:], top['global_atom_index'][:]
        cutoff = float(h['neighbors'].attrs['mm_all_atom_cutoff_angstrom'])
        lammps = 'lammps_timestep' in f
        atom_ids = top['atom_id'][:] if lammps else None
        monomers = []
        for identity,a,b in zip(ids,offsets[:-1],offsets[1:],strict=True):
            center = xyz[a:b].mean(axis=0)
            record = {'id':str(identity),'elements':elements[a:b].tolist(),
                      'coordinates_angstrom':(xyz[a:b]-center).tolist(),
                      'transport_position_angstrom':center.tolist(), 'global_atom_indices_zero_based':indices[a:b].tolist(),
                      'cluster':{'cutoff_angstrom':cutoff}}
            if lammps:
                record.update(lammps_atom_ids=atom_ids[a:b].tolist(), lammps_molecule_id=int(identity))
            monomers.append(record)
        sizes = dict(zip(map(int,ids),map(int,np.diff(offsets))))
        pairs = []
        g = h['neighbors/dimer']
        a,b = g['frame_ptr'][row:row+2]
        for (i,j),shift,d,com in zip(g['qm_mol_id'][a:b],g['qm_b_image'][a:b],g['core_min_angstrom'][a:b],
                                    g['core_com_angstrom'][a:b],strict=True):
            name = f'{i}_{j}' + ('_img_'+'_'.join(map(str,shift)) if any(shift) else '')
            pairs.append({'id':name,'mol_i':str(i),'mol_j':str(j),'image':shift.tolist(),'atom_count_a':sizes[i],
                          'core_min_angstrom':float(d),'core_com_angstrom':float(com)})
        data = {'schema':SCHEMA,'coordinate_unit':'Angstrom','frame':frame,'time_ps':float(f['time_ps'][row]),
                'frame_count':int(h.attrs['trajectory_frame_count']),'lattice_angstrom':f['cell'][row].tolist(),
                'core_atom_numbers_one_based':h['neighbors/core_local_atom_numbers_one_based'][:].tolist(),
                'cutoff_angstrom':float(h['neighbors'].attrs['dimer_core_cutoff_angstrom']),
                'provenance':json.loads(h.attrs['provenance']), 'molecule_count':len(monomers),'dimer_count':len(pairs),
                'monomers':monomers,'dimers':pairs}
        if lammps:
            data.update(lammps_timestep=int(f['lammps_timestep'][row]),lammps_box_origin_angstrom=f['origin'][row].tolist())
        return data
