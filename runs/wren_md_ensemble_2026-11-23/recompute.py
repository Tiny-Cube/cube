"""Recompute one real molecule on the eight reported frames from original MD files."""
import argparse
import hashlib
import json
from pathlib import Path

import MDAnalysis as mda
import numpy as np
import pyscf
from pyscf.data.nist import HARTREE2EV
from pyscf.scf import chkfile
from wren.workflow.dimers import file_hash
from wren.workflow.scf import run_single_point


def compute(args):
    paths = {'gro': args.gro, 'tpr': args.tpr, 'xtc': args.xtc}
    expected = json.loads((Path(__file__).parent/'orbital-levels.json').read_text())['sources_sha256']
    hashes = {k: file_hash(path) for k, path in paths.items()}
    if expected != hashes:
        raise ValueError('MD inputs do not match the reported frame sequence')
    settings = {'qc-engine': 'pyscf', 'basis-set': 'sto-3g', 'method': 'hf', 'auxbasis': 'weigend',
                'threads': args.threads, 'conv-tol': 1e-7, 'max-cycle': 100, 'grid-level': 3, 'max-memory-mb': 4000}
    universe = mda.Universe(str(args.tpr), str(args.xtc))
    gro = mda.Universe(str(args.gro))
    if not (np.array_equal(universe.atoms.names, gro.atoms.names)
            and np.array_equal(universe.atoms.resindices, gro.atoms.resindices)):
        raise ValueError('GRO and TPR atom ordering differs')
    molecule = universe.atoms.fragments[0]
    if len(molecule) != 86 or not np.array_equal(molecule.indices, np.arange(86)):
        raise ValueError('expected complete topology fragment 1 (86 atoms)')
    records = []
    for frame in range(1860, 2001, 20):
        ts = universe.trajectory[frame]
        molecule.unwrap(compound='fragments', reference=None)
        xyz = molecule.positions.astype(float)
        xyz -= xyz.mean(axis=0)
        atoms = list(zip(molecule.elements.tolist(), xyz.tolist()))
        status = run_single_point(atoms, settings, args.build/f'frame{frame}'/'1', resume=True)
        mol, orbitals = chkfile.load_scf(status['checkpoint'])
        if (not status['converged'] or file_hash(status['checkpoint']) != status['checkpoint_sha256']
                or not np.allclose(mol.atom_coords(unit='Angstrom'), xyz, atol=1e-8, rtol=0)):
            raise RuntimeError(f'failed to verify frame {frame} checkpoint')
        levels = np.asarray(orbitals['mo_energy'])*HARTREE2EV
        occ = np.asarray(orbitals['mo_occ'])
        records.append({'frame': frame, 'time_ps': float(ts.time), 'molecule_id': 1, 'qm_atom_count': len(atoms),
                        'coordinate_sha256': hashlib.sha256(json.dumps(atoms, separators=(',', ':')).encode()).hexdigest(),
                        'checkpoint_sha256': status['checkpoint_sha256'], 'scf_converged': True,
                        'orbital_energies_ev': levels.tolist(), 'occupations': occ.tolist(),
                        'homo_ev': float(levels[occ > 0].max()), 'lumo_ev': float(levels[occ == 0].min()),
                        'total_energy_hartree': float(orbitals['e_tot'])})
    output = {'schema': 'wren-md-sampled-orbitals-1', 'sources_sha256': hashes,
              'pyscf_version': pyscf.__version__, 'frames': records}
    (args.build/'orbital-levels.json').write_text(json.dumps(output, indent=2, allow_nan=False)+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ('gro', 'tpr', 'xtc', 'build'):
        parser.add_argument('--' + flag, required=True, type=Path)
    parser.add_argument('--threads', type=int, default=2)
    compute(parser.parse_args())
