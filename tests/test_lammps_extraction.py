"""LAMMPS text I/O, periodic image identity and independent geometry checks."""
import csv
import importlib.util
import itertools
import json
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

from wren.workflow.control import load_input
from wren.workflow.frame_geometry import load_frame, pair_geometry

EXAMPLE = Path(__file__).resolve().parents[1] / 'examples/lammps'


@unittest.skipUnless(importlib.util.find_spec('MDAnalysis'), 'workflow extra is required')
class LammpsExtractionTests(unittest.TestCase):
    def setUp(self):
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.root = Path(directory.name)
        for source in EXAMPLE.iterdir():
            if source.is_file():
                shutil.copy(source, self.root / source.name)
        self.input = self.root/'input'
        self.input.write_text(self.input.read_text().replace('../../build/lammps/mol', 'mol'))
        self.settings = load_input(self.input).sections['get_dimers_monomers']

    def extract(self):
        from wren.workflow.dimers import extract_dimers_monomers
        index = json.loads(extract_dimers_monomers(self.settings).read_text())
        return [load_frame(self.root/'mol'/f['manifest']) for f in index['frames']]

    def test_triclinic_two_frames_downstream_and_resume(self):
        import numpy as np
        from wren.dos.periodic import _full_cell
        from wren.workflow.dimers import extract_dimers_monomers
        from wren.workflow.embedding import load_manifest_environment, select_mm_sites
        frames = self.extract()
        for i, frame in enumerate(frames):
            self.assertEqual(frame['lammps_timestep'], 200*i)
            self.assertEqual(frame['time_ps'], .1*i)
            self.assertEqual(frame['dimer_count'], 1)
            self.assertEqual(frame['dimers'][0]['id'], '7_19_img_0_-1_0')
            self.assertAlmostEqual(frame['dimers'][0]['core_min_angstrom'], (.15+.01*i)*85**.5)
            self.assertEqual(frame['monomers'][0]['lammps_atom_ids'], [11, 14])
            atoms, cell = _full_cell(frame, frame['monomers'], 4)
            np.testing.assert_allclose(cell, [[10,0,0], [2,9,0], [-1,1,8]])
            np.testing.assert_allclose(atoms[0][1], [9.4,1.4,4.])
            np.testing.assert_allclose(np.subtract(atoms[1][1], atoms[0][1]), [.9,0,0])
            env = load_manifest_environment(frame)
            sites = select_mm_sites(env, frame['monomers'][0], 2.)
            with (self.root/f'mol/frame{i:04d}/clusters.csv').open() as table:
                self.assertEqual(list(csv.reader(table))[1], ['7','19'])
            offset = frame['monomers'][0]['cluster']['mm_image_offsets']['0']
            self.assertEqual(list(sites['provenance']['molecular_images'][0]['image']), offset)
            self.assertEqual(sites['charges_e'], [.4, -.4])
            self.assertEqual(sites['provenance']['molecular_images'][0]['molecule_id'], '19')
        self.assertEqual(extract_dimers_monomers(self.settings, resume=True).name, 'frames.json')
        with (self.root/'hf.data').open('a') as f:
            f.write('\n')
        with self.assertRaisesRegex(ValueError, 'changed since extraction'):
            load_manifest_environment(frames[0])

    def test_dimer_mm_lists_match_direct_all_atom_distances_in_triclinic_box(self):
        import numpy as np
        from scipy.spatial.distance import cdist
        self.settings.update(dimers_cutoff=10., **{'monomers-qmmm-cutoff': 11.})
        frames = self.extract()
        for frame in frames:
            cutoff = self.settings['monomers-qmmm-cutoff']
            lattice = np.array(frame['lattice_angstrom'])
            coords = {m['id']: np.asarray(m['coordinates_angstrom']) + m['transport_position_angstrom']
                      for m in frame['monomers']}
            with (self.root/f"mol/frame{frame['frame']:04d}/dimers.csv").open() as table:
                for row, pair in zip(csv.DictReader(table), frame['dimers'], strict=True):
                    a, b = row['dimer_qm'].split('-')
                    image = tuple(pair['image'])
                    qm = np.concatenate((coords[a], coords[b] + image @ lattice))
                    qm_sites = {(a, (0,0,0)), (b, image)}
                    expected = {(m, shift) for m, xyz in coords.items()
                                for shift in itertools.product(range(-3,4), repeat=3)
                                if (m, shift) not in qm_sites and cdist(qm, xyz + shift @ lattice).min() < cutoff}
                    listed = [(m, tuple(pair['mm_image_offsets'].get(str(i), (0,0,0))))
                              for i, m in enumerate(row['dimer_mm'].split('-'))]
                    self.assertCountEqual(expected, listed)
                    self.assertGreater(len(listed), len(coords))  # Several images of the same molecule.

    def test_cartesian_and_scaled_wrapped_and_unwrapped_with_metal_time(self):
        import numpy as np
        from wren.workflow.dimers import molecular_trajectory
        # Nonzero-origin orthogonal box. Molecule centres differ by 0.85 b;
        # inside each molecule both H-F bonds cross the a cell face.
        cases = {'xs ys zs':[(.97,.1,.5),(.06,.1,.5),(.97,.95,.5),(.06,.95,.5)],
                 'xsu ysu zsu':[(1.97,.1,.5),(2.06,.1,.5),(-.03,.95,.5),(.06,.95,.5)]}
        for scaled, cartesian in [('xs ys zs', 'x y z'), ('xsu ysu zsu', 'xu yu zu')]:
            cases[cartesian] = np.array(cases[scaled]) @ np.diag([10,9,8]) + [-3,1,-2]
        self.settings.update({'lammps-units':'metal','lammps-timestep':.002})
        for columns, rows in cases.items():
            with self.subTest(columns=columns):
                header = ('ITEM: TIMESTEP\n250\nITEM: NUMBER OF ATOMS\n4\n'
                          'ITEM: BOX BOUNDS pp pp pp\n-3 7\n1 10\n-2 6\nITEM: ATOMS id type '+columns+'\n')
                (self.root/'hf.dump').write_text(header+''.join(
                    f'{atom} {kind} '+ ' '.join(map(str,row))+'\n'
                    for atom, kind, row in zip([11,14,25,28], [1,2,1,2], rows)))
                with molecular_trajectory(self.settings) as (_, read, _, _):
                    molecules, cell, time, _ = read(0)
                    self.assertEqual(time, .5)
                    np.testing.assert_allclose(cell, np.diag([10,9,8]))
                    np.testing.assert_allclose(molecules[0].positions, [[9.7,.9,4],[10.6,.9,4]])
                    np.testing.assert_allclose(molecules[1].positions, [[9.7,8.55,4],[10.6,8.55,4]])

    def test_ambiguous_half_cell_bond(self):
        import numpy as np
        from wren.workflow.lammps import _bond_vectors
        with self.assertRaisesRegex(ValueError, 'equally short'):
            _bond_vectors(np.array([[5.,0.,0.]]), np.eye(3)*10)

    def test_diagnostic_errors(self):
        from wren.workflow.dimers import molecular_trajectory
        dump = self.root/'hf.dump'
        original = dump.read_text()
        for before, after, error in [
            ('ITEM: TIME\n0\n','', 'TIME or an explicit'),
            ('UNITS\nreal', 'UNITS\nmetal', 'unit style'),
            ('pp pp pp','pp ff pp','fully periodic'),
            ('14 7 2 0.06 0.10', '11 7 1 0.06 0.10','duplicate or unknown'),
            ('14 7 2 0.06 0.10','14 7 1 0.06 0.10','type or molecule ID')]:
            with self.subTest(error=error), self.assertRaisesRegex(ValueError,error):
                dump.write_text(original.replace(before, after))
                with molecular_trajectory(self.settings) as (_, read, _, _):
                    read(0)
        dump.write_text(original)
        data = self.root/'hf.data'
        original = data.read_text()
        for before, after, error in [('1 1 11 14','1 1 11 25','crosses molecule IDs'),
                                     ('Atoms # full','Atoms # atomic','atom_style full'),
                                     ('2 bonds','3 bonds','counts disagree')]:
            with self.subTest(error=error), self.assertRaisesRegex(ValueError,error):
                data.write_text(original.replace(before, after))
                with molecular_trajectory(self.settings):
                    pass
        data.write_text(original)
        (self.root/'elements.txt').write_text('1 H\n')
        with self.assertRaisesRegex(ValueError, 'type 2 needs'):
            with molecular_trajectory(self.settings):
                pass
        control = self.input.read_text()
        self.input.write_text(control.replace('lammps-data', 'tpr = hf.data\nlammps-data'))
        with self.assertRaisesRegex(ValueError, 'either LAMMPS data/dump or GROMACS'):
            load_input(self.input)
        self.input.write_text(control+('$reorganization-energy\ndo=1\nenvironment=qmmm\n'
                                       'mm-source=manifest\nmm-vdw=tpr\n$end\n'))
        with self.assertRaisesRegex(ValueError, 'requires a GROMACS TPR'):
            load_input(self.input)

    @unittest.skipUnless(shutil.which('gmx'), 'GROMACS grompp required for independent reader cross-check')
    def test_identical_gromacs_tpr_gro_xtc(self):
        import MDAnalysis as mda
        import numpy as np
        from wren.workflow.dimers import extract_dimers_monomers
        from wren.workflow.embedding import load_gromacs_environment, select_mm_sites
        # Independent input in nm. Atom order is H,F in each molecule.
        xyz_nm = [[.940,.140,.400],[.030,.140,.400],[1.110,.905,.400],[.200,.905,.400]]
        gro = 'Two HF molecules\n4\n'+''.join(
            f'{i//2+1:5d}{"MOL":<5}{symbol:>5}{i+1:5d}{x:8.3f}{y:8.3f}{z:8.3f}\n'
            for i,(symbol,(x,y,z)) in enumerate(zip(['H','F','H','F'],xyz_nm)))
        (self.root/'hf.gro').write_text(gro+'1.0 0.9 0.8 0 0 0.2 0 -0.1 0.1\n')
        (self.root/'hf.top').write_text("""[ defaults ]
1 2 yes 0.5 0.5
[ atomtypes ]
H 1 1.008 0 A 0.1 0.1
F 9 18.9984 0 A 0.1 0.1
[ moleculetype ]
MOL 1
[ atoms ]
1 H 1 MOL H 1 0.4 1.008
2 F 1 MOL F 1 -0.4 18.9984
[ bonds ]
1 2 1 0.09 100
[ system ]
small box
[ molecules ]
MOL 2
""")
        (self.root/'hf.mdp').write_text('integrator=md\nnsteps=0\ndt=.001\nrlist=.2\nrcoulomb=.2\nrvdw=.2\n')
        result = subprocess.run(['gmx','grompp','-f','hf.mdp','-p','hf.top','-c','hf.gro','-o','hf.tpr'],
                                cwd=self.root, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        reference = mda.Universe(str(self.root/'hf.gro'), dt=1.)
        with mda.Writer(str(self.root/'hf.xtc'), n_atoms=4) as writer:
            writer.write(reference)
        reference.trajectory.close()
        settings = {k:str(self.root/f'hf.{k}') for k in ('tpr','gro','xtc')}
        settings.update({'frame':0, 'fragment-file':self.settings['fragment-file'],
                         'output':str(self.root/'gromacs'), 'format':['xyz'],
                         'dimers_cutoff':2., 'monomers-qmmm-cutoff':2.})
        gromacs = load_frame(extract_dimers_monomers(settings))
        lammps = self.extract()[0]
        ga, la = [pair_geometry(f['dimers'][0], {m['id']:m for m in f['monomers']}, f['lattice_angstrom'])
                  for f in (gromacs, lammps)]
        self.assertEqual(ga['image'], la['image'])
        for key in ('core_min_angstrom','core_com_angstrom','coordinates_angstrom'):
            np.testing.assert_allclose(ga[key], la[key], rtol=0, atol=2e-5)
        frame = load_gromacs_environment(gromacs)
        gm = select_mm_sites(frame, gromacs['monomers'][0], 2.)
        from wren.workflow.embedding import load_manifest_environment
        lm = select_mm_sites(load_manifest_environment(lammps), lammps['monomers'][0], 2.)
        np.testing.assert_allclose(gm['charges_e'], lm['charges_e'], atol=1e-7)
        np.testing.assert_allclose(gm['coordinates_angstrom'], lm['coordinates_angstrom'], atol=2e-5)
