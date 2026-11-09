"""LJ units/derivatives, TPR table parsing and total QM/MM energy consistency."""
from dataclasses import replace
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

try:
    import numpy as np
    import pyscf
except ImportError:
    np = pyscf = None


@unittest.skipIf(np is None or pyscf is None, 'quantum dependencies required')
class LennardJonesTests(unittest.TestCase):
    def model(self):
        from universal_kmc.workflow.lj import LennardJones
        return LennardJones((0, 1), (1, 0), ((.002, .001), (.001, .003)),
                            ((2e-6, 1e-6), (1e-6, 4e-6)), {})

    def test_lj_units_finite_difference_and_translation(self):
        from scipy.constants import Avogadro, physical_constants
        from pyscf.data.nist import BOHR
        from universal_kmc.workflow.lj import LennardJones
        qm = np.array([[0., 0, 0], [.8, .2, 0]])
        mm = np.array([[3.5, .1, .3], [-3., .5, .1]])
        model = self.model()
        e, grad = model.evaluate(qm, mm)
        kj_per_hartree = physical_constants['Hartree energy'][0]*Avogadro/1000
        expected = sum(model.c12_kj_mol_nm12[a][b]/(np.linalg.norm(qm[i]-mm[j])/10)**12
                       -model.c6_kj_mol_nm6[a][b]/(np.linalg.norm(qm[i]-mm[j])/10)**6
                       for i, a in enumerate(model.qm_types) for j, b in enumerate(model.mm_types))/kj_per_hartree
        self.assertAlmostEqual(e, expected, places=14)
        for i in range(2):
            for axis in range(3):
                plus, minus = qm.copy(), qm.copy()
                plus[i, axis] += 1e-4*BOHR
                minus[i, axis] -= 1e-4*BOHR
                finite = (model.evaluate(plus, mm)[0]-model.evaluate(minus, mm)[0])/2e-4
                self.assertAlmostEqual(grad[i, axis], finite, places=8)
        moved, gradient = model.evaluate(qm+[4, -5, 3], mm+[4, -5, 3])
        self.assertAlmostEqual(e, moved, places=14)
        np.testing.assert_allclose(grad, gradient, atol=1e-12)
        for types in ((0, 2), (0., 1.), (-1, 0)):
            with self.assertRaises(ValueError):
                replace(model, qm_types=types)
        with self.assertRaises(ValueError):
            replace(model, c6_kj_mol_nm6=((.002, .003), (.001, .003)))
        with self.assertRaisesRegex(ValueError, 'coincident'):
            model.evaluate(qm, qm)

    def test_resolved_tpr_pair_table_and_atom_mapping_validation(self):
        from universal_kmc.workflow.lj import parse_gromacs_dumps
        parameter = '''atnr=2
functype[0]=LJ_SR, c6=1e-3, c12=1e-6
functype[1]=LJ_SR, c6=2e-3, c12=2e-6
functype[2]=LJ_SR, c6=2e-3, c12=2e-6
functype[3]=LJ_SR, c6=3e-3, c12=3e-6
reppow = 12
'''
        system = '''natoms = 2
atom[ 0]={type= 0, typeB= 0, ptype= Atom, q=0.1, qB=0.1, atomnumber=6}
atom[ 1]={type= 1, typeB= 1, ptype= Atom, q=-0.1, qB=-0.1, atomnumber=1}
'''
        result = parse_gromacs_dumps(parameter, system)
        self.assertEqual(result['types'], [0, 1])
        self.assertEqual(result['c6'], [[.001, .002], [.002, .003]])
        for p, s in ((parameter.replace('reppow = 12', 'reppow = 9'), system),
                     (parameter.replace('functype[2]=LJ_SR', 'functype[2]=BHAM'), system),
                     (parameter, system.replace('typeB= 1', 'typeB= 0')),
                     (parameter, system.replace('qB=0.1', 'qB=0.2')),
                     (parameter, system.replace('natoms = 2', 'natoms = 3'))):
            with self.assertRaises(ValueError):
                parse_gromacs_dumps(p, s)

    def test_optimizer_total_gradient_and_single_point_energy_include_lj(self):
        from pyscf.data.nist import BOHR
        from universal_kmc.workflow.control import load_input, scf_config
        from universal_kmc.workflow.embedding import read_charges
        from universal_kmc.workflow.lj import LennardJones, geometry_key, sites_key, read_lj_file
        from universal_kmc.workflow.reorganization import _optimize
        from universal_kmc.workflow.scf import run_single_point
        atoms = [('H', (-.4, 0., 0.)), ('F', (.5, 0., 0.))]
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root/'HF.xyz').write_text('2\nHF\nH -.4 0 0\nF .5 0 0\n')
            (root/'mm.pc').write_text('2 3 0 .2\n-1 -3 0 -.1\n')
            environment = read_charges(root/'mm.pc')
            lj = replace(self.model(), provenance={'qm_initial_geometry_fingerprint': geometry_key(atoms),
                                                  'mm_sites_fingerprint': sites_key(environment)})
            (root/'lj.json').write_text(json.dumps(lj.payload()))
            (root/'input').write_text('''$reorganization-energy
xyz = HF.xyz
qc-engine = pyscf
method = hf
basis-set = sto-3g
conv-tol = 1e-11
environment = qmmm
mm-charges = mm.pc
mm-vdw = file
mm-vdw-file = lj.json
$end
''')
            settings = load_input(root/'input').sections['reorganization-energy']
            environment = replace(environment, vdw=read_lj_file(root/'lj.json', atoms, environment))
            measured = {}
            def capture(adapter, **kwargs):
                measured['energy'], measured['gradient'] = adapter.nuc_grad_method().as_scanner()(adapter.mol)
                return True, adapter.mol
            with patch('pyscf.geomopt.geometric_solver.kernel', side_effect=capture):
                _optimize(atoms, scf_config(settings), 0, root, 5, embedding=environment)
            energies = []
            for i, offset in enumerate((0., 1e-3, -1e-3)):
                displaced = [atoms[0], ('F', (.5+offset*BOHR, 0., 0.))]
                point = run_single_point(displaced, settings, root/f'point{i}', embedding=environment)
                energies.append(point['energy_hartree'])
                self.assertAlmostEqual(point['energy_hartree'], point['electronic_energy_hartree']+point['qm_mm_vdw_energy_hartree'])
            self.assertAlmostEqual(measured['energy'], energies[0], places=9)
            self.assertLess(abs(measured['gradient'][1, 0]-(energies[1]-energies[2])/.002), 2e-5)
            with self.assertRaisesRegex(ValueError, 'different QM geometry'):
                read_lj_file(root/'lj.json', [('H', (0., 0., 0.)), atoms[1]], environment)
            self.assertTrue(environment.payload()['qm_mm_vdw_included'])


if __name__ == '__main__':
    unittest.main()
