"""Sparse KPM regression checks against dense spectra and independent quadrature."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

from wren.dos.control import normalize

HAVE_SCIPY = importlib.util.find_spec('scipy') is not None


@unittest.skipUnless(HAVE_SCIPY, 'KPM needs NumPy/SciPy')
class KernelPolynomialTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.network = {'schema': 'wren-dos-tb-1', 'coupling_convention': 'signed_real_consistent_orbital_phases',
                        'site_energies_ev': [.4, -.3, .1], 'pairs': [[0, 1, [1, 0, 0], .13],
                        [1, 2, [0, 0, 0], -.14], [0, 1, [0, 0, 0], -.21], [2, 2, [1, 0, 0], -.12]]}
        self.path = self.root/'network.json'
        self.path.write_text(json.dumps(self.network))

    def settings(self, **kwargs):
        return normalize({'mode': 'tight-binding', 'solver': 'kpm', 'input': self.path,
                          'output': 'out', 'points': 73, 'kpm-moments': 25, 'spin-degeneracy': 2,
                          'kmesh': [3, 1, 1], **kwargs}, self.root)

    def test_bloch_moments_and_bin_integrals_match_dense_phased_spectra(self):
        from contextlib import redirect_stdout
        import io
        import numpy as np
        from scipy.integrate import quad
        from wren.workflow.cli import main
        p = self.root/'input'
        p.write_text('$dos\ndo=1\nmode=tight-binding\nsolver=kpm\ninput=network.json\noutput=out\n'
                     'kpm-trace=exact\nkpm-moments=25\nkmesh=3,1,1\nspin-degeneracy=2\npoints=73\n$end\n')
        with redirect_stdout(io.StringIO()):
            self.assertEqual(main([str(p)]), 0)
            self.assertEqual(main([str(p), '--resume']), 0)
        result = json.loads((self.root/'out/result.json').read_text())
        a, b = result['kpm']['scale_ev'], result['kpm']['center_ev']
        traces = np.genfromtxt(self.root/'out/moments.csv', delimiter=',', names=True)
        bins = np.genfromtxt(self.root/'out/dos.csv', delimiter=',', names=True)
        order = np.arange(25)
        moments = np.zeros(len(order))
        for k in range(3):
            phase = np.exp(2j*np.pi*k/3)
            offdiag = -.21+.13*phase
            h = [[.4, offdiag, 0.], [offdiag.conjugate(), -.3, -.14], [0., -.14, .1-.24*phase.real]]
            angle = np.arccos((np.linalg.eigvalsh(h)-b)/a)
            moments += 2/3*np.cos(np.outer(order, angle)).sum(axis=1)
        np.testing.assert_allclose(traces['trace_moment'], moments, atol=2e-13)
        phi = np.pi/26
        damp = ((26-order)*np.cos(order*phi)+np.sin(order*phi)/np.tan(phi))/26
        def spectral_weight(theta):
            return (moments[0]+2*np.sum(damp[1:]*moments[1:]*np.cos(order[1:]*theta)))/np.pi
        for index in (0, 8, 36, 71, 72):
            lo, hi = bins['left_ev'][index], bins['right_ev'][index]
            angles = np.arccos(np.clip((np.array([hi, lo])-b)/a, -1, 1))
            mass = quad(spectral_weight, *angles)[0]
            self.assertAlmostEqual(bins['dos'][index]*(hi-lo), mass, places=12)
        self.assertAlmostEqual(float(bins['dos'] @ (bins['right_ev']-bins['left_ev'])), 6.)
        self.assertEqual(result['frame_count'], 1)
        self.assertFalse((self.root/'out/spectra.json').exists())

    def test_seeded_common_probes_errors_and_cropped_window(self):
        import numpy as np
        from wren.dos.runner import run
        s = self.settings(**{'kpm-vectors': 8, 'kpm-seed': 17, 'emin': -.12, 'emax': .3})
        report = run(s)
        out = self.root/'out'
        first = (out/'dos.csv').read_bytes()
        data = np.genfromtxt(out/'dos.csv', names=True, delimiter=',')
        self.assertGreater(data['trace_standard_error'].max(), 0.)
        self.assertGreater(report['finite_grid_integral'], 0.)
        self.assertLess(report['finite_grid_integral'], 6.)
        run(dict(s, output=str(self.root/'another')))
        self.assertEqual(first, (self.root/'another/dos.csv').read_bytes())
        # For a real probe, the first moment at 3 uniformly spaced k points
        # cancels image terms. The remaining two intracell hoppings give a
        # closed-form probe-by-probe variance (the same probes must span all k).
        signs = 2*np.random.default_rng(17).integers(0, 2, (8, 3), dtype=np.int8)-1
        diagonal = .4-.3+.1-3*report['kpm']['center_ev']
        direct = 2*(diagonal+2*(-.21*signs[:, 0]*signs[:, 1]-.14*signs[:, 1]*signs[:, 2]))/report['kpm']['scale_ev']
        rows = np.genfromtxt(out/'moments.csv', names=True, delimiter=',')
        self.assertAlmostEqual(rows['trace_moment'][1], direct.mean())
        self.assertAlmostEqual(rows['trace_standard_error'][1], direct.std(ddof=1)/np.sqrt(8))
        (out/'moments.csv').write_text('modified')
        with self.assertRaisesRegex(ValueError, 'missing or changed'):
            run(s, resume=True)

    def test_scaling_and_method_guards(self):
        from wren.dos.runner import run
        for options in ({'kpm-trace': 'exact', 'kpm-seed': 3}, {'sigma': .1}, {'kpm-vectors': 1},
                        {'kpm-moments': 1}, {'site-projections': True}, {'kpm-bounds-ev': [1, -1]},
                        {'solver': 'diagonalization', 'kpm-seed': 1}):
            with self.assertRaises(ValueError):
                self.settings(**options)
        with self.assertRaisesRegex(ValueError, 'enclose'):
            run(self.settings(**{'kpm-bounds-ev': [-.1, .1]}))
        self.path.write_text(json.dumps(dict(self.network, site_energies_ev=[-2.], pairs=[])))
        with self.assertRaisesRegex(ValueError, 'flat spectrum'):
            run(self.settings(output='flat'))
        result = run(self.settings(output='finite', **{'kpm-bounds-ev': [-3, -1]}))
        self.assertAlmostEqual(result['finite_grid_integral'], 2.)
        self.path.write_text(json.dumps({'schema': 'wren-dos-tb-2', 'site_energies_ev': [0., 0.],
                    'matrix_element_convention': 'raw_signed_H_and_S_in_same_orbital_gauge',
                    'hamiltonian_pairs_ev': [], 'overlap_pairs': [[0, 1, [0, 0, 0], .1]]}))
        with self.assertRaisesRegex(ValueError, 'orthonormal'):
            run(self.settings(output='overlap'))
