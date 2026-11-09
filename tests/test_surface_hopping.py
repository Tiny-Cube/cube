"""Physics invariants and independent representation checks for real FSSH."""
from contextlib import redirect_stdout
import io
import json
from pathlib import Path
import tempfile
import unittest

try:
    import numpy as np
except ImportError:
    np = None

if np is not None:
    from universal_kmc.surface_hopping.engine import (
        FSSH, AdiabaticFrame, adiabatic_frame, electronic_step, rescale_momentum, energy_decoherence,
    )
    from universal_kmc.surface_hopping.models import TullyModel, LinearVibronicModel
from universal_kmc.surface_hopping.control import load_control
from universal_kmc.surface_hopping.runner import run_control


@unittest.skipIf(np is None, 'install the dynamics extra for FSSH physics tests')
class SurfaceHoppingTests(unittest.TestCase):
    def test_all_tully_derivatives_and_adiabatic_gradients(self):
        for number in (1, 2, 3):
            model = TullyModel(number)
            for x in (-2., -.5, 0., .5, 2.):
                with self.subTest(model=number, x=x):
                    epsilon = 2e-7
                    h, dh = model.evaluate([x])
                    fd = (model.evaluate([x+epsilon])[0]-model.evaluate([x-epsilon])[0])/(2*epsilon)
                    np.testing.assert_allclose(dh[..., 0], fd, atol=2e-8, rtol=2e-6)
                    frame = adiabatic_frame(model, np.array([x]))
                    fd_e = (np.linalg.eigvalsh(model.evaluate([x+epsilon])[0])-np.linalg.eigvalsh(model.evaluate([x-epsilon])[0]))/(2*epsilon)
                    np.testing.assert_allclose(frame.gradients[:, 0], fd_e, atol=2e-8)
                    np.testing.assert_allclose(frame.nac, -frame.nac.swapaxes(0, 1), atol=1e-14)
        # Original parameters at the crossing and asymptotic limits.
        np.testing.assert_allclose(TullyModel(1).evaluate([0])[0], [[0,.005],[.005,0]])
        np.testing.assert_allclose(TullyModel(2).evaluate([0])[0], [[0,.015],[.015,-.05]])
        np.testing.assert_allclose(TullyModel(3).evaluate([0])[0], [[.0006,.1],[.1,-.0006]])

    def test_phase_alignment_and_degenerate_rejection(self):
        model = TullyModel()
        frame = adiabatic_frame(model, np.array([.1]))
        previous = frame.vectors*np.array([-1, 1])
        aligned = adiabatic_frame(model, np.array([.1]), previous)
        np.testing.assert_allclose(aligned.vectors, previous)
        np.testing.assert_allclose(aligned.nac, frame.nac*np.array([-1,1])[:,None,None]*np.array([-1,1])[None,:,None])
        degenerate = LinearVibronicModel(np.zeros((2,2)), np.zeros((2,2,1)), [0])
        with self.assertRaisesRegex(ValueError, 'degenerate'):adiabatic_frame(degenerate, [0])

    def test_unitary_flux_matches_exact_rotation_and_population_loss(self):
        frame = AdiabaticFrame(np.zeros(2), np.eye(2), np.zeros((2,1)), np.array([[[0],[.3]],[[-.3],[0]]]))
        c, probability = electronic_step(np.array([1,0], dtype=complex), frame, frame, [1], 0, 1., 100)
        np.testing.assert_allclose(c, [np.cos(.3), np.sin(.3)], atol=1e-13)
        self.assertAlmostEqual(probability[1], np.sin(.3)**2, delta=2e-7)
        self.assertAlmostEqual(probability[1], 1-abs(c[0])**2, delta=2e-7)
        self.assertAlmostEqual(float(np.vdot(c,c).real), 1, delta=1e-13)

    def test_electronic_adiabatic_vs_independent_diabatic_propagation(self):
        # Prescribed straight nuclear path; no stochastic hopping in this check.
        # Compare two different representations, independently constructing H(q).
        for number in (1, 2, 3):
            model = TullyModel(number)
            dt, velocity, q = .5, .01, -2.
            frame = adiabatic_frame(model, [q])
            c = np.array([1,0], dtype=complex)
            diabatic = frame.vectors @ c
            for _ in range(800):
                new = adiabatic_frame(model, [q+dt*velocity], frame.vectors)
                c, _ = electronic_step(c, frame, new, [velocity], 0, dt, 4)
                h, _ = model.evaluate([q+.5*dt*velocity])
                e, u = np.linalg.eigh(h)
                diabatic = u @ (np.exp(-1j*e*dt)*(u.T @ diabatic))
                q, frame = q+dt*velocity, new
            fidelity = abs(np.vdot(diabatic, frame.vectors @ c))**2
            self.assertGreater(fidelity, 1-1e-7, (number, fidelity))

    def test_mass_weighted_momentum_jump_and_frustration(self):
        p, m, d = np.array([3.,4.]), np.array([2.,5.]), np.array([.4,-.7])
        for delta in (-.7, .002, 0.):
            adjusted, accepted = rescale_momentum(p,m,d,delta)
            self.assertTrue(accepted)
            self.assertAlmostEqual(float(np.sum((adjusted**2-p**2)/(2*m))+delta), 0, delta=1e-14)
            np.testing.assert_allclose((adjusted-p)[0]*d[1], (adjusted-p)[1]*d[0], atol=1e-14)
            for scale in (-1e12, 1e-12):
                other, _ = rescale_momentum(p,m,scale*d,delta)
                np.testing.assert_allclose(adjusted,other,atol=1e-13)
        rejected, accepted = rescale_momentum(p,m,d,100.)
        self.assertFalse(accepted);np.testing.assert_array_equal(rejected,p)
        zero, accepted = rescale_momentum(p,m,[0,0],.2)
        self.assertFalse(accepted);np.testing.assert_array_equal(zero,p)

    def test_unequal_mass_multistate_lvc_and_verlet_convergence(self):
        provider = LinearVibronicModel(np.diag([0.,.2,.5]), np.zeros((3,3,2)), [.1,.2])
        def run(dt):
            engine = FSSH(provider,[2.,5.],dt=dt,electronic_substeps=2)
            state = engine.initialize([1.,-.7],[.3,.4],active=1)
            initial = engine.energy(state)
            rng = np.random.default_rng(22)
            drift=0.
            for _ in range(round(10/dt)):
                state,event=engine.step(state,rng)
                drift=max(drift,abs(engine.energy(state)-initial))
                self.assertEqual(event['event'],'none');self.assertEqual(state.active,1)
            return drift
        coarse, fine = run(.1), run(.05)
        self.assertLess(fine, coarse*.26)

    def test_edc_decays_inactive_amplitude_and_preserves_norm_phase(self):
        c = np.array([.6j,.8],dtype=complex)
        new = energy_decoherence(c,np.array([0.,.2]),0,.5,2.,.1)
        self.assertAlmostEqual(abs(new[1]),.8*np.exp(-2*.2/(1+.1/.5)))
        self.assertAlmostEqual(float(np.vdot(new,new).real),1.)
        self.assertAlmostEqual(new[0].real,0.);self.assertGreater(new[0].imag,0.)
        np.testing.assert_allclose(energy_decoherence(c,np.array([0.,.2]),0,0.,2.),c)

    def test_selected_frustrated_hop_keeps_momentum_and_reports_target(self):
        class Draw:
            def __init__(self,value):self.value=value
            def random(self):return self.value
        engine=FSSH(TullyModel(),[2000.],dt=.1,electronic_substeps=4)
        frame=adiabatic_frame(TullyModel(),[-.1])
        sign=np.sign(frame.nac[0,1,0])
        state=engine.initialize([-.1],[1.],amplitudes=[np.sqrt(.8),sign*np.sqrt(.2)])
        frustrated,event=engine.step(state,Draw(0.))
        unchanged,_=engine.step(state,Draw(.999999))
        self.assertEqual(event['event'],'frustrated');self.assertEqual(event['to'],1)
        self.assertEqual(frustrated.active,0)
        np.testing.assert_array_equal(frustrated.p,unchanged.p)
        np.testing.assert_array_equal(frustrated.c,unchanged.c)

    def test_overlarge_probabilities_are_rejected_without_renormalization(self):
        frame=AdiabaticFrame(np.zeros(2),np.eye(2),np.zeros((2,1)),np.array([[[0],[.3]],[[-.3],[0]]]))
        with self.assertRaisesRegex(RuntimeError,'exceeds one'):
            electronic_step([1,0],frame,frame,[1],0,5.,1)

    def test_invalid_providers_initial_states_and_controls_fail(self):
        class Bad:
            def evaluate(self,q):return np.array([[0.,1.],[0.,2.]]), np.zeros((2,2,1))
        with self.assertRaisesRegex(ValueError,'symmetric'):adiabatic_frame(Bad(),[0])
        with self.assertRaises(ValueError):FSSH(TullyModel(),[-1])
        engine=FSSH(TullyModel(),[2000])
        for q,p,a,c in [([0],[1],2,None),([0,1],[1],0,None),([0],[1],0,[0,1]),([0],[1],0,[1,1])]:
            with self.assertRaises(ValueError):engine.initialize(q,p,a,c)


SMALL='''[model]
type="tully-1"
[initial]
coordinates=[-0.5]
momenta=[20.0]
masses=[2000.0]
[dynamics]
trajectories=3
max-steps=30
seed=42
output-stride=10
'''


@unittest.skipIf(np is None, 'install the dynamics extra for FSSH execution tests')
class SurfaceControlTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name);self.path=self.root/'dynamics.toml';self.path.write_text(SMALL)

    def test_audited_resume_reproducibility_and_tamper_detection(self):
        with redirect_stdout(io.StringIO()):
            result=run_control(self.path,output=self.root/'a')
            second=run_control(self.path,output=self.root/'b')
            reused=run_control(self.path,output=self.root/'a',resume=True)
        self.assertEqual(result['final_states'],second['final_states']);self.assertTrue(reused['reused'])
        self.assertLess(result['max_norm_error'],1e-12)
        self.assertLess(result['max_energy_drift_hartree'],1e-7)
        self.assertEqual(sum(result['active_counts']),3)
        with self.assertRaisesRegex(ValueError,'existing'):run_control(self.path,output=self.root/'a')
        (self.root/'a/trajectories.jsonl').write_text('changed')
        with self.assertRaisesRegex(ValueError,'changed'):run_control(self.path,output=self.root/'a',resume=True)

    def test_strict_parameters_and_no_fallback_to_qc(self):
        for text in (SMALL.replace('trajectories=3','trajectories=true'), SMALL.replace('max-steps=30','max-steps=0'),
                     SMALL.replace('seed=42','seed=-1'), SMALL.replace('tully-1','gpu-pyscf'),
                     SMALL.replace('[model]','[model]\nunknown=1'), SMALL+'[scattering]\nleft=1\nright=2\n'):
            self.path.write_text(text)
            with self.assertRaises(ValueError):load_control(self.path)

    def test_workflow_gating_and_required_dynamics_before_science(self):
        from universal_kmc.workflow.control import load_input
        from universal_kmc.workflow.runner import run_workflow
        path=self.root/'input'
        path.write_text('$surface-hopping\ndo=1\n$end\n')
        with self.assertRaisesRegex(ValueError,'requires.*surface-hopping-dynamics'):load_input(path).execution_plan()
        path.write_text('$surface-hopping\ndo=0\n$end\n$surface-hopping-dynamics\ninput=absent.toml\n$end\n')
        control=load_input(path);self.assertEqual(control.plan()['tasks'],[])
        self.assertEqual(control.plan(['surface-hopping-dynamics'])['tasks'],[])
        path.write_text('$surface-hopping\ndo=1\n$end\n$surface-hopping-dynamics\ninput=dynamics.toml\noutput=run\n$end\n')
        with redirect_stdout(io.StringIO()):self.assertEqual(run_workflow(load_input(path)),0)
        report=json.loads((self.root/'cube-results/workflow.json').read_text())
        self.assertEqual(report['calculations']['surface-hopping']['scope'],'model_hamiltonian')
        self.assertEqual(report['status'],'completed')

    def test_unfinished_scattering_is_partial_and_not_resumable(self):
        self.path.write_text(SMALL+'[scattering]\nleft=-10\nright=8\n')
        with redirect_stdout(io.StringIO()):result=run_control(self.path)
        self.assertEqual(result['status'],'partial')
        self.assertEqual(result['exits']['max_steps'],3)
        with self.assertRaisesRegex(ValueError,'completed'):run_control(self.path,resume=True)
        from universal_kmc.workflow.control import load_input
        from universal_kmc.workflow.runner import run_workflow
        path=self.root/'input'
        path.write_text('$surface-hopping\ndo=1\n$end\n$surface-hopping-dynamics\ninput=dynamics.toml\noutput=other\n$end\n')
        with redirect_stdout(io.StringIO()):self.assertEqual(run_workflow(load_input(path)),2)
        report=json.loads((self.root/'cube-results/workflow.json').read_text())
        self.assertEqual(report['calculations']['surface-hopping']['status'],'partial')


if __name__=='__main__':unittest.main()
