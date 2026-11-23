"""Material CUDA SSA against the scalar Python engine with identical streams."""
from contextlib import redirect_stdout
import ctypes
import importlib.util
import io
import json
import math
import os
from pathlib import Path
import random
import shutil
import subprocess
import sys
import tempfile
import types
import unittest
from unittest.mock import patch


HAVE_NUMPY = importlib.util.find_spec('numpy') is not None
ROOT = Path(__file__).resolve().parents[1]


def settings(**overrides):
    return {'duration': 130., 'samples': 47, 'trajectories': 7, 'trace-trajectories': 3,
            'seed': 123459, 'max-events': 10000, **overrides}


def check_walk(case, cfg, edges, pi):
    import numpy as np
    from wren.transport import gpu, simulation
    with redirect_stdout(io.StringIO()):
        cpu = simulation.simulate(len(pi), edges, pi, cfg, 'hole')
        candidate = gpu.simulate(len(pi), edges, pi, cfg, 'hole')
    np.testing.assert_array_equal(candidate[0], cpu[0])
    np.testing.assert_allclose(candidate[1], cpu[1], atol=0, rtol=1e-14, equal_nan=True)
    case.assertEqual(candidate[3:], cpu[3:])
    case.assertEqual(len(candidate[2]), len(cpu[2]))
    for a, b in zip(candidate[2], cpu[2], strict=True):
        for field in ('initial_site', 'event_seed', 'trajectory'):
            case.assertEqual(a[field], b[field])
        case.assertEqual(len(a['records']), len(b['records']))
        for x, y in zip(a['records'], b['records'], strict=True):
            for field in ('step', 'event', 'site', 'displacement_m'):
                case.assertEqual(x[field], y[field])
            case.assertTrue(math.isclose(x['time_s'], y['time_s'], rel_tol=1e-13, abs_tol=1e-20))
    return candidate


def mixed_channels():
    from wren.models import Edge
    return [Edge(0, 1, 2., (1., 2., 0.)), Edge(0, 0, 0., (99., 0., 0.)),
            Edge(0, 1, 2.5, (-1., 0., 1.)), Edge(1, 0, 3., (-1., -2., 0.)),
            Edge(1, 0, 1., (1., 0., -1.)), Edge(0, 0, 1., (3., 0., -2.))]


@unittest.skipUnless(HAVE_NUMPY, 'requires NumPy')
class DeviceSelectionTests(unittest.TestCase):
    def test_missing_cupy_never_runs_python_as_gpu(self):
        from wren.transport.config import normalize
        from wren.transport.runner import run_network
        data = json.loads((ROOT/'examples/charge-carrier/chain.json').read_text())
        with tempfile.TemporaryDirectory() as temp:
            p = Path(temp)
            cfg = normalize({'backend': 'gpu'}, p)
            with patch.dict(sys.modules, {'cupy': None}), patch('wren.transport.simulation.simulate') as cpu:
                with self.assertRaisesRegex(RuntimeError, 'transport-cuda13'):
                    run_network(data, cfg, p/'test', {})
                cpu.assert_not_called()
            self.assertFalse((p/'test/result.json').exists())


@unittest.skipUnless(HAVE_NUMPY and shutil.which('c++'), 'requires NumPy and C++ compiler')
class NativeKernelOracleTests(unittest.TestCase):
    """Compile the same kernel arithmetic as host C++; no CUDA hardware implied."""
    @classmethod
    def setUpClass(cls):
        import numpy as np
        from wren.transport.gpu import KERNEL
        temp = tempfile.TemporaryDirectory()
        cls.addClassCleanup(temp.cleanup)
        libpath = Path(temp.name)/'walk.so'
        prefix = '''#include <cmath>
using std::log; using std::nan; using std::isfinite;
#define __device__
#define __global__
struct D { int x; }; D threadIdx = {0}, blockDim = {1}, blockIdx = {0};
extern "C" void set_lane(int n) { threadIdx.x = n; }
'''
        suffix = '''
extern "C" void mt_draws(unsigned int* words, double* out, int n) {
    int cursor=0; for(int i=0;i<n;++i) out[i]=positive_uniform(words,1,cursor);
}'''
        subprocess.run(['c++', '-x', 'c++', '-', '-std=c++11', '-O2', '-ffp-contract=off',
                        '-shared', '-fPIC', '-o', str(libpath)], input=prefix+KERNEL.read_text()+suffix,
                       text=True, check=True, capture_output=True)
        cls.native = ctypes.CDLL(str(libpath))
        cls.native.kmc_walk.restype = None
        def launch(blocks, threads, args):
            converted = [ctypes.c_void_p(a.ctypes.data) if isinstance(a, np.ndarray)
                         else ctypes.c_longlong(a) if isinstance(a, np.int64)
                         else ctypes.c_int(a) for a in args]
            for i in range(int(args[0])):
                cls.native.set_lane(i)
                cls.native.kmc_walk(*converted)
        cls.array_api = types.SimpleNamespace(
            RawKernel=lambda *a, **kw: launch, asarray=np.array, asnumpy=np.asarray,
            empty=np.empty, full=np.full, zeros=np.zeros,
            int32=np.int32, int64=np.int64, float64=np.float64)

    def test_seeded_uniforms_exactly_match_three_twist_cycles(self):
        import numpy as np
        for seed in (0, 9876543210123456789):
            rng = random.Random(seed)
            words = np.asarray(rng.getstate()[1][:-1], dtype=np.uint32)
            values = np.empty(1100)
            self.native.mt_draws(ctypes.c_void_p(words.ctypes.data), ctypes.c_void_p(values.ctypes.data), len(values))
            np.testing.assert_array_equal(values, [rng.random() for _ in values])

    def test_mixed_edges_absorption_event_cap_and_trace_indices(self):
        import numpy as np
        from wren.models import Edge
        with patch.dict(sys.modules, {'cupy': self.array_api}):
            result = check_walk(self, settings(), mixed_channels(), [.3, .7])
            self.assertGreater(min(result[3]), 312)
            result = check_walk(self, settings(**{'max-events': 2}), mixed_channels(), [.3, .7])
            self.assertEqual(result[4], list(range(7)))
            self.assertTrue(np.isnan(result[1][:, -1]).all())
            result = check_walk(self, settings(), [Edge(0, 0, 0., (0., 0., 0.))], [1.])
            self.assertTrue(np.all(result[1] == 0))
            self.assertIsNone(result[2][0]['records'][-1]['event'])

    def test_exact_observation_and_end_times_and_multiple_batches(self):
        import numpy as np
        from wren.models import Edge
        from wren.transport.simulation import IndependentRNG, trajectory_seed
        edge = Edge(0, 0, 1., (1., 0., 0.))
        t = -math.log(IndependentRNG(trajectory_seed(123459, 'hole', 0, 'events')).uniform())
        with patch.dict(sys.modules, {'cupy': self.array_api}):
            r = check_walk(self, settings(duration=t, samples=5, trajectories=1), [edge], [1.])
            self.assertEqual((r[3], r[2][0]['records'][-1]['event']), ([1], 0))
            np.testing.assert_array_equal(r[1][0,:,0], [0., 0., 0., 0., 1.])
            r = check_walk(self, settings(duration=4*t, samples=5, trajectories=1), [edge], [1.])
            self.assertEqual(r[1][0,1,0], 1.)
            with patch('wren.transport.gpu.MAX_BATCH_TRAJECTORIES', 3):
                check_walk(self, settings(duration=.2, samples=5, trajectories=10,
                                          **{'trace-trajectories': 5}), mixed_channels(), [.3, .7])


@unittest.skipUnless(HAVE_NUMPY and os.environ.get('WREN_TEST_GPU_KMC') == '1', 'explicit CUDA hardware test')
class ActualGpuTests(unittest.TestCase):
    def test_device_walk_and_public_workflow(self):
        import numpy as np
        from wren.transport.gpu import device_info
        from wren.workflow.cli import main
        print('KMC device', json.dumps(device_info(), sort_keys=True))
        check_walk(self, settings(), mixed_channels(), [.3, .7])
        check_walk(self, settings(**{'max-events': 1}), mixed_channels(), [.3, .7])
        check_walk(self, settings(), [], [1.])
        with tempfile.TemporaryDirectory() as temp:
            p = Path(temp)
            (p/'network.json').write_bytes((ROOT/'examples/charge-carrier/chain.json').read_bytes())
            records = {}
            for backend in ('python', 'gpu'):
                cfg = p/('input-'+backend)
                cfg.write_text(f'$kmc\ndo=1\nbackend={backend}\ninput=network.json\noutput={backend}\n'
                               'carriers=hole\ntrajectories=256\nseed=43\nsamples=31\n'
                               'duration=2e-10\nrelative-tolerance=10\n$end\n')
                self.assertEqual(main([str(cfg)]), 0)
                self.assertEqual(main([str(cfg), '--resume']), 0)
                records[backend] = json.loads((p/backend/'result.json').read_text())
            self.assertIn('gpu_runtime', records['gpu']['input'])
            cpu, gpu = (records[b]['carriers']['hole'] for b in ('python', 'gpu'))
            self.assertEqual(cpu['events_per_trajectory'], gpu['events_per_trajectory'])
            with np.load(p/'gpu/hole_samples.npz') as a, np.load(p/'python/hole_samples.npz') as b:
                np.testing.assert_array_equal(a['time_s'], b['time_s'])
                np.testing.assert_allclose(a['displacement_m'], b['displacement_m'], atol=0, rtol=1e-14)
            np.testing.assert_allclose(cpu['estimate']['diffusion_tensor_m2_s'], gpu['estimate']['diffusion_tensor_m2_s'])
            print('KMC hardware comparison: events', gpu['total_events'], 'Dxx reference',
                  gpu['exact_reference']['diffusion_tensor_m2_s'][0][0], 'Dxx estimated',
                  gpu['estimate']['diffusion_tensor_m2_s'][0][0], 'single-run seconds',
                  {name: record['elapsed_seconds'] for name, record in records.items()})
