import importlib.util
import math
from pathlib import Path
import sys
import unittest
from universal_kmc import Engine

spec = importlib.util.spec_from_file_location('generic_models', Path(__file__).resolve().parents[1]/'examples/generic.py')
generic = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = generic
spec.loader.exec_module(generic)


class GenericPhysicsTests(unittest.TestCase):
    def test_ising_detailed_balance(self):
        def energy(spins):
            return -sum(spins[i]*spins[(i+1)%len(spins)] for i in range(len(spins)))
        # Enumerate every state and transition of a three-spin ring.
        for bits in range(8):
            initial = [1 if bits & (1 << i) else -1 for i in range(3)]
            for event in range(3):
                model = generic.ising_model(3)
                model.data.spins = initial.copy()
                before = energy(model.data.spins)
                forward = model.propensities()[event]
                model.fire(event)
                delta = energy(model.data.spins)-before
                reverse = model.propensities()[event]
                self.assertAlmostEqual(forward/reverse, math.exp(-model.data.beta*delta), places=12)

    def test_queue_stationary_time_average(self):
        model = generic.queue_model()
        engine = Engine(model, seed=982341)
        weighted, duration = 0.0, 0.0
        # M/M/1 with arrival=2, service=3 has stationary E[queue length]=2.
        for step in range(100000):
            before = engine.time
            length = len(model.data.jobs)
            engine.step()
            if step >= 10000:
                weighted += length*(engine.time-before)
                duration += engine.time-before
        self.assertAlmostEqual(weighted/duration, 2.0, delta=.18)
        self.assertEqual(model.data.completed + len(model.data.jobs), model.data.next_id)


if __name__ == '__main__':
    unittest.main()
