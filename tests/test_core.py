import math
import unittest
from universal_kmc import Engine, Event, EventModel, Graph, Edge, PortableRNG, Reaction, ReactionNetwork
from universal_kmc.core import MAX_COUNT
from universal_kmc.io import loads
from universal_kmc.rates import KB, arrhenius, marcus, miller_abrahams


class Fixed:
    def __init__(self, *values):
        self.values = iter(values)
    def uniform(self):
        return next(self.values)


class CoreTests(unittest.TestCase):
    def test_minimal_model_needs_only_two_methods(self):
        class Custom:
            def __init__(self):
                self.state = {'history': [], 'energy': 1.25}
                self.position = 'domain-specific label'
                self.data = self.state
            def propensities(self):
                return [2.0] if len(self.data['history']) < 2 else []
            def fire(self, event):
                self.data['history'].append('changed')
                self.data['energy'] *= 0.5
        model = Custom()
        engine = Engine(model)
        rows = list(engine.run(10))
        self.assertEqual(len(rows), 3)
        self.assertEqual(engine.status, 'absorbed')
        self.assertEqual(model.data['energy'], 0.3125)
        self.assertEqual(rows[-1].state, ())

    def test_arbitrary_state_and_dynamic_events(self):
        data = {'pending': ['alpha', 'beta', 'gamma'], 'done': []}
        def events(data):
            return [Event(lambda s: 1.0,
                          lambda s, name=name: (s['pending'].remove(name), s['done'].append(name)))
                    for name in data['pending']]
        model = EventModel(data, events)
        engine = Engine(model)
        list(engine.run(10))
        self.assertEqual(engine.steps, 3)
        self.assertEqual(engine.status, 'absorbed')
        self.assertEqual(sorted(data['done']), ['alpha', 'beta', 'gamma'])
        self.assertEqual(data['pending'], [])

    def test_disabled_callback_not_evaluated(self):
        model = EventModel({}, [Event(lambda s: math.nan, lambda s: None, enabled=lambda s: False)])
        self.assertEqual(model.propensities(), (0.0,))
        with self.assertRaises(ValueError): model.fire(0)

    def test_exact_clock_selection_and_terminal_observation(self):
        model = Graph(1, [Edge(0, 0, 2), Edge(0, 0, 3)])
        engine = Engine(model, rng=Fixed(math.exp(-1), 0.8, math.exp(-1)))
        rows = list(engine.run(10, until=0.3))
        self.assertEqual([r.event for r in rows], [-1, 1, -1])
        self.assertEqual([r.step for r in rows], [0, 1, 1])
        self.assertAlmostEqual(rows[1].time, 0.2)
        self.assertEqual(rows[-1].time, 0.3)
        self.assertEqual(engine.status, 'time_limit')

    def test_absorbed_horizon_zero_steps_and_boundary(self):
        engine = Engine(Graph(1, []))
        self.assertEqual([r.time for r in engine.run(5, 3)], [0, 3])
        self.assertEqual(engine.status, 'absorbed')
        engine = Engine(Graph(1, []))
        self.assertEqual(len(list(engine.run(0, 10))), 1)
        self.assertEqual(engine.time, 0)
        engine = Engine(Graph(1, [Edge(0, 0, 1)]), rng=Fixed(math.exp(-1), 0.5))
        self.assertIsNotNone(engine.step(1.0))  # event exactly at horizon is included
        self.assertIsNone(engine.step(1.0))

    def test_invalid_rates_and_rng(self):
        class Model:
            def propensities(self): return [math.inf]
            def fire(self, event): raise AssertionError('should not fire')
        with self.assertRaises(ValueError): Engine(Model()).step()
        for rate in [-1, math.nan, math.inf]:
            with self.assertRaises(ValueError): Graph(1, [Edge(0, 0, rate)])
        for u in [0, 1, math.nan, math.inf, -1]:
            engine = Engine(Graph(1, [Edge(0, 0, 1)]), rng=Fixed(u))
            with self.assertRaises(ValueError): engine.step()
        with self.assertRaises(ValueError): Engine(Graph(1, [Edge(0, 0, 1e308), Edge(0, 0, 1e308)])).step()
        engine = Engine(Graph(1, [Edge(0, 0, 1)]), rng=Fixed(0.5))
        engine.time = 1e300
        with self.assertRaises(ValueError): engine.step()
        with self.assertRaises(ValueError): PortableRNG(-1)
        with self.assertRaises(ValueError): list(Engine(Graph(1, [])).run(-1))
        with self.assertRaises(ValueError): list(Engine(Graph(1, [])).run(1, math.nan))

    def test_mass_action_and_conservation(self):
        model = ReactionNetwork([10, 0], [Reaction(0.1, (2, 0), (0, 1)), Reaction(0.3, (0, 1), (2, 0))])
        self.assertAlmostEqual(model.propensities()[0], 4.5)
        for row in Engine(model, 99).run(2000):
            self.assertEqual(row.state[0] + 2 * row.state[1], 10)
            self.assertTrue(all(n >= 0 for n in row.state))
        model = ReactionNetwork([MAX_COUNT], [Reaction(1, (0,), (1,))])
        with self.assertRaises(ValueError): model.fire(0)
        self.assertEqual(model.state(), (MAX_COUNT,))
        model = ReactionNetwork([2], [Reaction(1e308, (3,), (0,))])
        self.assertEqual(model.propensities(), [0.0])

    def test_rate_laws_and_detailed_balance(self):
        self.assertAlmostEqual(marcus(.01, .2, .05) / marcus(.01, .2, -.05), math.exp(-.05/(KB*300)), places=13)
        self.assertEqual(marcus(0, .2, .1), 0)
        self.assertAlmostEqual(arrhenius(1e12, 0), 1e12)
        self.assertAlmostEqual(miller_abrahams(1e12, 1, 2, -.1), 1e12/math.e, places=3)
        for args in [(.01, 0, 0, 300), (.01, .2, 0, 0), (math.nan, .2, 0, 300)]:
            with self.assertRaises(ValueError): marcus(*args)

    def test_parser_validation(self):
        for text in ['bad', 'KMC1 graph 1 1 0 0 1 2 0 0 0',
                     'KMC1 graph 1 0 0 trailing', 'KMC1 graph 1 1 0 0 0 nan 0 0 0',
                     'KMC1 graph 1 1 0 0 0 -1 0 0 0', 'KMC1 reaction 0 0',
                     'KMC1 reaction 1 1 3 1 2', 'KMC1 graph +1 0 0']:
            with self.assertRaises(ValueError, msg=text): loads(text)
        self.assertEqual(loads('KMC1 graph # comment\n 1 0 0').state(), (0,))


class StatisticalTests(unittest.TestCase):
    def test_exponential_waiting_and_branch_probabilities(self):
        model = Graph(1, [Edge(0, 0, 2), Edge(0, 0, 8)])
        engine = Engine(model, 8273)
        n = 50000
        first, total, squares = 0, 0.0, 0.0
        for _ in range(n):
            before = engine.time
            row = engine.step()
            dt = row.time - before
            total += dt
            squares += dt * dt
            first += row.event == 0
        self.assertAlmostEqual(total/n, .1, delta=.002)
        self.assertAlmostEqual(squares/n - (total/n)**2, .01, delta=.0007)
        self.assertAlmostEqual(first/n, .2, delta=.008)

    def test_two_state_time_weighted_occupation(self):
        engine = Engine(Graph(2, [Edge(0, 1, 2), Edge(1, 0, 1)]), 9273)
        occupied, prev_time, prev_state = 0.0, 0.0, 0
        for row in engine.run(50000):
            occupied += (row.time - prev_time) * (prev_state == 1)
            prev_time, prev_state = row.time, row.state[0]
        self.assertAlmostEqual(occupied/engine.time, 2/3, delta=.012)

    def test_fixed_time_random_walk_and_birth_death(self):
        # One continuous RNG stream avoids treating adjacent seeds as independent streams.
        rng = PortableRNG(123456)
        positions = []
        populations = []
        for _ in range(2500):
            model = Graph(1, [Edge(0, 0, 1, (1, 0, 0)), Edge(0, 0, 1, (-1, 0, 0))])
            engine = Engine(model, rng=rng)
            list(engine.run(10000, 5))
            positions.append(model.position()[0])
            model = ReactionNetwork([0], [Reaction(4, (0,), (1,)), Reaction(.5, (1,), (0,))])
            engine = Engine(model, rng=rng)
            list(engine.run(10000, 3))
            populations.append(model.state()[0])
        mean = sum(positions)/len(positions)
        variance = sum((x-mean)**2 for x in positions)/len(positions)
        self.assertAlmostEqual(mean, 0, delta=.22)
        self.assertAlmostEqual(variance, 10, delta=.8)
        self.assertAlmostEqual(sum(populations)/len(populations), 8*(1-math.exp(-1.5)), delta=.18)


if __name__ == '__main__':
    unittest.main()
