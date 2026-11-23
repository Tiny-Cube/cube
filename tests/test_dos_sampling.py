import importlib.util
import unittest


@unittest.skipUnless(importlib.util.find_spec('numpy'), 'DOS NumPy dependency required')
class DosSamplingTests(unittest.TestCase):
    def test_preflight_is_json_serializable_and_matches_analysis(self):
        import json
        import numpy as np
        from wren.dos.sampling import block_curve_statistics, validate_sampling

        # Unsigned frame indices must not wrap into apparently positive strides.
        with self.assertRaisesRegex(ValueError, 'strictly increasing'):
            validate_sampling(np.array([4, 3, 2, 1], dtype=np.uint32), 2)
        frames = np.arange(10, 22, 3, dtype=np.int32)
        preflight = validate_sampling(frames, np.int64(2), times_ps=[0., .3, .6, .9])
        self.assertEqual(json.loads(json.dumps(preflight)), preflight)
        result = block_curve_statistics([[0]]*4, frames, 2, times_ps=[0., .3, .6, .9])
        self.assertEqual(result['metadata'], preflight)

    def test_exact_block_errors_and_excluded_tail(self):
        import numpy as np
        from wren.dos.sampling import block_curve_statistics

        # At E1 the three complete block means are 1, 4, 7, with sample
        # variance 9. E2 is constant. E3 has means 6, 0, 0 (variance 12).
        curves = [[0, 2, 5], [2, 2, 7], [3, 2, 0], [5, 2, 0],
                  [6, 2, 0], [8, 2, 0], [100, 2, 100]]
        result = block_curve_statistics(curves, range(20, 34, 2), 2, times_ps=np.arange(7)*.25+4.20)
        np.testing.assert_array_equal(result['block_means'], [[1, 2, 6], [4, 2, 0], [7, 2, 0]])
        np.testing.assert_array_equal(result['retained_mean'], [4, 2, 2])
        np.testing.assert_allclose(result['standard_error'], [np.sqrt(3), 0, 2])
        metadata = result['metadata']
        self.assertEqual((metadata['full_blocks'], metadata['used_frames'], metadata['discarded_tail_frames']),
                         (3, 6, 1))
        self.assertEqual((metadata['first_frame'], metadata['last_used_frame'], metadata['last_input_frame']),
                         (20, 30, 32))
        self.assertEqual(metadata['degrees_of_freedom'], 2)
        self.assertAlmostEqual(metadata['sample_spacing_ps'], .25)
        self.assertFalse(metadata['block_independence_verified'])
        self.assertFalse(metadata['material_convergence_verified'])

    def test_block_size_one_and_correlated_adjacent_frames(self):
        from wren.dos.sampling import block_curve_statistics

        curves = [[0], [0], [2], [2]]
        plain = block_curve_statistics(curves, range(4), 1)
        paired = block_curve_statistics(curves, range(4), 2)
        self.assertAlmostEqual(plain['standard_error'][0], 1/3**.5)
        self.assertAlmostEqual(paired['standard_error'][0], 1)
        self.assertEqual(plain['retained_mean'][0], paired['retained_mean'][0])
        self.assertEqual(paired['metadata']['discarded_tail_frames'], 0)
        self.assertIsNone(paired['metadata']['sample_spacing_ps'])

    def test_reject_unusable_blocks_irregular_frames_and_nonfinite_data(self):
        from wren.dos.sampling import block_curve_statistics

        data = [[1], [2], [3], [4]]
        for size in (True, 0, -1, 1.5, 3):
            with self.subTest(size=size), self.assertRaises(ValueError):
                block_curve_statistics(data, range(4), size)
        for indices in ([0, 1, 3, 4], [3, 2, 1, 0], [0, 0, 0, 0],
                        [0, 1, 2], [True, 2, 3, 4], [0., 1., 2., 3.]):
            with self.subTest(indices=indices), self.assertRaises(ValueError):
                block_curve_statistics(data, indices, 2)
        for times in ([0, .5, 1, 1.7], [0, 1, 1, 2], [0, 1, None, 3], [0, 1, 2]):
            with self.subTest(times=times), self.assertRaises(ValueError):
                block_curve_statistics(data, range(4), 2, times_ps=times)
        for curves in ([1, 2, 3, 4], [[float('nan')]]*4, [[1j]]*4, [[]]*4):
            with self.subTest(curves=curves), self.assertRaises(ValueError):
                block_curve_statistics(curves, range(4), 2)
