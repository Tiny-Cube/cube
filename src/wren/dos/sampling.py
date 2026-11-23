"""Fixed-size, nonoverlapping batch means for a regularly sampled MD DOS.

See Flyvbjerg & Petersen, J. Chem. Phys. 91, 461–466 (1989),
https://doi.org/10.1063/1.457480. Selecting a single block length does
not demonstrate the error plateau studied in their blocking analysis.
"""
import math
from numbers import Integral

import numpy as np


def validate_sampling(frame_indices, block_size, *, times_ps=None):
    """Validate frame order, optional times, and the number of complete blocks."""
    if not isinstance(block_size, Integral) or isinstance(block_size, bool) or block_size < 1:
        raise ValueError('block_size must be a positive integer number of frames')
    block_size = int(block_size)
    try:
        indices = list(frame_indices)
    except TypeError as exc:
        raise ValueError('frame indices must be a sequence of nonnegative integers') from exc
    if any(not isinstance(i, Integral) or isinstance(i, bool) or i < 0 for i in indices):
        raise ValueError('frame indices must be nonnegative integers')
    indices = [int(i) for i in indices]
    blocks, tail = divmod(len(indices), block_size)
    if blocks < 2:
        raise ValueError('block_size must leave at least two complete blocks')
    stride = indices[1] - indices[0]
    if stride <= 0 or any(b - a != stride for a, b in zip(indices, indices[1:])):
        raise ValueError('frame indices must be strictly increasing and equally spaced')

    dt = None
    if times_ps is not None:
        times = np.asarray(times_ps)
        if times.shape != (len(indices),) or times.dtype.kind not in 'fiu' or not np.isfinite(times).all():
            raise ValueError('times_ps must contain one finite real time for each frame')
        intervals = np.diff(times.astype(float))
        dt = float(np.mean(intervals))
        if np.any(intervals <= 0) or not np.allclose(intervals, dt, rtol=1e-5, atol=1e-10):
            raise ValueError('times_ps must be strictly increasing and equally spaced')

    used = blocks * block_size
    return {
        'method': 'nonoverlapping_block_means',
        'error_scope': 'pointwise_standard_error_of_complete_block_mean',
        'block_size_frames': block_size, 'full_blocks': blocks, 'degrees_of_freedom': blocks - 1,
        'input_frames': len(indices), 'used_frames': used, 'discarded_tail_frames': tail,
        'first_frame': indices[0], 'last_used_frame': indices[used - 1],
        'last_input_frame': indices[-1], 'frame_stride': stride, 'sample_spacing_ps': dt,
        'block_independence_verified': False, 'material_convergence_verified': False,
    }


def block_curve_statistics(curves, frame_indices, block_size, *, times_ps=None):
    """Estimate pointwise standard errors of the mean over complete blocks.

    Rows are ordered trajectory frames and columns are points on one common
    energy grid. Each frame has equal statistical weight. For B complete
    blocks of the requested length, SE = std(block_means, ddof=1)/sqrt(B).
    The trailing incomplete block is excluded from both retained_mean and SE.
    These errors describe the sampled trajectory conditional on the chosen
    blocking; sufficiently decorrelated blocks require separate verification.
    """
    meta = validate_sampling(frame_indices, block_size, times_ps=times_ps)
    values = np.asarray(curves)
    if (values.ndim != 2 or not values.size or values.dtype.kind not in 'fiu'
            or not np.isfinite(values).all()):
        raise ValueError('curves must be a finite real frames-by-energy array')
    if values.shape[0] != meta['input_frames']:
        raise ValueError('each curve needs exactly one frame index')
    blocks = meta['full_blocks']
    means = values[:meta['used_frames']].reshape(blocks, block_size, values.shape[1]).mean(axis=1, dtype=float)
    return {
        'block_means': means,
        'retained_mean': means.mean(axis=0),
        'standard_error': means.std(axis=0, ddof=1) / math.sqrt(blocks),
        'metadata': meta,
    }
