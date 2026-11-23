"""CUDA direct SSA for independent carriers on a static directed network."""
from pathlib import Path
import random

import numpy as np

from .simulation import IndependentRNG, trajectory_seed

KERNEL = Path(__file__).with_name('walk.cu')
MAX_BATCH_TRAJECTORIES = 16384


def device_info():
    try:
        import cupy as cp
    except ImportError as exc:
        raise RuntimeError('KMC backend=gpu requires CuPy; install uv sync --locked --extra transport-cuda13') from exc
    if cp.cuda.runtime.is_hip or cp.cuda.runtime.getDeviceCount() < 1:
        raise RuntimeError('KMC backend=gpu requires an available NVIDIA CUDA device')
    device = cp.cuda.Device()
    name = cp.cuda.runtime.getDeviceProperties(device.id)['name']
    return {'device': name.decode() if isinstance(name, bytes) else str(name),
            'device_id': device.id, 'compute_capability': device.compute_capability,
            'cupy': cp.__version__, 'cuda_runtime': cp.cuda.runtime.runtimeGetVersion(),
            'cuda_driver': cp.cuda.runtime.driverGetVersion()}


def _graph(sites, edges):
    if max(sites, len(edges)) >= 2**31:
        raise ValueError('GPU KMC CSR indices exceed int32')
    outgoing = [[] for _ in range(sites)]
    for index, edge in enumerate(edges):
        if edge.rate > 0:
            outgoing[edge.source].append((index, edge))
    row, targets, ids, cdf, delta = [0], [], [], [], []
    for channels in outgoing:
        total = 0.
        for index, edge in channels:
            total += edge.rate  # Match the scalar addition order of wren.Engine.
            cdf.append(total)
            targets.append(edge.target)
            ids.append(index)
            delta.append(edge.displacement)
        if not np.isfinite(total):
            raise ValueError('total propensity overflow')
        row.append(len(ids))
    return (np.asarray(row, np.int32), np.asarray(targets, np.int32), np.asarray(ids, np.int32),
            np.asarray(cdf), np.array(delta, dtype=float).reshape(-1, 3))


def _seeds(first, width, equilibrium_cdf, settings, carrier):
    initial = np.empty(width, np.int32)
    words = np.empty((624, width), np.uint32)
    for local in range(width):
        trajectory = first + local
        draw = IndependentRNG(trajectory_seed(settings['seed'], carrier, trajectory, 'initial')).uniform()
        initial[local] = min(np.searchsorted(equilibrium_cdf, draw, side='right'), len(equilibrium_cdf)-1)
        seed = trajectory_seed(settings['seed'], carrier, trajectory, 'events')
        version, state, _ = random.Random(seed).getstate()
        if version != 3 or len(state) != 625 or state[-1] != 624:
            raise RuntimeError('GPU KMC requires the Python MT19937 624-word initial state')
        words[:, local] = state[:-1]
    return initial, words


def _trace(first, initial, edges, seed, times, events):
    site = int(initial)
    xyz = (0., 0., 0.)
    step = 0
    records = [{'step': 0, 'time_s': 0., 'site': site, 'event': None, 'displacement_m': list(xyz)}]
    for t, index in zip(times, events, strict=True):
        if index == -2:  # Last reserved entry unused after an event limit or exact final hop.
            break
        if index >= 0:
            step += 1
            edge = edges[index]
            site = edge.target
            xyz = tuple(a+b for a, b in zip(xyz, edge.displacement))
        records.append({'step': step, 'time_s': float(t), 'site': site,
                        'event': None if index == -1 else int(index), 'displacement_m': list(xyz)})
    return {'trajectory': first, 'initial_site': int(initial), 'event_seed': seed, 'records': records}


def _walk(kernel, graph, grid, max_events, initial, words, offsets=None):
    import cupy as cp
    n = len(initial)
    track = offsets is not None
    capacity = int(offsets[-1]) if track else 1
    if offsets is None:
        offsets = np.zeros(1, np.int64)
    samples = cp.empty(1 if track else (n, len(grid), 3), dtype=cp.float64)
    steps = cp.empty(n, dtype=cp.int64)
    limited = cp.empty(n, dtype=cp.int32)
    errors = cp.zeros(n, dtype=cp.int32)
    times = cp.empty(capacity, dtype=cp.float64)
    events = cp.full(capacity, -2, dtype=cp.int32)
    args = (np.int32(n), *graph, cp.asarray(initial), cp.asarray(np.ascontiguousarray(words)),
            np.int32(len(grid)), grid, np.int64(max_events), samples, steps,
            limited, errors, np.int32(track), cp.asarray(offsets), times, events)
    kernel(((n+127)//128,), (128,), args)
    error = int(errors.max().item())
    if error:
        raise RuntimeError('GPU KMC '+('waiting time outside floating-point range' if error == 1
                                      else 'trace replay exceeded first-pass event count'))
    return (None if track else cp.asnumpy(samples), cp.asnumpy(steps), cp.asnumpy(limited),
            cp.asnumpy(times) if track else None, cp.asnumpy(events) if track else None)


def simulate(sites, edges, pi, settings, carrier):
    import cupy as cp
    if settings['max-events'] > np.iinfo(np.int64).max:
        raise ValueError('GPU max-events exceeds int64')
    kernel = cp.RawKernel(KERNEL.read_text(), 'kmc_walk', options=('--std=c++11', '--fmad=false'))
    graph = tuple(cp.asarray(a) for a in _graph(sites, edges))
    grid = np.linspace(0., settings['duration'], settings['samples'])
    gpu_grid = cp.asarray(grid)
    total = settings['trajectories']
    positions = np.empty((total, len(grid), 3))
    counts, capped, traces = [], [], []
    equilibrium = np.cumsum(pi)
    batch = min(MAX_BATCH_TRAJECTORIES, max(1, 1_000_000//len(grid)))
    for first in range(0, total, batch):
        width = min(batch, total-first)
        initial, words = _seeds(first, width, equilibrium, settings, carrier)
        samples, steps, limited, _, _ = _walk(kernel, graph, gpu_grid, settings['max-events'], initial, words)
        positions[first:first+width] = samples
        counts.extend(steps.tolist())
        capped.extend((first+np.flatnonzero(limited)).tolist())
        tracked = min(width, max(0, settings['trace-trajectories']-first))
        if tracked:
            offsets = np.concatenate(([0], np.cumsum(steps[:tracked]+1, dtype=np.int64)))
            _, replay_steps, replay_limits, times, events = _walk(
                kernel, graph, gpu_grid, settings['max-events'], initial[:tracked], words[:, :tracked], offsets)
            if not (np.array_equal(steps[:tracked], replay_steps) and np.array_equal(limited[:tracked], replay_limits)):
                raise RuntimeError('GPU KMC trace replay disagrees with first pass')
            for local in range(tracked):
                part = slice(offsets[local], offsets[local+1])
                trajectory = first+local
                seed = trajectory_seed(settings['seed'], carrier, trajectory, 'events')
                traces.append(_trace(trajectory, initial[local], edges, seed, times[part], events[part]))
    return grid, positions, traces, counts, capped
