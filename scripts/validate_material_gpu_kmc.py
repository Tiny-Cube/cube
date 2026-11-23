#!/usr/bin/env python3
"""Run both material KMC backends against the archived e_006 transfer graph."""
import argparse
from collections import deque
import csv
import json
import math
from pathlib import Path

import numpy as np

from wren.io import load
from wren.transport.network import build_rates
from wren.workflow.cli import main
from wren.workflow.dimers import file_hash


def archived_network(directory):
    audit = json.loads((directory/'input_audit.json').read_text())
    with (directory/'rates_and_geometry.csv').open() as stream:
        rows = list(csv.DictReader(stream))
    count = audit['molecules']
    box = np.array(audit['box_nm'])*10
    pairs, neighbors = [], [[] for _ in range(count)]
    for i, row in enumerate(rows):
        a, b = int(row['mol_i'])-1, int(row['mol_j'])-1
        image = [int(row['image_'+axis]) for axis in 'abc']
        delta = np.array([float(row['d'+axis+'_nm'])*10 for axis in 'xyz'])-box*image
        neighbors[a].append((b, delta))
        neighbors[b].append((a, -delta))
        pairs.append({'id': str(i+1), 'source': str(a+1), 'target': str(b+1), 'image': image,
                      'coupling_ev': {'hole': float(row['VH_meV'])/1000, 'electron': float(row['VL_meV'])/1000},
                      'reorganization_ev': audit['lambda_eV']})
    # Only position differences enter this graph. Recover a consistent set of
    # centers from its measured displacements, checking every closed loop.
    centers = {0: np.zeros(3)}
    queue = deque([0])
    while queue:
        i = queue.popleft()
        for j, delta in neighbors[i]:
            target = centers[i]+delta
            if j in centers:
                np.testing.assert_allclose(centers[j], target, rtol=0, atol=1e-10)
            else:
                centers[j] = target
                queue.append(j)
    if len(centers) != count:
        raise ValueError('incomplete archived graph')
    sites = [{'id': str(i+1), 'position_angstrom': centers[i].tolist(),
              'free_energy_ev': {'hole': 0., 'electron': 0.}} for i in range(count)]
    return {'schema': 'cube-charge-network-1', 'lattice_angstrom': np.diag(box).tolist(),
            'sites': sites, 'pairs': pairs,
            'provenance': {'system': directory.name, 'temperature_k': audit['temperature_K'],
                           'frame': 2000, 'uniform_site_energies_ev': 0.}}


def run(root, archived, backends, trajectories, duration):
    root.mkdir(parents=True, exist_ok=True)
    network = archived_network(archived)
    for carrier in ('hole', 'electron'):
        old = load(archived/f'{carrier}.kmc')
        _, edges, _, _, _ = build_rates(network, carrier, 300.)
        assert len(edges) == len(old.edges)
        for a, b in zip(edges, old.edges, strict=True):
            assert (a.source, a.target) == (b.source, b.target)
            np.testing.assert_allclose(a.displacement, b.displacement, rtol=0, atol=1e-21)
            if not math.isclose(a.rate, b.rate, rel_tol=1e-12):
                raise ValueError('recomputed Marcus rate differs from archived graph')
    path = root/'network.json'
    path.write_text(json.dumps(network, allow_nan=False))
    report = {'system': archived.name, 'sites': len(network['sites']), 'pairs': len(network['pairs']),
              'network_sha256': file_hash(path), 'duration_seconds': duration, 'trajectories_per_carrier': trajectories,
              'samples': 51, 'seed': 20261109, 'runs': {}}
    for backend in backends:
        control = root/(backend+'.input')
        control.write_text(f'$kmc\ndo=1\ninput=network.json\noutput=run-{backend}\nbackend={backend}\n'
                           f'carriers=hole,electron\ntemperature=300\ntrajectories={trajectories}\nduration={duration}\n'
                           'samples=51\ntrace-trajectories=3\nseed=20261109\n$end\n')
        code = main([str(control)])
        out = root/f'run-{backend}'
        result = json.loads((out/'result.json').read_text())
        if code not in (0, 2) or result['status'] not in ('completed', 'partial'):
            raise RuntimeError(f'{backend} material KMC failed')
        summary = {'elapsed_seconds': result['elapsed_seconds'], 'status': result['status'], 'carriers': {}}
        for carrier, data in result['carriers'].items():
            estimate = data['estimate']
            summary['carriers'][carrier] = {
                'total_events': data['total_events'], 'sampler': data['sampler'],
                'mobility_tensor_cm2_Vs': estimate['einstein_mobility_tensor_cm2_Vs'],
                'exact_mobility_tensor_cm2_Vs': data['exact_reference']['einstein_mobility_tensor_cm2_Vs'],
                'sampling_diagnostics': estimate['reason']}
        report['runs'][backend] = summary
    if set(backends) == {'python', 'gpu'}:
        comparisons = {}
        for carrier in ('hole', 'electron'):
            outputs = [root/f'run-{b}' for b in ('python', 'gpu')]
            results = [json.loads((p/'result.json').read_text())['carriers'][carrier] for p in outputs]
            assert results[0]['events_per_trajectory'] == results[1]['events_per_trajectory']
            with np.load(outputs[0]/f'{carrier}_samples.npz') as a, np.load(outputs[1]/f'{carrier}_samples.npz') as b:
                np.testing.assert_array_equal(a['time_s'], b['time_s'])
                np.testing.assert_array_equal(a['displacement_m'], b['displacement_m'])
            left, right = [json.loads((p/f'{carrier}_traces.json').read_text()) for p in outputs]
            max_dt, records = 0., 0
            for one, two in zip(left, right, strict=True):
                assert one['initial_site'] == two['initial_site']
                for x, y in zip(one['records'], two['records'], strict=True):
                    for key in ('step','event','site','displacement_m'):
                        assert x[key] == y[key]
                    max_dt = max(max_dt, abs(x['time_s']-y['time_s']))
                    records += 1
            if max_dt > duration * 1e-11:
                raise ValueError('event times differ beyond roundoff')
            comparisons[carrier] = {'exact_sampled_displacements': True, 'exact_per_trajectory_counts': True,
                                    'trace_records_checked': records, 'max_abs_event_time_difference_s': max_dt}
        report['comparisons'] = comparisons
    (root/'comparison.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    parser.add_argument('--cpu-only', action='store_true')
    parser.add_argument('--trajectories', type=int, default=256)
    parser.add_argument('--duration', type=float, default=2e-9)
    args = parser.parse_args()
    archived = Path(__file__).resolve().parents[1]/'runs/e_006_ether_12-17_0'
    run(args.output.resolve(), archived, ['python'] if args.cpu_only else ['python', 'gpu'], args.trajectories, args.duration)
