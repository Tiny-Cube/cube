#!/usr/bin/env python3
"""Geometric examples for the cube design discussion, not a dimer extractor.

Coordinates and distances are in nm. Needs NumPy and Matplotlib.
"""
import argparse
import json
import math
import os
from pathlib import Path

import numpy as np


def measures(a, b):
    ca, cb = a.mean(axis=0), b.mean(axis=0)
    _, _, va = np.linalg.svd(a-ca)
    _, _, vb = np.linalg.svd(b-cb)
    na, nb = va[-1], vb[-1]
    cosine = float(np.clip(abs(na @ nb), 0, 1))
    delta = cb-ca
    h_a, h_b = abs(float(delta @ na)), abs(float(delta @ nb))
    angle = math.degrees(math.acos(cosine))
    return {
        'minimum_selected_atom_distance_nm': float(np.linalg.norm(a[:, None]-b[None, :], axis=2).min()),
        'geometric_center_distance_nm': float(np.linalg.norm(delta)),
        'center_B_distance_to_plane_A_nm': h_a,
        'center_A_distance_to_plane_B_nm': h_b,
        'lateral_offset_relative_to_plane_A_nm': math.sqrt(max(0., float(delta @ delta)-h_a*h_a)),
        'acute_plane_angle_degrees': angle,
        'infinite_plane_distance_nm': h_a if angle < 1e-8 else 0.,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('docs/research/dimer-distance-examples'))
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    theta = np.arange(5)*2*math.pi/5
    a = np.column_stack([.14*np.cos(theta), .14*np.sin(theta), np.zeros(5)])
    angle = math.pi/3
    rotate = np.array([[math.cos(angle), 0, math.sin(angle)], [0, 1, 0],
                       [-math.sin(angle), 0, math.cos(angle)]])
    cases = [
        ('stacked', 'Facing regions', a+np.array([.05, 0, .35])),
        ('slipped', 'Large lateral offset', a+np.array([1.20, 0, .35])),
        ('tilted', 'Tilted regions', a@rotate.T+np.array([.05, 0, .35])),
    ]
    result = {'purpose': 'synthetic five-point geometry examples; not a molecule or production screening rule',
              'length_unit': 'nm', 'example_cutoff_nm': .6, 'cases': []}
    for name, title, b in cases:
        m = measures(a, b)
        # Check geometric invariance under a common rigid rotation and translation.
        other = measures(a@rotate.T+1., b@rotate.T+1.)
        for key, value in m.items():
            assert math.isclose(value, other[key], rel_tol=1e-10, abs_tol=2e-6), (name,key)
        result['cases'].append({'name': name, 'region_A_nm': a.tolist(), 'region_B_nm': b.tolist(),
                                **m, 'passes_example_min_atom_cutoff': m['minimum_selected_atom_distance_nm'] <= .6})
    assert math.isclose(result['cases'][0]['center_B_distance_to_plane_A_nm'], .35, abs_tol=1e-14)
    assert math.isclose(result['cases'][1]['center_B_distance_to_plane_A_nm'], .35, abs_tol=1e-14)
    assert not result['cases'][1]['passes_example_min_atom_cutoff']
    assert math.isclose(result['cases'][2]['acute_plane_angle_degrees'], 60., abs_tol=1e-12)
    assert result['cases'][2]['infinite_plane_distance_nm'] == 0.
    (args.output/'measurements.json').write_text(json.dumps(result, indent=2)+'\n')

    os.environ.setdefault('MPLCONFIGDIR', '/tmp/cube-mpl-cache')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(12, 5.7), constrained_layout=True)
    panels = fig.add_gridspec(2, 3, height_ratios=[4, 1])
    for number, ((name, title, b), row) in enumerate(zip(cases, result['cases']), 1):
        ax = fig.add_subplot(panels[0, number-1], projection='3d')
        for pts, color, label in [(a, '#1677b8', 'Region A'), (b, '#d26316', 'Region B')]:
            closed = np.vstack([pts, pts[:1]])
            ax.plot(*closed.T, color=color, lw=1.7, label=label)
            ax.scatter(*pts.T, color=color, s=24)
        for i, pt in enumerate(a, 1): ax.text(*pt, str(i), fontsize=7, color='#125580')
        ax.set(xlabel='x (nm)', ylabel='y (nm)', zlabel='z (nm)', title=title)
        if name == 'slipped':
            ax.set_xlim(-.2, 1.4); ax.set_ylim(-.25, .25); ax.set_zlim(-.05, .55)
            ax.set_box_aspect((1.6, .5, .6))
        else:
            ax.set_xlim(-.25, .35); ax.set_ylim(-.25, .25); ax.set_zlim(-.05, .55)
            ax.set_box_aspect((.6, .5, .6))
        ax.view_init(elev=18, azim=-65)
        notes = fig.add_subplot(panels[1, number-1])
        notes.axis('off')
        notes.text(.02, .95,
                  f"Nearest atoms: {row['minimum_selected_atom_distance_nm']:.3f} nm\n"
                  f"Centers: {row['geometric_center_distance_nm']:.3f} nm; tilt: {row['acute_plane_angle_degrees']:.0f}°\n"
                  f"Center B to plane A: {row['center_B_distance_to_plane_A_nm']:.3f} nm",
                  transform=notes.transAxes, fontsize=9, va='top')
    fig.suptitle('Distances between selected conjugated regions · synthetic geometry')
    fig.savefig(args.output/'geometry.png', dpi=160)
    fig.savefig(args.output/'geometry.svg')
    plt.close(fig)
    for row in result['cases']:
        print(row['name'], json.dumps({k:v for k,v in row.items() if not k.startswith('region_')}))


if __name__ == '__main__':
    main()
