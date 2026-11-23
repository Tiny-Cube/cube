"""Plot the archived zero-onsite comparison on each channel's relative axis."""
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


root = Path(__file__).resolve().parent
fig, axes = plt.subplots(1, 2, figsize=(9, 3.8), layout='constrained', sharey=True)
for ax, name in zip(axes, ('homo', 'lumo')):
    data = np.genfromtxt(root/f'{name}.csv', delimiter=',', names=True)
    ax.plot(data['energy_ev'], data['signed_k4'], label='Archived signs', color='#18558e', lw=2)
    ax.plot(data['energy_ev'], data['absolute_k4'], label='All-positive test', color='#ae5726',
            dashes=(3,2), lw=1.5)
    visible = data['signed_k4']+data['absolute_k4'] > 1.e-3
    ax.set_xlim(data['energy_ev'][visible].min(), data['energy_ev'][visible].max())
    ax.set_title(name.upper()+' manifold')
    ax.set_xlabel('Energy relative to channel zero (eV)')
    ax.grid(axis='y', alpha=.22)
    ax.spines[['top', 'right']].set_visible(False)
axes[0].set_ylabel('DOS (states / eV / molecule)')
axes[0].legend(frameon=False, fontsize=9)
fig.suptitle('100 molecules  ·  uniform site energies  ·  4 × 4 × 4 k points')
fig.savefig(root/'frontier.png', dpi=180)
