# Eight periodic naphthalene MD frames

The initial crystal is the 298 K, 0.4 GPa structure used in the parent directory. GFN-FF optimization with xTB 6.7.1 relaxed the atom positions and lattice, changing the volume from 337.65 to 301.99 Å³. `relaxed.coord` contains the optimized 36-atom cell. The subsequent fixed-volume run used a Berendsen thermostat with a 298 K target, a 0.5 fs timestep, unscaled hydrogen masses, and 400 steps. The measured mean kinetic temperature over this short run was 327 K. Use the archived frames to reproduce the electronic calculation. A new MD run draws new initial velocities. [Experimental coordinates](https://doi.org/10.1107/S0108768106026814) · [GFN-FF](https://doi.org/10.1002/anie.202004239) · [Thermostat](https://doi.org/10.1063/1.448118)

The trajectory was written every 5 fs. `frames.json` selects snapshot indices 4, 9, …, 39 at 0.020–0.195 ps. Each manifest gives the complete lattice and both molecules. The SHA-256 of the original `xtb.trj` is `0ac957d73a58ab022845e876d88f8d5ee4eeda94ce4ca02e5b010575235ff14a`. xTB generated the force-field topology from the optimized structure. It issued a missing saved-topology warning during startup and completed the run.

From the repository root, run

```bash
wren-dos examples/dos/naphthalene-md/dos.toml
wren-dos examples/dos/naphthalene-md/dos.toml --block-study 1 2 4
```

The control uses periodic PBE/GTH-SZV, the internal GTH cell-mean reference, and a 2×1×1 mesh. The primary output is an eight-frame trajectory average, with four two-frame blocks. It demonstrates the calculation and error estimator on a short, correlated piece of dynamics. Equilibrium material spectra require convergence with sampling, cell size, basis and k mesh.
