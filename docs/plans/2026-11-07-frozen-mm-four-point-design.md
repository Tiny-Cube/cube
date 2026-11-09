# Frozen-MM electrostatic four-point reorganization

User request: implement QM/MM four-point reorganization with all MM atoms frozen. Reuse cube's CPU/GPU SCF and geometry workflow; keep vacuum inputs valid.

## Model and alternatives

Implement PySCF/GPU4PySCF fixed point-charge electrostatic embedding. Both electronic and nuclear QM–MM interactions contribute to every energy and QM gradient. This gives a concrete CPU/GPU implementation with analytic gradients. A full additive force-field model would also need QM–MM van der Waals parameters; polarizable embedding additionally needs an environment response model. Those are separate models and are not silently approximated here. MM-only energy is constant and cancels in each relaxation difference; omitted QM–MM van der Waals energy is not constant and must be explicitly reported as absent.

For one target QM molecule, choose a finite cluster once at the initial MD geometry. Include whole neighboring molecules within mm-cutoff (Angstrom), using consistent whole-molecule triclinic images. Include noncentral periodic replicas, but never central QM atoms. There is no Ewald sum or charge re-fitting. Alternatively load explicit x y z q point charges (Angstrom, elementary charge) in the same coordinate frame as a standalone QM XYZ. No covalent QM/MM boundary/link atoms are implemented.

## Interface and data flow

`$reorganization-energy`: `environment = qmmm`, `mm-frozen = true`, either `mm-source = manifest` with `mm-cutoff = 12` or `mm-charges = environment.pc`. Read GROMACS source coordinates/charges through the extraction manifest, verify original hashes and atom mapping, and translate the MM coordinates into the centered QM frame. Build the selected cluster once, never rebuild during optimization. An explicit file permits other MD engines after exporting charges and coordinates.

Optimize R0 and Rq with the same frozen environment M. Compute E0(R0;M), Eq(R0;M), Eq(Rq;M), E0(Rq;M), then the existing signed four-point formula. There are no MM coordinates in the optimization vector. Verify the actual SCF MM coordinates and charges at each evaluation, snapshot/fingerprint the environment, and include it in all resume keys. Embedded checkpoints use a distinct marker and cannot be translated or consumed as vacuum transfer checkpoints.

## Validation

Test parsing and invalid configurations; periodic/full-molecule selection and exclusion of central QM charges; real embedded Hamiltonian and finite-difference QM gradients; nonzero environment effects, common-translation invariance, resume invalidation on MM edits, and rejection of embedded checkpoints by vacuum consumers. Run a small CPU/GPU four-point case with fixed MM snapshots, both carriers, failure handling and repeat/resume checks. Preserve actual hardware provenance and do not describe the full material workload as completed.

Sources: https://pyscf.org/user/qmmm.html and installed PySCF/GPU4PySCF qmmm/itrf.py (PySCF 2.15.0, GPU4PySCF 1.9.1).
