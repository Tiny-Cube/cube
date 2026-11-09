#!/usr/bin/env python3
"""Stage 1: arbitrary monomer/dimer XYZ files -> PySCF .chk and .log files.

Example (each XYZ keeps its coordinates in the same dimer frame):
  uv run --extra quantum python scripts/run_pyscf.py A.xyz B.xyz AB.xyz \
      --output-dir qc_results --basis def2-svp --xc b3lyp --backend cpu

Use --extra quantum-cuda13 and --backend gpu for the validated CUDA 13 stack.
This script does NOT compute transfer integrals or optimize structures.
"""
import argparse
from pathlib import Path
import sys

from universal_kmc.quantum.config import SCFConfig
from universal_kmc.quantum.control import validate_output
from universal_kmc.quantum.geometry import read_xyz


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('xyz', nargs='+', type=Path, help='one or more monomer/dimer XYZ files')
    parser.add_argument('--output-dir', required=True, type=Path)
    parser.add_argument('--backend', choices=('cpu', 'gpu'), default='cpu')
    parser.add_argument('--basis', default='def2-svp')
    parser.add_argument('--xc', default='b3lyp', help='RKS functional or hf for RHF')
    parser.add_argument('--auxbasis', default='def2-universal-jkfit')
    parser.add_argument('--unit', choices=('Angstrom', 'Bohr'), default='Angstrom')
    parser.add_argument('--threads', type=int, default=4)
    parser.add_argument('--grid-level', type=int, default=3)
    parser.add_argument('--conv-tol', type=float, default=1e-10)
    parser.add_argument('--max-cycle', type=int, default=100)
    parser.add_argument('--max-memory-mb', type=int, default=4000)
    args = parser.parse_args(argv)
    try:
        config = SCFConfig(**{key: value for key, value in vars(args).items() if key not in ('xyz', 'output_dir')})
        config.validate()
        sources = [path.expanduser().resolve() for path in args.xyz]
        if len({p.stem for p in sources}) != len(sources):
            raise ValueError('XYZ basenames must be unique because they name the output files')
        atoms = [read_xyz(path) for path in sources]
        destinations = [[validate_output((args.output_dir.expanduser()/f'{path.stem}{suffix}').resolve(), sources)
                         for suffix in ('.chk', '.log', '.status.json')] for path in sources]
        from pyscf import gto, lib
        from universal_kmc.quantum.checkpoints import save_transfer_checkpoint
        from universal_kmc.quantum.cli import _write_result
        from dataclasses import asdict
        lib.num_threads(config.threads)
        for source, geometry, (checkpoint, log, status) in zip(sources, atoms, destinations):
            # An old successful checkpoint must never survive a failed rerun.
            checkpoint.unlink(missing_ok=True)
            checkpoint.parent.mkdir(parents=True, exist_ok=True)
            record = {'source': str(source), 'checkpoint': str(checkpoint), 'log': str(log),
                      'settings': asdict(config), 'status': 'running'}
            _write_result(status, record)
            try:
                mol = gto.M(atom=geometry, basis=config.basis, unit=config.unit,
                            charge=0, spin=0, cart=False, output=str(log), verbose=4,
                            max_memory=config.max_memory_mb)
                mf = mol.RHF() if config.xc.lower() == 'hf' else mol.RKS(xc=config.xc)
                mf = mf.density_fit(auxbasis=config.auxbasis)
                mf.conv_tol, mf.max_cycle, mf.direct_scf_tol = config.conv_tol, config.max_cycle, 1e-13
                if hasattr(mf, 'grids'):
                    mf.grids.level = config.grid_level
                if config.backend == 'gpu':
                    import gpu4pyscf
                    mf = mf.to_gpu()
                mf.kernel()
                # Reusable output template for your own PySCF scripts:
                # after a converged mf.kernel(), call this function.
                save_transfer_checkpoint(mf, checkpoint)
                record.update(status='completed', energy_hartree=float(mf.e_tot), converged=True,
                              source_class=type(mf).__module__+'.'+type(mf).__name__)
                _write_result(status, record)
                print(f'{source.name} -> {checkpoint}')
            except (ImportError, ValueError, RuntimeError, OSError) as exc:
                checkpoint.unlink(missing_ok=True)
                record.update(status='failed', error={'type': type(exc).__name__, 'message': str(exc)})
                _write_result(status, record)
                raise
        return 0
    except (ImportError, ValueError, RuntimeError, OSError) as exc:
        print(f'PySCF calculation failed: {exc}', file=sys.stderr)
        return 1


if __name__ == '__main__':
    raise SystemExit(main())
