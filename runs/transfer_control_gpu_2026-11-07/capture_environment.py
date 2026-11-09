"""Capture the tested wheel, source hashes and hardware without host credentials."""
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform
import subprocess
import sys

import cupy as cp
import numpy as np
import universal_kmc.quantum.checkpoints as checkpoints

root = Path.cwd()
output = Path(sys.argv[1])
installed = Path(checkpoints.__file__).parent
hashes = {}
for file in sorted((root/'src/universal_kmc/quantum').glob('*.py')):
    source_hash = hashlib.sha256(file.read_bytes()).hexdigest()
    installed_hash = hashlib.sha256((installed/file.name).read_bytes()).hexdigest()
    assert source_hash == installed_hash, file.name
    hashes[str(file.relative_to(root))] = source_hash
for name in ('scripts/run_pyscf.py', 'scripts/calculate_transfer.py', 'scripts/validate_transfer_control.py',
             'pyproject.toml', 'uv.lock', 'examples/quantum/ethylene_a.xyz',
             'examples/quantum/ethylene_b.xyz', 'examples/quantum/ethylene_ab.xyz'):
    hashes[name] = hashlib.sha256((root/name).read_bytes()).hexdigest()
properties = cp.cuda.runtime.getDeviceProperties(cp.cuda.runtime.getDevice())
name = properties['name']
hardware = {'device': name.decode() if isinstance(name, bytes) else str(name),
            'compute_capability': [properties['major'], properties['minor']],
            'memory_bytes': properties['totalGlobalMem'],
            'driver_api_version': cp.cuda.runtime.driverGetVersion(),
            'runtime_api_version': cp.cuda.runtime.runtimeGetVersion(),
            'nvidia_smi': subprocess.check_output(['nvidia-smi', '--query-gpu=name,driver_version,memory.total',
                                                  '--format=csv,noheader'], text=True).strip()}
wheel = root/'validation/control-20261107/universal_kmc-0.1.0-py3-none-any.whl'
distribution = importlib.metadata.distribution('universal-kmc')
report = {'python': platform.python_version(), 'platform': platform.platform(), 'hardware': hardware,
          'installed_package_path': str(installed),
          'installed_direct_url': json.loads(distribution.read_text('direct_url.json')),
          'wheel_sha256': hashlib.sha256(wheel.read_bytes()).hexdigest(),
          'installed_quantum_source_matches_uploaded_source': True, 'source_sha256': hashes,
          'versions': {d.metadata['Name']: d.version for d in importlib.metadata.distributions()}}

def couplings(result):
    return np.array([result['results'][c]['methods'][m]['scalar_magnitude_ev']
                     for c in ('hole', 'electron') for m in ('dipro', 'esid', 'orbital_gmh')])

comparisons = {}
for source in ('cpu', 'gpu'):
    cpu = json.loads((output/f'{source}-source-cpu-projection.json').read_text())
    gpu = json.loads((output/f'{source}-source-gpu-projection.json').read_text())
    differences = abs(couplings(cpu)-couplings(gpu))
    assert float(differences.max()) < 1e-11
    comparisons[source+'_checkpoints'] = {'maximum_cpu_gpu_projection_difference_ev': float(differences.max()),
                                         'couplings_ev': couplings(gpu).tolist()}
report['same_checkpoint_projection_comparisons'] = comparisons
(output/'environment.json').write_text(json.dumps(report, indent=2)+'\n')
print(json.dumps(report, indent=2))
