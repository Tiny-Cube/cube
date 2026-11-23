import json, sys, time
from pathlib import Path
from wren.dos.runner import run
from wren.quantum.cli import _write_result
root=Path('/home/codex/wren-validation-20261123')
gpu_root=root/'gpu-check-1/md-monomer'
gpu=json.loads((gpu_root/'result.json').read_text())
settings={**gpu['request']['settings'],'qc-engine':'pyscf','dos-backend':'cpu','output':str(root/'md-cpu')}
cpu=run(settings)
def values(p):
 return json.loads((p/'spectra.json').read_text())['frames'][0]['molecules'][0]
a,b=values(root/'md-cpu'),values(gpu_root)
import numpy as np
error=float(abs(np.asarray(a['orbital_energies_ev'])-b['orbital_energies_ev']).max())
report={'schema':'wren-dos-md-device-comparison-1','status':'completed','same_host':True,'threads':8,'method':'HF/STO-3G, Weigend auxiliary fitting, conv_tol=1e-7',
'cpu_wall_seconds':cpu['elapsed_seconds'],'gpu_wall_seconds':gpu['elapsed_seconds'],'max_orbital_error_ev':error,
'cpu_total_energy_hartree':a['scf']['energy_hartree'],'gpu_total_energy_hartree':b['scf']['energy_hartree'],
'cpu_homo_ev':a['homo_ev'],'gpu_homo_ev':b['homo_ev'],'cpu_lumo_ev':a['lumo_ev'],'gpu_lumo_ev':b['lumo_ev'],
'note':'One repetition; GPU warmed by H2 checks; implementation timing only, not material convergence or statistical benchmark'}
_write_result(root/'md-comparison.json',report)
print(json.dumps(report,indent=2))
