"""Validate KMC settings without importing NumPy or any QC/MD backend."""
import math
from pathlib import Path

FIELDS={'do','input','manifest','couplings','reorganization','output','carriers','temperature','duration',
        'trajectories','samples','max-events','seed','fit-start','relative-tolerance','trace-trajectories',
        'site-energy-model','pair-reorganization-model','exact-max-sites','backend'}


def normalize(settings,base):
    s=settings
    if set(s)-FIELDS:raise ValueError('unknown KMC settings')
    def real(name,default,positive=True):
        value=float(s.get(name,default))
        if not math.isfinite(value) or (positive and value<=0):raise ValueError(f'kmc {name} must be finite and positive')
        return value
    def integer(name,default,minimum=1):
        value=s.get(name,str(default))
        if not value.isdecimal() or int(value)<minimum:raise ValueError(f'kmc {name} must be an integer >= {minimum}')
        return int(value)
    n={'temperature':real('temperature',300),'duration':real('duration',1e-9),
       'trajectories':integer('trajectories',128),'samples':integer('samples',51,5),
       'max-events':integer('max-events',1000000),'seed':integer('seed',42,0),
       'fit-start':real('fit-start',.5),'relative-tolerance':real('relative-tolerance',.25),
       'trace-trajectories':integer('trace-trajectories',3,0),'exact-max-sites':integer('exact-max-sites',512),
       'backend':s.get('backend','python').lower(),
       'carriers':[v.strip().lower() for v in s.get('carriers','hole,electron').split(',')]}
    if n['backend']!='python':raise ValueError('material KMC currently supports backend=python; native generic KMC engines remain separately available')
    if n['fit-start']>=.9:raise ValueError('kmc fit-start must be >0 and <0.9')
    if not n['carriers'] or len(set(n['carriers']))!=len(n['carriers']) or not set(n['carriers'])<={'hole','electron'}:
        raise ValueError('kmc carriers must select hole and/or electron without duplicates')
    if n['trajectories']*n['samples']>10_000_000:raise ValueError('too many KMC observation points; reduce trajectories or samples')
    for key in ('input','manifest','couplings','reorganization'):
        if s.get(key,'').strip():
            path=(base/Path(s[key]).expanduser()).resolve()
            if not path.is_file():raise ValueError(f'kmc {key} file not found: {path}')
            n[key]=str(path)
    if 'input' in n and any(k in n for k in ('manifest','couplings','reorganization')):
        raise ValueError('kmc input network and manifest/couplings/reorganization sources are mutually exclusive')
    model=s.get('site-energy-model')
    if model is not None and model not in ('uniform','four-point-adiabatic'):
        raise ValueError('site-energy-model must be uniform or four-point-adiabatic')
    if model is not None:n['site-energy-model']=model
    n['pair-reorganization-model']=s.get('pair-reorganization-model','four-point-equivalent')
    if n['pair-reorganization-model']!='four-point-equivalent':
        raise ValueError('computed QC assembly currently supports four-point-equivalent only; supply explicit pair lambdas in a network for inequivalent sites')
    n['output']=str((base/Path(s.get('output','cube-results/charge-carrier')).expanduser()).resolve())
    if any(n['output']==n.get(key) for key in ('input','manifest','couplings','reorganization')):
        raise ValueError('KMC output must not overwrite a source file')
    return n


def preflight(settings,tasks):
    if 'input' in settings:return
    for source,module in (('manifest','get_dimers_monomers'),('couplings','two-mol'),('reorganization','reorganization-energy')):
        if source not in settings and module not in tasks:
            raise ValueError(f'$kmc needs {source}=existing_file or an enabled ${module} stage; no MD/QC tasks were started')
    if 'site-energy-model' not in settings:
        raise ValueError('$kmc QC assembly requires an explicit site-energy-model; no MD/QC tasks were started')
