"""KMC task dispatcher, source audit, independent snapshots and reusable artifacts."""
import csv
import hashlib
import json
from pathlib import Path
import platform
import time

from universal_kmc.quantum.cli import _write_result


def sha256(path):
    digest=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):digest.update(block)
    return digest.hexdigest()


def read_source(path,sources,json_data=True):
    path=Path(path).resolve()
    sources[str(path)]=sha256(path)
    return json.loads(path.read_text()) if json_data else None


def run_network(data,settings,directory,sources,*,resume=False):
    import numpy as np
    from . import network,simulation,assembly
    from universal_kmc import core,rates
    from .network import build_rates,exact_diffusion,validate_network
    from .simulation import simulate,analyze
    validate_network(data,settings['carriers'])
    directory=Path(directory);output=directory/'result.json'
    artifact_names=['network.json','result.json']+[carrier+suffix for carrier in settings['carriers']
        for suffix in ('_rates.json','.kmc','_samples.npz','_traces.json','_observations.csv')]
    if any(str((directory/name).resolve()) in sources for name in artifact_names):
        raise ValueError('KMC output artifacts must not overwrite a source file; choose another output directory')
    request={'algorithm':'cube-static-single-carrier-1','settings':{k:v for k,v in settings.items() if k!='output'},
             'network':data,'sources':sources,'numpy_version':np.__version__,'python_version':platform.python_version(),
             'implementation_sha256':{Path(m.__file__).name:sha256(m.__file__) for m in (network,simulation,assembly,core,rates)},
             'runner_sha256':sha256(__file__)}
    key=hashlib.sha256(json.dumps(request,sort_keys=True,allow_nan=False).encode()).hexdigest()
    if output.exists():
        old=json.loads(output.read_text())
        if not resume or old.get('fingerprint')!=key or old.get('status') not in ('completed','partial'):
            raise ValueError('existing KMC task; use --resume with identical finished inputs or select a new output directory')
        for relative,expected in old['output_sha256'].items():
            path=directory/relative
            if not path.is_file() or sha256(path)!=expected:raise ValueError('cannot resume: KMC output missing/changed')
        return {**old,'output':str(output),'reused':True}
    if any((directory/name).exists() for name in artifact_names):
        raise ValueError('untracked/partial KMC artifacts already exist; choose another output directory')
    directory.mkdir(parents=True,exist_ok=True)
    # Validate all rates/connectivity before writing scientific output.
    prepared={carrier:build_rates(data,carrier,settings['temperature']) for carrier in settings['carriers']}
    record={'schema':'cube-charge-transport-1','status':'running','fingerprint':key,
            'scope':'static_snapshot_zero_field_independent_single_carrier_Marcus',
            'input':request,'carriers':{},'limitations':['no exclusion or carrier-carrier interactions',
                'static snapshot, not a time-dependent trajectory network','no outer-sphere relaxation inferred',
                'CPU Python application backend; GPU and native material ensembles are not yet connected']}
    _write_result(output,record);artifacts=[];started=time.perf_counter()
    try:
        _write_result(directory/'network.json',data);artifacts.append('network.json')
        for carrier,(ids,edges,rows,pi,rank) in prepared.items():
            name=f'{carrier}_rates.json';_write_result(directory/name,rows);artifacts.append(name)
            name=f'{carrier}.kmc'
            with (directory/name).open('w') as stream:
                stream.write(f'KMC1 graph\n{len(ids)} {len(edges)} 0\n')
                for e in edges:stream.write(f'{e.source} {e.target} {e.rate:.17g} '+' '.join(f'{v:.17g}' for v in e.displacement)+'\n')
            artifacts.append(name)
            reference=exact_diffusion(len(ids),edges,pi,settings['exact-max-sites'])
            if reference['status']=='completed':
                reference['einstein_mobility_tensor_cm2_Vs']=(np.array(reference['diffusion_tensor_m2_s'])*1e4/rates.KB/settings['temperature']).tolist()
            grid,positions,traces,counts,capped=simulate(len(ids),edges,pi,settings,carrier)
            name=f'{carrier}_samples.npz';np.savez_compressed(directory/name,time_s=grid,displacement_m=positions);artifacts.append(name)
            name=f'{carrier}_traces.json';_write_result(directory/name,traces);artifacts.append(name)
            estimate=({'status':'partial','reason':'event_cap_reached','capped_trajectories':capped}
                      if capped else analyze(grid,positions,settings['temperature'],settings['fit-start'],settings['relative-tolerance'],rank))
            if 'covariance_m2' in estimate:
                name=f'{carrier}_observations.csv'
                with (directory/name).open('w',newline='') as f:
                    writer=csv.writer(f);writer.writerow(['time_s','mean_dx_m','mean_dy_m','mean_dz_m','cov_xx_m2','cov_yy_m2','cov_zz_m2','cov_xy_m2','cov_xz_m2','cov_yz_m2'])
                    for t,mean,cov in zip(grid,estimate.pop('mean_displacement_m'),estimate.pop('covariance_m2')):
                        writer.writerow([t,*mean,cov[0][0],cov[1][1],cov[2][2],cov[0][1],cov[0][2],cov[1][2]])
                artifacts.append(name)
            record['carriers'][carrier]={'status':estimate['status'],'site_ids':ids,'directed_channels':len(edges),
                'zero_coupling_channels':sum(r['zero_coupling'] for r in rows),'underflow_channels':sum(r['underflow'] for r in rows),
                'periodic_transport_rank':rank,'exact_reference':reference,'estimate':estimate,
                'trajectories':settings['trajectories'],'duration_s':settings['duration'],
                'total_events':sum(counts),'events_per_trajectory':counts,
                'sampler':'universal_kmc.Engine + outgoing edges + independent SHA256-seeded random.Random streams',
                'initial_distribution':'zero-field Boltzmann equilibrium over the connected network'}
            _write_result(output,record)
        record.update(status='completed' if all(v['status']=='completed' for v in record['carriers'].values()) else 'partial',
                      elapsed_seconds=time.perf_counter()-started,reused=False,
                      output_sha256={name:sha256(directory/name) for name in artifacts})
        _write_result(output,record)
    except Exception as exc:
        record.update(status='failed',error={'type':type(exc).__name__,'message':str(exc)})
        _write_result(output,record);raise
    return {**record,'output':str(output)}


def run_task(settings,preparation_report,extracted,*,resume=False):
    from .assembly import assemble
    jobs=[]
    if 'input' in settings:
        sources={};data=read_source(settings['input'],sources)
        jobs.append(('network',data,sources))
    else:
        from universal_kmc.workflow.runner import _dataset_sources
        manifest=settings.get('manifest') or extracted
        if manifest is None:raise ValueError('KMC has no extraction manifest')
        paths=_dataset_sources({'input':str(manifest)},None)
        if len(paths)>1 and any(k in settings for k in ('couplings','reorganization')):
            raise ValueError('multiple frames need matching workflow outputs per frame; external source files must be run frame by frame')
        for path in paths:
            sources={};data=read_source(path,sources)
            sources[str(Path(manifest).resolve())]=sha256(manifest)
            label=f"frame{data['frame']:04d}" if 'frame' in data else 'dataset'
            def task(name):
                result=preparation_report['tasks'].get(name,{})
                return result.get('frames',{}).get(label,{}) if 'frames' in result else result
            coupling=settings.get('couplings') or task('two-mol').get('csv')
            reorganization=settings.get('reorganization') or task('reorganization-energy').get('outputs')
            if not coupling or not reorganization:raise ValueError('KMC is missing completed couplings or reorganization results')
            network=assemble(path,coupling,reorganization,settings,sources)
            jobs.append((label,network,sources))
    results={}
    for label,data,sources in jobs:
        directory=Path(settings['output'])/(label if len(jobs)>1 else '')
        result=run_network(data,settings,directory,sources,resume=resume)
        results[label]={'status':result['status'],'output':result['output'],'reused':result.get('reused',False)}
    return {'status':'completed' if all(r['status']=='completed' for r in results.values()) else 'partial',
            'scope':'independent_static_snapshots','count':len(results),'outputs':results}
