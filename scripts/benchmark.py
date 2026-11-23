#!/usr/bin/env python3
"""Reproducible serial benchmark of the existing four-language implementations."""
import argparse
from datetime import date
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import random
import re
import statistics
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
BUILD = ROOT / 'build/benchmark'
LANGS = ['python', 'cpp', 'rust', 'fortran']

def run(command, **kwargs):
    return subprocess.run(list(map(str, command)), cwd=ROOT, check=True, **kwargs)

def build():
    run(['cmake', '-S', ROOT, '-B', BUILD, '-DCMAKE_BUILD_TYPE=Release',
         '-DWREN_BUILD_PYTHON=OFF', '-DBUILD_TESTING=OFF'])
    run(['cmake', '--build', BUILD, '--config', 'Release', '--parallel', '--target', 'wren_benchmarks'])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=Path('runs/performance'))
    parser.add_argument('--repeats',type=int,default=5)
    parser.add_argument('--python',default=str(ROOT/'.venv/bin/python'))
    args=parser.parse_args()
    if args.repeats<2: parser.error('at least two repeats required')
    output=args.output.resolve(); output.mkdir(parents=True,exist_ok=True)
    if (output/'summary.json').exists(): parser.error('use a new output directory')
    build()
    bin_dir = BUILD/'bin'
    if (bin_dir/'Release').is_dir():
        bin_dir = bin_dir/'Release'
    affinity=sorted(os.sched_getaffinity(0)) if hasattr(os,'sched_getaffinity') else []
    if affinity: os.sched_setaffinity(0,{affinity[0]})
    env=os.environ.copy(); env['PYTHONPATH']=str(ROOT/'src')
    kernels={'python':[args.python,'benchmarks/kernel.py'],
             **{lang:[str(bin_dir/f'kernel-{lang}')] for lang in LANGS[1:]}}
    clis={'python':[args.python,'-m','wren'],
          'cpp':[str(bin_dir/'ukmc-cpp')], 'rust':[str(bin_dir/'ukmc-rust')],
          'fortran':[str(bin_dir/'ukmc-fortran')]}
    workloads=[('two_state','examples/two_state.kmc',200000),
               ('birth_death','examples/birth_death.kmc',200000),
               ('hole','runs/e_006_ether_12-17_0/hole.kmc',20000),
               ('electron','runs/e_006_ether_12-17_0/electron.kmc',20000)]
    def version(cmd): return run(cmd,capture_output=True,text=True).stdout.splitlines()[0]
    cpu=next((line.split(':',1)[1].strip() for line in Path('/proc/cpuinfo').read_text().splitlines()
              if line.startswith('model name')),platform.processor())
    commit = subprocess.run(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True, capture_output=True) if (ROOT/'.git').exists() else None
    result={'date':date.today().isoformat(),'cpu':cpu,'platform':platform.platform(),
            'pinned_cpu':affinity[0] if affinity else None,'repeats':args.repeats,'seed':42,
            'versions':{k:version(v) for k,v in {'python':[args.python,'--version'],'cpp':['c++','--version'],
                         'fortran':['gfortran','--version'],'rust':['rustc','--version']}.items()},
            'build':'CMake Release / Cargo release; see CMakeCache.txt for configured compilers and flags',
            'source_commit':commit.stdout.strip() if commit and commit.returncode == 0 else None,
            'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
                for folder in ['src','cpp','rust/src','fortran','benchmarks'] for p in (ROOT/folder).rglob('*')
                if p.is_file() and p.suffix in ['.py','.cpp','.hpp','.rs','.f90']},'results':[]}
    rng=random.Random(20261102)
    for name,model,steps in workloads:
        reference=None
        for mode,commands in [('kernel',kernels),('cli',clis)]:
            samples={lang:[] for lang in LANGS}
            # One full unrecorded warmup per language, then randomized, interleaved repeats.
            for rep in range(-1,args.repeats):
                order=LANGS.copy(); rng.shuffle(order)
                for lang in order:
                    command=[*commands[lang],model,str(steps),'42']
                    start=time.perf_counter()
                    if mode=='kernel':
                        proc=run(command,env=env,capture_output=True,text=True)
                        values=list(map(float,proc.stdout.split()))
                        duration=values[0]; signature=values[1:]
                        assert signature[0]==steps
                        if reference is None: reference=signature
                        assert len(signature)==len(reference)
                        assert all(math.isclose(a,b,rel_tol=3e-13,abs_tol=1e-25)
                                   for a,b in zip(signature,reference)),(name,lang,signature,reference)
                    else:
                        with open(os.devnull,'w') as sink:
                            proc=run(command,env=env,stdout=sink,stderr=subprocess.PIPE,text=True)
                        duration=time.perf_counter()-start
                        match=re.search(r'events=(\d+) time=\s*(\S+)',proc.stderr)
                        assert match and int(match[1])==steps,proc.stderr
                        assert math.isclose(float(match[2]),reference[1],rel_tol=3e-13)
                    assert duration>0
                    if rep>=0: samples[lang].append(duration)
                print(f'{name} {mode}: {"warmup" if rep<0 else str(rep+1)+"/"+str(args.repeats)} complete',flush=True)
            for lang in LANGS:
                med=statistics.median(samples[lang])
                result['results'].append({'model':name,'input':model,'events':steps,'mode':mode,'language':lang,
                    'seconds':samples[lang],'median_s':med,'min_s':min(samples[lang]),'max_s':max(samples[lang]),
                    'events_per_second':steps/med,'speedup_vs_python':statistics.median(samples['python'])/med,
                    'final_signature':reference,'correctness_passed':True})
            (output/'partial.json').write_text(json.dumps(result,indent=2)+'\n')
    (output/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    (output/'partial.json').unlink()
    print(json.dumps(result['results'],indent=2))

if __name__=='__main__': main()
