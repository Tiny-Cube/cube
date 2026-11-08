#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2.0,<2.6"]
# ///
"""Standalone local GRO cluster viewer; uv run view_clusters.py structure.gro.

Only NumPy is required. With GRO alone, one contiguous residue is treated as
one molecule; bonds are inferred for display. Optional --tpr uses actual bonded
fragments (requires MDAnalysis: uv run --with MDAnalysis view_clusters.py ...).
Selection always uses all atoms, including hydrogen, and whole-molecule PBC
images. No QM/MM energies, charges or force-field parameters are inferred.
"""
import argparse
from collections import OrderedDict, deque
from dataclasses import dataclass
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import itertools
import json
import math
from pathlib import Path
import re
import sys
import threading
from urllib.parse import parse_qs, urlsplit
import webbrowser

import numpy as np


@dataclass
class Molecule:
    indices: np.ndarray
    coordinates: np.ndarray
    elements: list
    bonds: list
    label: str


def read_gro(filename):
    """Parse the first GRO frame, retaining record identities; lengths -> Å."""
    path=Path(filename)
    lines=path.read_text(encoding='utf-8-sig').splitlines()
    if len(lines)<3:
        raise ValueError('GRO needs a title, atom count, atom records and box')
    try:
        count=int(lines[1].strip())
    except ValueError as exc:
        raise ValueError('invalid GRO atom count') from exc
    if count<1 or len(lines)<count+3:
        raise ValueError('GRO atom count is zero or its records are incomplete')
    if any(line.strip() for line in lines[count+3:]):
        raise ValueError('provide a single GRO frame (concatenated frames are not supported)')
    atoms=[];blocks=[];last=None
    for index,line in enumerate(lines[2:count+2]):
        try:
            residue=int(line[:5]);resname=line[5:10].strip();name=line[10:15].strip();serial=int(line[15:20])
            dots=[m.start() for m in re.finditer(r'\.',line[20:])]
            if len(dots)<3:
                raise ValueError('coordinates must be fixed-width floating-point numbers')
            width=dots[1]-dots[0]
            if width<5 or dots[2]-dots[1]!=width or len(line)<20+3*width:
                raise ValueError('invalid coordinate precision/field width')
            xyz=[float(line[20+k*width:20+(k+1)*width])*10 for k in range(3)]
            if not name or not resname or not np.isfinite(xyz).all():
                raise ValueError('empty labels or nonfinite coordinates')
        except (ValueError,IndexError) as exc:
            raise ValueError(f'GRO atom line {index+3}: {exc}') from exc
        key=(residue,resname)
        if key!=last:
            blocks.append([]);last=key
        blocks[-1].append(index)
        atoms.append({'residue':residue,'resname':resname,'name':name,'serial':serial,
                      'resindex':len(blocks)-1,'coordinates':xyz})
    try:
        box=np.array([float(v) for v in lines[count+2].split()],dtype=float)
    except ValueError as exc:
        raise ValueError('invalid GRO box') from exc
    if len(box) not in (3,9) or not np.isfinite(box).all():
        raise ValueError('GRO box must contain 3 or 9 finite numbers')
    cell=np.diag(box[:3])
    if len(box)==9:
        cell[0,1],cell[0,2],cell[1,0],cell[1,2],cell[2,0],cell[2,1]=box[3:]
    cell*=10
    if not np.any(cell):
        cell=None
    elif np.linalg.det(cell)<=1e-8:
        raise ValueError('GRO periodic box must be right-handed and nondegenerate')
    return {'path':str(path.resolve()),'filename':path.name,'title':lines[0],
            'atoms':atoms,'blocks':blocks,'cell':cell}


def minimum_images(delta,cell):
    """Exact closest images, also for non-reduced triclinic cells."""
    delta=np.asarray(delta,dtype=float)
    if cell is None:
        return delta.copy()
    inverse=np.linalg.inv(cell)
    base=delta-np.rint(delta@inverse)@cell
    best=base.copy();best2=np.sum(best*best,axis=-1)
    radius=float(np.sqrt(best2.max()))
    extent=np.ceil(radius*np.linalg.norm(inverse,axis=0)+.5+1e-12).astype(int)
    if math.prod(2*int(n)+1 for n in extent)>100000:
        raise ValueError('periodic box is too skewed for a reliable image search')
    for image in itertools.product(*(range(-int(n),int(n)+1) for n in extent)):
        if image==(0,0,0):continue
        moved=base+np.array(image)@cell
        r2=np.sum(moved*moved,axis=-1)
        use=r2<best2-1e-14
        best[use]=moved[use];best2[use]=r2[use]
    return best


def unwrap_molecule(points,cell,bonds=None):
    """Lift a full molecule using bonded edges, or a geometric minimum tree."""
    points=np.asarray(points,dtype=float)
    if cell is None or len(points)==1:return points.copy()
    if bonds is None and len(points)>5000:
        raise ValueError('a residue has over 5000 atoms; provide --tpr to identify molecules or split the input')
    if bonds is not None:
        graph=[[] for _ in points]
        for a,b in bonds:graph[a].append(b);graph[b].append(a)
        result=np.empty_like(points);result[0]=points[0];seen={0};pending=deque([0])
        while pending:
            a=pending.popleft()
            neighbors=graph[a]
            if not neighbors:continue
            moves=minimum_images(points[neighbors]-points[a],cell)
            for b,move in zip(neighbors,moves):
                proposed=result[a]+move
                if b in seen:
                    if not np.allclose(result[b],proposed,rtol=0,atol=1e-6):
                        raise ValueError('bonded molecule winds through the periodic cell; no finite unwrapped cluster exists')
                else:
                    result[b]=proposed;seen.add(b);pending.append(b)
        if len(seen)!=len(points):raise ValueError('topology molecule is not connected')
        return result
    delta=minimum_images(points[None,:,:]-points[:,None,:],cell)
    distances=np.sum(delta*delta,axis=-1)
    result=np.empty_like(points);result[0]=points[0]
    visited=np.zeros(len(points),dtype=bool);visited[0]=True
    best=distances[0].copy();parents=np.zeros(len(points),dtype=int)
    for _ in range(1,len(points)):
        choices=np.where(visited,np.inf,best);b=int(np.argmin(choices));a=int(parents[b])
        result[b]=result[a]+delta[a,b];visited[b]=True
        improve=distances[b]<best
        best[improve]=distances[b,improve];parents[improve]=b
    return result


RADII={'H':.31,'C':.76,'N':.71,'O':.66,'F':.57,'P':1.07,'S':1.05,'Cl':1.02,
       'Br':1.20,'I':1.39,'Si':1.11,'B':.85,'Na':1.66,'Mg':1.41,'K':2.03,'Ca':1.76,'Fe':1.32,'Zn':1.22}


def infer_element(name,resname):
    letters=re.sub('[^A-Za-z]','',name)
    if letters in RADII:return letters
    if letters.upper() in ('CL','BR'):return letters.title()[:2]
    if letters.upper() in ('NA','MG','CA','ZN','FE') and resname.upper().strip('+-0123456789')==letters.upper():
        return letters.title()
    return letters[0].upper() if letters else 'X'


def inferred_bonds(points,elements):
    delta=points[:,None,:]-points[None,:,:]
    r2=np.sum(delta*delta,axis=-1)
    radii=np.array([RADII.get(e,.77) for e in elements])
    mask=np.triu((r2>.35**2)&(r2<(radii[:,None]+radii[None,:]+.4)**2),k=1)
    return np.argwhere(mask).tolist()


class Structure:
    def __init__(self,gro,*,tpr=None,pbc=True,max_cutoff=30.):
        self.data=read_gro(gro);self.cell=self.data['cell'] if pbc else None
        self.max_cutoff=float(max_cutoff);self.molecules=[]
        if not math.isfinite(self.max_cutoff) or self.max_cutoff<=0:
            raise ValueError('max-cutoff must be positive and finite')
        atoms=self.data['atoms'];positions=np.array([a['coordinates'] for a in atoms])
        topology=None
        if tpr is not None:
            try:import MDAnalysis as mda
            except ImportError as exc:raise ValueError('--tpr needs MDAnalysis; run uv run --with MDAnalysis scripts/view_clusters.py ...') from exc
            topology=mda.Universe(str(tpr))
            if len(topology.atoms)!=len(atoms) or list(topology.atoms.names)!=[a['name'] for a in atoms]:
                raise ValueError('TPR atom count/names do not match GRO')
            resindices=np.array([a['resindex'] for a in atoms])
            if not np.array_equal(topology.atoms.resindices,resindices):
                raise ValueError('TPR residue ordering does not match GRO')
            blocks=[m.indices.tolist() for m in topology.atoms.fragments]
            bonds=topology.bonds.indices
            all_elements=topology.atoms.elements.tolist()
            self.grouping='TPR 键拓扑：每个完整拓扑片段为一个分子'
            self.bond_source='TPR 真实化学键'
        else:
            blocks=self.data['blocks'];bonds=None
            all_elements=[infer_element(a['name'],a['resname']) for a in atoms]
            self.grouping='GRO 残基分组：假设每个残基是一个完整分子'
            self.bond_source='按距离推断的化学键，仅用于显示'
        for ids in blocks:
            lookup={int(g):i for i,g in enumerate(ids)}
            local_bonds=None if bonds is None else [[lookup[int(a)],lookup[int(b)]] for a,b in bonds if int(a) in lookup and int(b) in lookup]
            xyz=unwrap_molecule(positions[ids],self.cell,local_bonds)
            elements=[all_elements[i] for i in ids]
            if local_bonds is None:local_bonds=inferred_bonds(xyz,elements)
            residue_labels=list(dict.fromkeys(f"{atoms[i]['residue']}{atoms[i]['resname']}" for i in ids))
            label=', '.join(residue_labels[:3])+(' …' if len(residue_labels)>3 else '')
            self.molecules.append(Molecule(np.array(ids),xyz,elements,local_bonds,label))
        self._cache=OrderedDict();self._lock=threading.Lock()

    def metadata(self):
        return {'filename':self.data['filename'],'title':self.data['title'],'atom_count':len(self.data['atoms']),
                'molecule_count':len(self.molecules),'grouping':self.grouping,'bond_source':self.bond_source,
                'periodic':self.cell is not None,'max_cutoff':self.max_cutoff,
                'molecules':[{'id':i,'label':m.label,'atoms':len(m.indices)} for i,m in enumerate(self.molecules,1)]}

    def _candidates(self,center):
        with self._lock:
            if center in self._cache:
                self._cache.move_to_end(center);return self._cache[center]
        target=self.molecules[center-1].coordinates;c=target.mean(0)
        radius=np.linalg.norm(target-c,axis=1).max();result=[]
        inverse=None if self.cell is None else np.linalg.inv(self.cell)
        for index,molecule in enumerate(self.molecules,1):
            xyz=molecule.coordinates;m=xyz.mean(0)
            r=radius+np.linalg.norm(xyz-m,axis=1).max()+self.max_cutoff
            if self.cell is None:
                images=[(0,0,0)]
            else:
                fraction=(m-c)@inverse;bounds=r*np.linalg.norm(inverse,axis=0)
                lower=np.ceil(-fraction-bounds-1e-12).astype(int);upper=np.floor(-fraction+bounds+1e-12).astype(int)
                if math.prod(max(0,int(b-a+1)) for a,b in zip(lower,upper))>100000:
                    raise ValueError('too many periodic images; reduce max-cutoff or inspect box')
                images=itertools.product(*(range(int(a),int(b)+1) for a,b in zip(lower,upper)))
            for image in images:
                if index==center and image==(0,0,0):continue
                moved=xyz if self.cell is None else xyz+np.array(image)@self.cell
                if np.linalg.norm(moved.mean(0)-c)>r:continue
                minimum2=min(float(np.sum((target[start:start+64,None,:]-moved[None,:,:])**2,axis=-1).min())
                             for start in range(0,len(target),64))
                if minimum2<self.max_cutoff**2:
                    result.append((math.sqrt(minimum2),index,tuple(map(int,image))))
        result.sort(key=lambda v:(v[0],v[1],v[2]))
        with self._lock:
            self._cache[center]=result
            while len(self._cache)>8:self._cache.popitem(last=False)
        return result

    def cluster(self,center,cutoff):
        if type(center) is not int or not 1<=center<=len(self.molecules):
            raise ValueError('central molecule index is out of range')
        if not math.isfinite(cutoff) or not 0<=cutoff<=self.max_cutoff:
            raise ValueError(f'cutoff must be between 0 and {self.max_cutoff:g} Å')
        selected=[r for r in self._candidates(center) if r[0]<cutoff] if cutoff>0 else []
        if sum(len(self.molecules[i-1].indices) for d,i,image in selected)>100000:
            raise ValueError('cluster exceeds 100000 environment atoms; reduce cutoff')
        origin=self.molecules[center-1].coordinates.mean(0)
        coordinates=[];elements=[];names=[];serials=[];global_ids=[];groups=[];atom_groups=[];bonds=[]
        for distance,index,image in [(0.,center,(0,0,0)),*selected]:
            molecule=self.molecules[index-1];start=len(coordinates)
            shift=np.zeros(3) if self.cell is None else np.array(image)@self.cell
            coordinates.extend((molecule.coordinates+shift-origin).tolist());elements.extend(molecule.elements)
            atom_groups.extend([len(groups)]*len(molecule.indices));global_ids.extend((molecule.indices+1).tolist())
            names.extend(self.data['atoms'][i]['name'] for i in molecule.indices)
            serials.extend(self.data['atoms'][i]['serial'] for i in molecule.indices)
            bonds.extend([[start+a,start+b] for a,b in molecule.bonds])
            groups.append({'molecule':index,'label':molecule.label,'image':image,'start':start,
                           'count':len(molecule.indices),'minimum_distance':distance})
        return {'center':center,'cutoff':cutoff,'coordinates':coordinates,'elements':elements,'names':names,
                'serials':serials,'global_indices':global_ids,'groups':groups,'atom_groups':atom_groups,'bonds':bonds,
                'central_atoms':len(self.molecules[center-1].indices),'neighbor_atoms':len(coordinates)-len(self.molecules[center-1].indices),
                'neighbor_images':len(selected),'periodic_images':sum(any(image) for d,i,image in selected),
                'minimum_distance':selected[0][0] if selected else None,'coordinate_unit':'Angstrom'}


def xyz_text(cluster):
    return (f'{len(cluster["coordinates"])}\nGRO cluster; central={cluster["center"]}; cutoff={cluster["cutoff"]:g} Angstrom; '
            f'first {cluster["central_atoms"]} atoms central\n'+''.join(
                f'{e} {r[0]:.10f} {r[1]:.10f} {r[2]:.10f}\n' for e,r in zip(cluster['elements'],cluster['coordinates'])))


def make_server(structure,port=0,*,center=1,cutoff=6.):
    class Handler(BaseHTTPRequestHandler):
        def log_message(self,*args):pass

        def respond(self,data,content_type,status=200,filename=None):
            self.send_response(status);self.send_header('Content-Type',content_type)
            self.send_header('Cache-Control','no-store');self.send_header('X-Content-Type-Options','nosniff')
            self.send_header('Content-Length',str(len(data)))
            if filename:self.send_header('Content-Disposition',f'attachment; filename="{filename}"')
            self.end_headers()
            try:self.wfile.write(data)
            except (BrokenPipeError,ConnectionResetError):pass

        def do_GET(self):
            if self.headers.get('Host','').split(':')[0] not in ('127.0.0.1','localhost') or self.headers.get('Sec-Fetch-Site')=='cross-site':
                self.respond(b'Forbidden','text/plain',403);return
            url=urlsplit(self.path)
            try:
                if url.path=='/':
                    page=HTML.replace('__INITIAL__',json.dumps({'center':center,'cutoff':cutoff}))
                    self.respond(page.encode(),'text/html; charset=utf-8')
                elif url.path=='/api/structure':
                    self.respond(json.dumps(structure.metadata(),ensure_ascii=False).encode(),'application/json; charset=utf-8')
                elif url.path in ('/api/cluster','/api/cluster.xyz'):
                    query=parse_qs(url.query)
                    if set(query)-{'center','cutoff'} or any(len(v)!=1 for v in query.values()):raise ValueError('invalid query parameters')
                    central=int(query.get('center',[str(center)])[0]);distance=float(query.get('cutoff',[str(cutoff)])[0])
                    value=structure.cluster(central,distance)
                    if url.path.endswith('.xyz'):
                        self.respond(xyz_text(value).encode(),'chemical/x-xyz; charset=utf-8',filename=f'cluster_{central}_{distance:g}A.xyz')
                    else:self.respond(json.dumps(value,ensure_ascii=False,allow_nan=False,separators=(',',':')).encode(),'application/json; charset=utf-8')
                else:self.respond(b'Not found','text/plain',404)
            except (ValueError,OverflowError) as exc:
                self.respond(json.dumps({'error':str(exc)},ensure_ascii=False).encode(),'application/json; charset=utf-8',400)
    server=ThreadingHTTPServer(('127.0.0.1',port),Handler);server.daemon_threads=True
    return server


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('gro',type=Path)
    parser.add_argument('--tpr',type=Path)
    parser.add_argument('--cutoff',type=float,default=6.,help='initial all-atom minimum distance, Å')
    parser.add_argument('--max-cutoff',type=float,default=30.,help='maximum adjustable threshold, Å')
    parser.add_argument('--molecule',type=int,default=1,help='one-based central molecule index')
    parser.add_argument('--no-pbc',action='store_true')
    parser.add_argument('--port',type=int,default=0,help='0 selects a free local port')
    parser.add_argument('--open',action='store_true',help='open the local viewer in the default browser')
    args=parser.parse_args(argv)
    try:
        structure=Structure(args.gro,tpr=args.tpr,pbc=not args.no_pbc,max_cutoff=args.max_cutoff)
        if not 1<=args.molecule<=len(structure.molecules):raise ValueError('molecule index out of range')
        if not math.isfinite(args.cutoff) or not 0<=args.cutoff<=args.max_cutoff:raise ValueError('cutoff outside allowed range')
        if not 0<=args.port<=65535:raise ValueError('invalid port')
        server=make_server(structure,args.port,center=args.molecule,cutoff=args.cutoff)
    except (OSError,ValueError) as exc:
        print(f'view_clusters: {exc}',file=sys.stderr);return 1
    url=f'http://127.0.0.1:{server.server_port}/'
    print(f'{structure.data["filename"]}: {len(structure.data["atoms"])} atoms, {len(structure.molecules)} molecules',flush=True)
    print(f'{structure.grouping}; {structure.bond_source}',flush=True)
    print(url,flush=True);print('Ctrl+C to stop.',flush=True)
    if args.open:webbrowser.open(url)
    try:server.serve_forever()
    except KeyboardInterrupt:pass
    finally:server.server_close()
    return 0


HTML = r'''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>cube · GRO 簇查看器</title>
<style>
*{box-sizing:border-box}body{margin:0;background:#f4f7fa;color:#233949;font:14px/1.6 system-ui,sans-serif}header{padding:18px 24px 12px}h1{margin:0;font-size:25px}header p{margin:4px 0;color:#607585}.toolbar{padding:12px 24px;background:white;border-block:1px solid #dce5eb;display:flex;align-items:center;gap:12px;flex-wrap:wrap}.toolbar label{display:flex;align-items:center;gap:7px}button,select,input[type=number],.download{font:inherit;color:#233949;background:white;border:1px solid #b9cbd7;border-radius:6px;padding:6px 10px}button,select,.download{cursor:pointer}select{max-width:300px}input[type=range]{width:210px;accent-color:#db3357}input[type=number]{width:85px}.download{display:inline-block;text-decoration:none}button[aria-pressed=true]{color:white;background:#233949}button:focus-visible,input:focus-visible,select:focus-visible,a:focus-visible{outline:3px solid #70b9df;outline-offset:2px}main{display:grid;grid-template-columns:minmax(0,1fr) 300px;min-height:540px;height:calc(100vh - 205px)}.stage{position:relative;min-width:0;min-height:500px}canvas{width:100%;height:100%;display:block;cursor:grab;touch-action:none}.hint{position:absolute;bottom:12px;left:16px;font-size:12px;background:#ffffffd9;padding:5px 9px;border-radius:5px;pointer-events:none;color:#607585}aside{background:white;border-left:1px solid #dce5eb;padding:18px;overflow:auto}h2{font-size:18px;margin:0 0 8px}h3{font-size:14px;margin:18px 0 6px}.muted{color:#607585}.legend{display:flex;align-items:center;gap:8px;margin:4px 0}.dot{width:10px;height:10px;border-radius:50%}.neighbor{padding:8px 0;border-bottom:1px solid #edf1f5;font-size:13px}.mono{font:12px ui-monospace,monospace}.info{font-size:12px;color:#607585;background:#f4f7fa;padding:9px;border-radius:6px}.status{padding:5px 24px;font-size:12px;color:#607585;min-height:30px}.status.error{color:#a1253d;background:#fff1f2}#tip{position:absolute;background:#213747ed;color:white;font-size:12px;padding:8px 10px;border-radius:5px;pointer-events:none;display:none;max-width:270px}#loading{position:absolute;top:16px;left:18px;background:#fffffff0;padding:7px 13px;border-radius:6px;color:#607585;display:none}
@media(max-width:850px){main{grid-template-columns:1fr;height:auto}.stage{height:65vh}aside{border-left:0;border-top:1px solid #dce5eb;max-height:500px}.toolbar,header{padding-inline:15px}input[type=range]{width:150px}}
</style>
<header><h1>cube · GRO 簇查看器</h1><p id="file">正在读取结构…</p></header>
<div class="toolbar"><label>中央分子 <select id="center" aria-label="中央分子"></select></label><label>最短原子距离 <input id="slider" aria-label="距离滑块" type="range" min="0" max="30" step="0.1"><input id="cutoff" aria-label="距离阈值" type="number" min="0" max="30" step="0.1"> Å</label><button id="all" aria-pressed="true">完整环境</button><button id="near" aria-pressed="false">最近 4 个邻居</button><label><input id="hydrogen" type="checkbox">显示氢</label><label><input id="autofit" type="checkbox" checked>自动适配</label><button id="reset">重置视角</button><a id="export" class="download" href="#">导出簇 XYZ</a></div>
<div id="status" class="status" role="status" aria-live="polite">按两完整分子的全原子最短距离筛选，包含氢，不使用 core；阈值单位 Å。</div>
<main><div class="stage"><canvas id="canvas" aria-label="三维分子簇"></canvas><div id="loading">正在计算全原子邻域…</div><div class="hint">拖动旋转 · 滚轮缩放 · Shift＋拖动平移 · 悬停查看原子</div><div id="tip"></div></div><aside><h2 id="heading">中央分子</h2><div id="stats" class="muted"></div><h3>颜色</h3><div class="legend"><i class="dot" style="background:#db3357"></i>中央完整分子</div><div class="legend"><i class="dot" style="background:#7b9aa8"></i>邻居镜像</div><h3>分子与化学键来源</h3><div id="source" class="info"></div><h3>筛选定义</h3><div class="info">入选条件：两个完整分子的所有原子对中，最短距离严格小于当前阈值。氢的显示开关不改变筛选结果。周围分子整体保留，所以簇外缘可以超过阈值。</div><h3>最近邻居</h3><div id="neighbors"></div></aside></main>
<script>
const INITIAL=__INITIAL__;
const $=id=>document.getElementById(id),canvas=$('canvas'),ctx=canvas.getContext('2d');
let structure=null,data=null,mode='all',hydrogen=false,rx=.38,ry=-.5,zoom=1,panX=0,panY=0,points=[],visible=[],pending=false,W=0,H=0,radius=1,request=0,controller=null,timer=null;
function error(message){$('status').textContent=message;$('status').classList.add('error')}
function esc(value){return String(value).replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]))}
async function loadCluster(){
 const center=Number($('center').value),cutoff=Number($('cutoff').value);
 if(!Number.isFinite(cutoff)||cutoff<0||cutoff>structure.max_cutoff){error(`阈值需要在 0–${structure.max_cutoff} Å 之间`);return}
 const sequence=++request;if(controller)controller.abort();controller=new AbortController();$('loading').style.display='block';
 try{const query=new URLSearchParams({center:String(center),cutoff:String(cutoff)}),response=await fetch('/api/cluster?'+query,{signal:controller.signal});const value=await response.json();if(!response.ok)throw new Error(value.error||'无法读取簇');if(sequence!==request)return;
 data=value;if($('autofit').checked){zoom=1;panX=panY=0}metadata();schedule();
 $('status').classList.remove('error');$('status').textContent=`当前阈值 ${data.cutoff} Å · 全原子严格最短距离（含氢，不使用 core） · ${structure.periodic?'周期镜像已启用':'非周期模式'}`;
 }catch(e){if(e.name!=='AbortError'&&sequence===request)error('未更新结构：'+e.message)}finally{if(sequence===request)$('loading').style.display='none'}
}
function metadata(){
 $('heading').textContent=`分子 ${data.center} · ${data.groups[0].label}`;
 $('stats').innerHTML=`中央 ${data.central_atoms} 个原子<br>邻居 ${data.neighbor_atoms} 个原子 / ${data.neighbor_images} 个完整分子镜像<br>其中 ${data.periodic_images} 个非零周期镜像<br>最近全原子距离：<b>${data.minimum_distance===null?'无邻居':data.minimum_distance.toFixed(3)+' Å'}</b><br><span id="visible-count"></span>`;
 $('neighbors').innerHTML=data.groups.slice(1,21).map(g=>`<div class="neighbor">分子 ${g.molecule} · ${esc(g.label)} · ${g.minimum_distance.toFixed(3)} Å<br><span class="mono">image [${g.image.join(', ')}]</span></div>`).join('')||(data.cutoff===0?'仅显示中央分子':'该阈值下没有邻居');
 if(data.neighbor_images>20)$('neighbors').insertAdjacentHTML('beforeend',`<div class="muted">列出最近 20 个；共 ${data.neighbor_images} 个</div>`);
 $('export').href='/api/cluster.xyz?'+new URLSearchParams({center:String(data.center),cutoff:String(data.cutoff)});
}
function schedule(){if(!pending){pending=true;requestAnimationFrame(()=>{pending=false;draw()})}}
function draw(){
 const bounds=canvas.getBoundingClientRect();W=bounds.width;H=bounds.height;const dpr=Math.min(devicePixelRatio||1,2);
 if(canvas.width!==Math.round(W*dpr)||canvas.height!==Math.round(H*dpr)){canvas.width=Math.round(W*dpr);canvas.height=Math.round(H*dpr)}ctx.setTransform(dpr,0,0,dpr,0,0);ctx.fillStyle='#f4f7fa';ctx.fillRect(0,0,W,H);if(!data)return;
 visible=data.elements.map((e,i)=>(hydrogen||e!=='H')&&(mode==='all'||data.atom_groups[i]<=4));let fit=1;for(let i=0;i<data.coordinates.length;i++)if(visible[i])fit=Math.max(fit,Math.hypot(...data.coordinates[i]));if($('autofit').checked||!radius)radius=fit;
 const scale=Math.min(W,H)*.44/radius*zoom,cx=Math.cos(rx),sx=Math.sin(rx),cy=Math.cos(ry),sy=Math.sin(ry);
 function rotate(p){const x=p[0]*cy+p[2]*sy,z=-p[0]*sy+p[2]*cy;return [x,p[1]*cx-z*sx,p[1]*sx+z*cx]}
 points=data.coordinates.map(p=>{const r=rotate(p);return {x:W/2+panX+r[0]*scale,y:H/2+panY-r[1]*scale,z:r[2]}});
 const primitives=[];for(const [a,b] of data.bonds)if(visible[a]&&visible[b])primitives.push({kind:0,a,b,z:(points[a].z+points[b].z)/2});for(let i=0;i<visible.length;i++)if(visible[i])primitives.push({kind:1,i,z:points[i].z});primitives.sort((a,b)=>a.z-b.z);
 for(const p of primitives){const g=data.atom_groups[p.kind===0?p.a:p.i],color=g===0?'#db3357':'#7b9aa8';ctx.globalAlpha=g===0?1:mode==='all'?.42:.8;if(p.kind===0){const a=points[p.a],b=points[p.b];ctx.strokeStyle=color;ctx.lineWidth=g===0?2.2:mode==='all'?.65:1.1;ctx.beginPath();ctx.moveTo(a.x,a.y);ctx.lineTo(b.x,b.y);ctx.stroke()}else{const a=points[p.i],size=(g===0?3:mode==='all'?1.1:1.8)*(data.elements[p.i]==='H'?.65:1)*Math.sqrt(zoom);ctx.fillStyle=color;ctx.beginPath();ctx.arc(a.x,a.y,Math.max(.5,size),0,Math.PI*2);ctx.fill()}}
 ctx.globalAlpha=1;const base=Math.pow(10,Math.floor(Math.log10(90/scale)));let bar=base;for(const m of [1,2,5])if(base*m*scale<=110)bar=base*m;ctx.strokeStyle='#233949';ctx.lineWidth=2;ctx.beginPath();ctx.moveTo(20,32);ctx.lineTo(20+bar*scale,32);ctx.stroke();ctx.fillStyle='#233949';ctx.font='12px system-ui';ctx.fillText(`${bar} Å`,20,23);
 const origin=[W-60,H-65];for(const [label,p,color] of [['x',[1,0,0],'#aa4355'],['y',[0,1,0],'#448568'],['z',[0,0,1],'#4479a8']]){const r=rotate(p);ctx.strokeStyle=color;ctx.beginPath();ctx.moveTo(...origin);ctx.lineTo(origin[0]+r[0]*32,origin[1]-r[1]*32);ctx.stroke();ctx.fillStyle=color;ctx.fillText(label,origin[0]+r[0]*42-3,origin[1]-r[1]*42+3)}
 const count=visible.filter(Boolean).length;$('visible-count').textContent=`当前显示 ${count} 个原子${mode==='near'?'（中央＋最近至多 4 个邻居）':''}`;if(!count)ctx.fillText('所有原子均被隐藏，请启用“显示氢”',Math.max(20,W/2-130),H/2);
 canvas.dataset.rendered='true';canvas.dataset.visibleAtoms=String(count);canvas.dataset.cutoff=String(data.cutoff);canvas.dataset.neighborImages=String(data.neighbor_images);canvas.dataset.center=String(data.center);
}
function reset(){rx=.38;ry=-.5;zoom=1;panX=panY=0;$('tip').style.display='none';radius=0;schedule()}
$('center').onchange=()=>{reset();clearTimeout(timer);loadCluster()};
function changeDistance(event){if(event.target.id==='slider')$('cutoff').value=event.target.value;else if(Number.isFinite(Number(event.target.value)))$('slider').value=event.target.value;clearTimeout(timer);timer=setTimeout(loadCluster,180)}
$('slider').oninput=changeDistance;$('cutoff').oninput=changeDistance;$('cutoff').onchange=()=>{clearTimeout(timer);loadCluster()};
for(const m of ['all','near'])$(m).onclick=()=>{mode=m;for(const b of ['all','near'])$(b).setAttribute('aria-pressed',String(b===m));schedule()};
$('hydrogen').onchange=e=>{hydrogen=e.target.checked;schedule()};$('autofit').onchange=schedule;$('reset').onclick=reset;
let drag=null;canvas.onpointerdown=e=>{drag={x:e.clientX,y:e.clientY};canvas.setPointerCapture(e.pointerId);$('tip').style.display='none';canvas.style.cursor='grabbing'};canvas.onpointerup=()=>{drag=null;canvas.style.cursor='grab'};canvas.onpointercancel=canvas.onpointerup;
canvas.onpointermove=e=>{if(drag){const dx=e.clientX-drag.x,dy=e.clientY-drag.y;drag={x:e.clientX,y:e.clientY};if(e.shiftKey){panX+=dx;panY+=dy}else{ry+=dx*.008;rx+=dy*.008}schedule();return}if(!data)return;const box=canvas.getBoundingClientRect(),x=e.clientX-box.left,y=e.clientY-box.top;let closest=-1,distance=64;for(let i=0;i<points.length;i++)if(visible[i]){const d=(points[i].x-x)**2+(points[i].y-y)**2;if(d<distance){closest=i;distance=d}}if(closest<0){$('tip').style.display='none';return}const group=data.groups[data.atom_groups[closest]];$('tip').innerHTML=`${group.molecule===data.center&&data.atom_groups[closest]===0?'中央':'邻居'} · 分子 ${group.molecule} · ${esc(data.elements[closest])}<br>${esc(data.names[closest])} · 分子内原子 ${closest-group.start+1}<br>GRO 原子记录 ${data.global_indices[closest]} / 编号 ${data.serials[closest]}<br>image [${group.image.join(', ')}]`;$('tip').style.display='block';$('tip').style.left=Math.max(5,Math.min(x+12,W-260))+'px';$('tip').style.top=Math.max(5,y-85)+'px'};canvas.onpointerleave=()=>{$('tip').style.display='none'};
canvas.addEventListener('wheel',e=>{e.preventDefault();zoom=Math.max(.2,Math.min(10,zoom*Math.exp(-e.deltaY*.001)));schedule()},{passive:false});new ResizeObserver(schedule).observe(canvas);
async function boot(){try{const response=await fetch('/api/structure');if(!response.ok)throw new Error('无法读取 GRO');structure=await response.json();$('file').textContent=`${structure.filename} · ${structure.atom_count} 个原子 · ${structure.molecule_count} 个分子 · ${structure.periodic?'周期盒':'非周期结构'}`;$('source').textContent=structure.grouping+'；'+structure.bond_source+'。本页只显示几何，不从 GRO 推断电荷或力场参数。';for(const m of structure.molecules){const o=document.createElement('option');o.value=m.id;o.textContent=`${m.id} · ${m.label} · ${m.atoms} atoms`;$('center').append(o)}$('center').value=INITIAL.center;$('slider').max=$('cutoff').max=structure.max_cutoff;$('slider').value=$('cutoff').value=INITIAL.cutoff;await loadCluster()}catch(e){error(e.message)}}boot();
</script></html>
'''


if __name__=='__main__':
    raise SystemExit(main())
