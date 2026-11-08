// Exercise the embedded browser renderer without installing browser dependencies.
// Canvas calls are counted, not rasterized: timings here are JS overhead, not FPS.
const vm=require('node:vm'),assert=require('node:assert/strict'),{performance}=require('node:perf_hooks');
function renderer(script,cluster){
 const nodes=new Map(),calls={},frames=new Map(),timers=new Map();let id=0;
 const ctx={};for(const name of ['setTransform','fillRect','beginPath','moveTo','lineTo','stroke','arc','rect','fill','fillText'])ctx[name]=()=>{calls[name]=(calls[name]||0)+1};
 function node(key){if(!nodes.has(key))nodes.set(key,{value:'',checked:key==='autofit',hidden:true,width:0,height:0,dataset:{},style:{},classList:{add(){},remove(){}},textContent:'',innerHTML:'',getContext:()=>ctx,getBoundingClientRect:()=>({width:1000,height:700,left:0,top:0}),addEventListener(type,fn){this[type]=fn},setAttribute(){},setPointerCapture(){},insertAdjacentHTML(){},append(){}});return nodes.get(key)}
 const sandbox={console,performance,URLSearchParams,document:{getElementById:node,createElement:()=>({})},devicePixelRatio:2,requestAnimationFrame:fn=>{frames.set(++id,fn);return id},cancelAnimationFrame:n=>frames.delete(n),setTimeout:fn=>{timers.set(++id,fn);return id},clearTimeout:n=>timers.delete(n),ResizeObserver:class{observe(){}},input:cluster};
 vm.createContext(sandbox);vm.runInContext(script.replace(/boot\(\);\s*$/,''),sandbox);
 const run=code=>vm.runInContext(code,sandbox);
 run('data=input;metadata();if(typeof prepareGeometry!=="undefined")prepareGeometry();draw()');
 return {run,node,calls,frames,timers,flush(){const pending=[...frames.values()];frames.clear();for(const fn of pending)fn()},reset(){for(const name in calls)delete calls[name]}};
}
function test(script,cluster){
 const r=renderer(script,cluster),canvas=r.node('canvas'),n=cluster.elements.length,original=JSON.stringify(cluster);
 const heavy=cluster.elements.filter(e=>e!=='H').length;
 assert.equal(+canvas.dataset.visibleAtoms,heavy);
 r.run('hydrogen=true;prepareVisibility();draw()');assert.equal(+canvas.dataset.drawnAtoms,n);
 assert.equal(canvas.dataset.quality,'full');
 const before=JSON.parse(r.run('JSON.stringify(Array.from(points.slice(0,3)))'));
 r.run('ry+=.4;draw()');const after=JSON.parse(r.run('JSON.stringify(Array.from(points.slice(0,3)))'));
 if(Math.abs(cluster.coordinates[0][0])+Math.abs(cluster.coordinates[0][2])>1e-6)assert.notEqual(before[0],after[0]);
 // Independently calculate one projected coordinate after changing the view.
 const view=JSON.parse(r.run('JSON.stringify({rx,ry,radius,zoom})')),p=cluster.coordinates[n-1],scale=700*.44/view.radius*view.zoom;
 const x=p[0]*Math.cos(view.ry)+p[2]*Math.sin(view.ry),z=-p[0]*Math.sin(view.ry)+p[2]*Math.cos(view.ry);
 assert.ok(Math.abs(r.run(`points[${(n-1)*3}]`)-(500+x*scale))<1e-3);
 assert.ok(Math.abs(r.run(`points[${(n-1)*3+1}]`)-(350-(p[1]*Math.cos(view.rx)-z*Math.sin(view.rx))*scale))<1e-3);
 r.reset();r.node('canvas').onpointerdown({clientX:10,clientY:10,pointerId:1});r.flush();
 assert.equal(canvas.dataset.quality,n>4000?'preview':'full');
 assert.ok(+canvas.dataset.drawnAtoms<=Math.min(n,cluster.central_atoms+12000));
 assert.equal(+canvas.dataset.visibleAtoms,n); // displayed count describes full selected scene
 assert.ok(r.run('render.groups[0].atoms.every(i=>render.quick.includes(i))'));
 r.node('canvas').onpointermove({clientX:20,clientY:30});r.flush();
 r.node('canvas').onpointercancel();r.flush();assert.equal(canvas.dataset.quality,'full');assert.equal(+canvas.dataset.drawnAtoms,n);
 r.node('canvas').wheel({preventDefault(){},deltaY:10});r.flush();assert.equal(r.timers.size,1);
 for(const timer of [...r.timers.values()])timer();r.timers.clear();r.flush();assert.equal(canvas.dataset.quality,'full');
 r.run('mode="near";prepareVisibility();draw()');
 const near=cluster.atom_groups.filter(g=>g<=4).length;assert.equal(+canvas.dataset.visibleAtoms,near);
 assert.equal(r.node('export').href,`/api/cluster.xyz?center=${cluster.center}&cutoff=${cluster.cutoff}`);
 assert.equal(JSON.stringify(cluster),original); // preview/toggles must never mutate source/export data
 r.node('canvas').onpointermove({clientX:1,clientY:1});r.node('canvas').onpointermove({clientX:2,clientY:2});assert.equal(r.frames.size,1);r.flush();
 // Hydrogen-only empty display, then restoration, uses no stale selection/bounds.
 const h={...cluster,coordinates:[[1,0,0]],elements:['H'],bonds:[],atom_groups:[0],groups:[{molecule:1,label:'H',image:[0,0,0],start:0,count:1}],names:['H'],serials:[1],global_indices:[1],central_atoms:1,neighbor_atoms:0,neighbor_images:0,periodic_images:0,minimum_distance:null};
 const empty=renderer(script,h);assert.equal(+empty.node('canvas').dataset.drawnAtoms,0);empty.run('hydrogen=true;prepareVisibility();draw()');assert.equal(+empty.node('canvas').dataset.drawnAtoms,1);
 return {atoms:n,full_scene_unchanged:true,preview_and_cancel_restore:true,projection_matches:true};
}
function benchmark(script,cluster,preview=false){
 const r=renderer(script,cluster);r.run(`hydrogen=true;if(typeof prepareVisibility!=="undefined")prepareVisibility();${preview?'motion=true;':''}`);
 for(let i=0;i<5;i++)r.run('rx+=.002;draw()');r.reset();
 const samples=[];for(let i=0;i<25;i++){const start=performance.now();r.run('rx+=.002;draw()');samples.push(performance.now()-start)}samples.sort((a,b)=>a-b);
 return {js_only_median_ms:samples[12],js_only_p95_ms:samples[23],stroke_calls_per_frame:(r.calls.stroke||0)/25,fill_calls_per_frame:(r.calls.fill||0)/25,visible_atoms:+r.node('canvas').dataset.visibleAtoms,drawn_atoms:+(r.node('canvas').dataset.drawnAtoms||r.node('canvas').dataset.visibleAtoms)};
}
if(require.main===module){
 const input=JSON.parse(require('node:fs').readFileSync(0,'utf8'));
 const result={checks:test(input.script,input.cluster)};
 if(input.baseline_script)result.benchmark={baseline:benchmark(input.baseline_script,input.cluster),optimized_full:benchmark(input.script,input.cluster),optimized_preview:benchmark(input.script,input.cluster,true)};
 console.log(JSON.stringify(result));
}
module.exports={renderer,test,benchmark};
