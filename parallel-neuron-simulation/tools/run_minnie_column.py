"""Simulate Minnie wiring with recorded skeletons or published M1 shape surrogates."""

import argparse
from collections import Counter
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import runpy
import sys
import time
import platform

from minnie_network import PROFILES, APPROXIMATE_TYPES, receptor_components, graph_summary
from run_mouse_m1 import non_thread_safe_mechanisms, rss_bytes

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'vendor/m1_260015'


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--threads',type=int,default=32)
    parser.add_argument('--duration',type=float,default=100.)
    parser.add_argument('--weight-e',type=float,default=.001,help='Assumed base conductance per excitatory contact, uS')
    parser.add_argument('--weight-i',type=float,default=.001,help='Assumed base conductance per inhibitory contact, uS')
    parser.add_argument('--weight-factor',type=float,default=1.,help='0 creates a matched zero-transmission control')
    parser.add_argument('--stimulate-cells',type=int,default=32)
    parser.add_argument('--stimulus-profile',default='IT5B_reduced',choices=sorted({item[0] for item in PROFILES.values()}))
    parser.add_argument('--amplitude',type=float,default=.5,help='Controlled soma input, nA')
    parser.add_argument('--morphologies',type=Path,help='Prepared per-neuron skeleton directory, with manifest.json and swc/')
    parser.add_argument('--electrical-profiles',type=Path,help='Exported M1 density priors; defaults to the morphology directory')
    parser.add_argument('--d-lambda',type=float,default=.1,help='Maximum segment electrotonic fraction at the specified frequency')
    parser.add_argument('--frequency',type=float,default=100.,help='Spatial discretization frequency, Hz')
    parser.add_argument('--paired-control',action='store_true',help='Reinitialize the same constructed cells with zero transmission after the coupled run')
    args=parser.parse_args()
    if args.threads<1 or not 0<args.duration<500 or min(args.weight_e,args.weight_i,args.weight_factor,args.stimulate_cells)<0 or args.d_lambda<=0 or args.frequency<=0:
        parser.error('Invalid threads, duration or conductance/stimulation parameters')
    if args.paired_control and args.weight_factor<=0:parser.error('Paired control requires a nonzero initial weight factor')
    data=args.data.resolve();output=args.output.resolve();output.mkdir(parents=True,exist_ok=True)
    morphology_dir=args.morphologies.resolve() if args.morphologies else None
    electrical_path=args.electrical_profiles.resolve() if args.electrical_profiles else morphology_dir/'electrical_profiles.json' if morphology_dir else None
    os.environ['NEURON_MODULE_OPTIONS']=os.environ.get('NEURON_MODULE_OPTIONS','')+' -nogui'
    os.environ['MPLBACKEND']='Agg';os.environ['MPLCONFIGDIR']=str(output/'matplotlib-cache')
    os.environ['OMP_NUM_THREADS']='1';os.environ['OPENBLAS_NUM_THREADS']='1'
    import numpy as np
    import neuron
    from neuron import h
    import netpyne
    from netpyne import sim

    started=time.perf_counter()
    with (data/'nodes.csv').open() as stream:nodes=list(csv.DictReader(stream))
    with (data/'edges.csv').open() as stream:edges=list(csv.DictReader(stream))
    lookup={row['root_id']:row for row in nodes}
    if len(lookup)!=len(nodes):raise RuntimeError('Duplicate neuron identities')
    graph=graph_summary(nodes,edges)
    (output/'graph.json').write_text(json.dumps(graph,indent=2)+'\n')
    manifest=json.loads((SOURCE/'SOURCE_MANIFEST.json').read_text())
    for item in manifest['files']:
        if hashlib.sha256((SOURCE/item['path']).read_bytes()).hexdigest()!=item['sha256']:
            raise RuntimeError('Published M1 source changed')
    source_dir=SOURCE/'sim';neuron.load_mechanisms(str(source_dir));os.chdir(source_dir);sys.path.insert(0,str(source_dir))
    global cfg
    cfg=runpy.run_path(str(source_dir/'cfg.py'))['cfg']
    cfg.scale=1.;cfg.duration=args.duration
    cfg.addConn=0;cfg.addLongConn=0;cfg.addSubConn=0;cfg.addIClamp=0
    cfg.recordLFP=[];cfg.recordDipole=False;cfg.recordDipoles=False
    cfg.recordTraces={};cfg.recordCells=[];cfg.recordTime=False;cfg.analysis={}
    cfg.printRunTime=False;cfg.verbose=False;cfg.savePickle=False;cfg.saveJson=False
    parameters=runpy.run_path(str(source_dir/'netParams.py'))['netParams']
    parameters.popParams.clear();parameters.connParams.clear();parameters.subConnParams.clear()
    parameters.stimSourceParams.clear();parameters.stimTargetParams.clear()
    groups={}
    for node in nodes:
        label,cell_type,cell_model=PROFILES[node['manual_cell_type']]
        groups.setdefault(label,[]).append(node)
    for label in list(parameters.cellParams):
        if label not in groups:del parameters.cellParams[label]
        else:parameters.cellParams[label]['conds']={'pop':label}
    for label,members in sorted(groups.items()):
        _,cell_type,cell_model=PROFILES[members[0]['manual_cell_type']]
        parameters.popParams[label]={
            'cellType':cell_type,'cellModel':cell_model,
            'cellsList':[{'root_id':n['root_id'],'manual_cell_type':n['manual_cell_type'],'ei_class':n['ei_class'],
                          **{a:float(n[a+'_um']) for a in 'xyz'}} for n in sorted(members,key=lambda n:n['root_id'])],
        }
    sim.initialize(simConfig=cfg,netParams=parameters);sim.pc.nthread(1)
    sim.net.createPops()
    geometry_records={}
    if morphology_dir:
        from reconstructed_cell import ReconstructedCell, ElectricalProfile
        h.load_file('stdrun.hoc');h.load_file('import3d.hoc')
        geometry_manifest=json.loads((morphology_dir/'manifest.json').read_text())
        geometry_records={row['root_id']:row for row in geometry_manifest['records']}
        electrical_data=json.loads(electrical_path.read_text())
        if electrical_data['source_commit']!=manifest['commit']:raise RuntimeError('Electrical reference source mismatch')
        density_profiles={label:ElectricalProfile(reference) for label,reference in electrical_data['profiles'].items()}
        sim.net.cells=[]
        for label,members in sorted(groups.items()):
            _,kind,model=PROFILES[members[0]['manual_cell_type']]
            for node in sorted(members,key=lambda n:n['root_id']):
                path=morphology_dir/'swc'/(node['root_id']+'.swc')
                record=geometry_records[node['root_id']]
                if record['nucleus_id']!=node['nucleus_id'] or hashlib.sha256(path.read_bytes()).hexdigest()!=record['import_sha256']:
                    raise RuntimeError('Reconstruction identity or checksum mismatch')
                tagged={**node,'profile':label,'cell_type':kind,'cell_model':model}
                cell=ReconstructedCell(path,tagged,len(sim.net.cells),density_profiles[label],sim.pc,args.d_lambda,args.frequency)
                sim.net.cells.append(cell)
                if len(sim.net.cells)%25==0:print(f'Imported recorded branches {len(sim.net.cells)}/{len(nodes)}; elapsed {time.perf_counter()-started:.1f} s',flush=True)
        geometry_info={'mode':'individual_recorded_skeletons','reconstruction_version':geometry_manifest['reconstruction_version'],
                       'wiring_version':943,'manifest_sha256':hashlib.sha256((morphology_dir/'manifest.json').read_bytes()).hexdigest(),
                       'electrical_profile_sha256':hashlib.sha256(electrical_path.read_bytes()).hexdigest(),
                       'source_skeleton_points':sum(geometry_records[node['root_id']]['nodes'] for node in nodes),
                       'with_recorded_axon':sum(geometry_records[node['root_id']]['has_axon'] for node in nodes),
                       'd_lambda':args.d_lambda,'discretization_frequency_hz':args.frequency,
                       'density_transfer':'Compartment class and physical soma path distance, with endpoint clamping and median over coincident reference distances',
                       'spike_source':'Somatic upward crossing of 0 mV for every reconstructed neuron'}
    else:
        sim.net.createCells()
        geometry_info={'mode':'surrogate_M1_templates'}
    cells={cell.tags['root_id']:cell for cell in sim.net.cells}
    if set(cells)!=set(lookup):raise RuntimeError('Model cells and observed neuron IDs disagree')
    print(f'Created {len(cells)} typed Minnie cells in {time.perf_counter()-started:.2f} s',flush=True)
    rng=np.random.default_rng(20261027)
    synapses=[];connections=[];anchors={};placement_cache={};receptor_counts=Counter()
    delays=[];base_weights=[]
    def placement_options(root,kind):
        key=(root,kind)
        if key not in placement_cache:
            cell=cells[root]
            names=list(cell.secLists.get(kind,[]))
            names=[name for name in names if name in cell.secs] or ['soma']
            weights=np.array([max(cell.section_areas[name] if morphology_dir else cell.secs[name]['hObj'].L*cell.secs[name]['hObj'](.5).diam,1e-12) for name in names])
            placement_cache[key]=(names,weights/weights.sum())
        return placement_cache[key]
    with gzip.open(output/'modeled_sites.csv.gz','wt',compresslevel=1,newline='') as stream:
        writer=csv.writer(stream);writer.writerow(['pre_root_id','post_root_id','contact_ordinal','model_section','model_x','base_weight_us','target_normalization','delay_ms'])
        for index,edge in enumerate(edges,1):
            pre,post=edge['pre_pt_root_id'],edge['post_pt_root_id'];pre_node,post_node=lookup[pre],lookup[post]
            components=receptor_components(pre_node,post_node)
            if pre_node['ei_class']=='E':kind='spiny'
            elif pre_node['manual_cell_type']=='BC' or post_node['ei_class']=='I':kind='perisom'
            else:kind='alldend'
            names,probabilities=placement_options(post,kind)
            distance=sum((float(pre_node[a+'_um'])-float(post_node[a+'_um']))**2 for a in 'xyz')**.5
            delay=2.+distance/500.
            base_weight=args.weight_e if pre_node['ei_class']=='E' else args.weight_i
            for ordinal in range(int(edge['n_syn'])):
                section_name=names[int(rng.choice(len(names),p=probabilities))]
                section_info=cells[post].secs[section_name];section=section_info['hObj'];x=float(rng.uniform(.1,.9))
                # CLI weights are already conductances in uS. The published
                # weightNorm converts nominal somatic EPSPs (mV) to uS and
                # must not be applied again to these physical conductances.
                norm=1.
                for label,fraction in components:
                    specification=parameters.synMechParams[label];synapse=getattr(h,specification['mod'])(section(x))
                    for name,value in specification.items():
                        if name!='mod':setattr(synapse,name,value)
                    connection=sim.pc.gid_connect(cells[pre].gid,synapse)
                    connection.delay=delay;connection.weight[0]=base_weight*fraction*norm*args.weight_factor
                    connections.append(connection);synapses.append(synapse);anchors.setdefault(post,synapse);receptor_counts[label]+=1
                writer.writerow([pre,post,ordinal,section_name,x,base_weight,norm,delay])
            delays.append(delay);base_weights.append(base_weight)
            if index%10000==0:print(f'Connected {index:,}/{len(edges):,} observed pairs',flush=True)
    hoc_count=sum(int(h.List(name).count()) for name in ['MyExp2SynBB','MyExp2SynNMDABB'])
    if hoc_count!=len(synapses):raise RuntimeError('Synaptic object count mismatch')
    sections=[s['hObj'] for cell in cells.values() for s in cell.secs.values()]
    unsafe=non_thread_safe_mechanisms(source_dir);partitions=[h.SectionList() for _ in range(args.threads)];loads=[0]*args.threads;pinned=[]
    for root,cell in sorted(cells.items(),key=lambda item:(-sum(s['hObj'].nseg for s in item[1].secs.values()),item[0])):
        restricted={m for s in cell.secs.values() for m in s.get('mechs',{})}&unsafe
        thread=0 if restricted else min(range(args.threads),key=lambda i:loads[i])
        if restricted:pinned.append(root)
        for info in cell.secs.values():
            if info['hObj'].parentseg() is None:partitions[thread].append(sec=info['hObj'])
        loads[thread]+=sum(s['hObj'].nseg for s in cell.secs.values())
    sim.pc.nthread(args.threads)
    for thread,roots in enumerate(partitions):sim.pc.partition(thread,roots)
    h('objref minnie_thread_context');h.minnie_thread_context=sim.pc
    h('func minnie_thread_id() { return minnie_thread_context.sec_in_thread() }')
    stimulated=sorted((n['root_id'] for n in nodes if n['ei_class']=='E' and PROFILES[n['manual_cell_type']][0]==args.stimulus_profile),key=lambda root:(-graph['outgoing_contacts'].get(root,0),root))[:args.stimulate_cells]
    clamps=[]
    for root in stimulated:
        clamp=h.IClamp(cells[root].secs['soma']['hObj'](.5));clamp.amp=args.amplitude;clamp.delay=10.;clamp.dur=50.;clamps.append(clamp)
    ordered=sorted(cells)
    records=[]
    for root in ordered:
        section=cells[root].secs['soma']['hObj']
        if root not in anchors:
            anchor=h.IClamp(section(.5));anchor.amp=anchor.dur=0.;anchors[root]=anchor
        records.append(h.Vector().record(anchors[root],section(.5)._ref_v))
    spike_times,spike_gids=h.Vector(),h.Vector();sim.pc.spike_record(-1,spike_times,spike_gids)
    with (output/'cell_models.csv').open('w',newline='') as stream:
        writer=csv.writer(stream);writer.writerow(['root_id','gid','manual_type','ei_class','M1_profile','coarse_subtype_surrogate','stimulated','geometry_mode','skeleton_points','recorded_axon'])
        for root in ordered:writer.writerow([root,cells[root].gid,lookup[root]['manual_cell_type'],lookup[root]['ei_class'],PROFILES[lookup[root]['manual_cell_type']][0],lookup[root]['manual_cell_type'] in APPROXIMATE_TYPES,root in stimulated,geometry_info['mode'],geometry_records.get(root,{}).get('nodes'),geometry_records.get(root,{}).get('has_axon')])
    summary={
        'scope':'Observed Minnie wiring with individual recorded skeletons and transferred electrical priors' if morphology_dir else 'Observed Minnie v943 within-column wiring with approximate M1 cell templates; not a calibrated digital twin',
        'geometry':geometry_info,
        'hostname':platform.node(),'neuron_version':neuron.__version__,'netpyne_version':netpyne.__version__,
        'anatomical_source_version':943,'M1_source_commit':manifest['commit'],
        'input_hashes':{name:hashlib.sha256((data/name).read_bytes()).hexdigest() for name in ['nodes.csv','edges.csv','summary.json']},
        'neurons':len(cells),'ei_classes':dict(Counter(n['ei_class'] for n in nodes)),
        'anatomical_contact_count':graph['anatomical_contacts'],'directed_pairs':len(edges),
        'synaptic_receptor_objects':len(synapses),'receptor_counts':dict(receptor_counts),
        'sections':len(sections),'segments':sum(s.nseg for s in sections),
        'threads':args.threads,'mpi_ranks':int(sim.pc.nhost()),'cells_restricted_to_thread_zero':len(pinned),
        'duration_ms':args.duration,'dt_ms':float(cfg.dt),'temperature_c':34,
        'base_weight_e_us':args.weight_e,'base_weight_i_us':args.weight_i,'weight_factor':args.weight_factor,
        'synaptic_weight_convention':'Raw per-contact conductance in uS; no EPSP-to-conductance normalization',
        'delay_ms_range':[min(delays),max(delays)],'placement_seed':20261027,
        'stimulus':{'cells':stimulated,'electrical_profile':args.stimulus_profile,'amplitude_na':args.amplitude,'onset_ms':10.,'duration_ms':50.},
        'assumptions':['Individual recorded v661 geometries are linked to v943 wiring by stable nucleus ID; electrical-density profiles are transferred from M1, not measured Minnie conductances.' if morphology_dir else 'Published M1 electrical cell templates represent Minnie V1 neurons; actual Minnie skeletons are not used as model geometry.',
                       'Morphological classes are mapped to a smaller family of electrical templates; some subtype mappings are coarse surrogates.',
                       'Synaptic kinetics and receptor mixtures follow M1; per-contact conductances are assumed in uS and are not EPSP-normalized.',
                       'Anatomical contact counts are observed, but locations are sampled on model compartments, not measured synapse coordinates.',
                       'Delays use 2 ms plus soma separation / 500 um/ms; delays are assumptions.',
                       'The first experiment uses controlled soma currents; activity from outside the selected column is omitted.'],
        'construction_seconds':time.perf_counter()-started,
    }
    print(json.dumps({k:summary[k] for k in ['neurons','anatomical_contact_count','synaptic_receptor_objects','segments','construction_seconds']},indent=2),flush=True)
    pre_run_done=False
    def simulate(destination,factor,reused):
        nonlocal pre_run_done
        destination.mkdir(parents=True,exist_ok=True)
        result=dict(summary);result['weight_factor']=factor
        if reused:
            for connection in connections:connection.weight[0]=0.
            result['construction_seconds']=0.;result['shared_construction_from']=str(output)
        spike_times.resize(0);spike_gids.resize(0)
        init=time.perf_counter()
        if not pre_run_done:sim.preRun();pre_run_done=True
        h.finitialize(float(cfg.hParams['v_init']))
        result['initialization_seconds']=time.perf_counter()-init
        actual=[0]*args.threads
        for section in sections:actual[int(h.minnie_thread_id(sec=section))]+=section.nseg
        if actual!=loads:raise RuntimeError('Thread partition verification failed')
        result['segments_per_thread']=actual
        solve=time.perf_counter();stop=0.
        while stop<args.duration:
            stop=min(args.duration,stop+10.);sim.pc.psolve(stop)
            print(f'{destination.name}: {stop:g}/{args.duration:g} ms; solve {time.perf_counter()-solve:.2f} s',flush=True)
        result['solve_seconds']=time.perf_counter()-solve;result['actual_final_time_ms']=float(h.t)
        voltages=np.asarray([np.asarray(v).copy()[::4] for v in records]);times=np.arange(voltages.shape[1])*float(cfg.dt)*4
        final=np.fromiter((segment.v for section in sections for segment in section),dtype=float)
        if not np.isfinite(voltages).all() or not np.isfinite(final).all() or abs(float(h.t)-args.duration)>1e-7:raise RuntimeError('Invalid simulation output')
        event_times=np.asarray(spike_times).copy();event_gids=np.asarray(spike_gids,dtype=np.int64).copy();event_order=np.lexsort((event_gids,event_times))
        event_times,event_gids=event_times[event_order],event_gids[event_order]
        gid_to_root={cell.gid:root for root,cell in cells.items()};spike_roots=np.asarray([gid_to_root[int(g)] for g in event_gids],dtype='U18')
        result['recorded_neuronal_spikes']=len(event_times);result['cells_with_spikes']=len(set(spike_roots.tolist()))
        result['stimulated_cell_spikes']=sum(root in set(stimulated) for root in spike_roots)
        result['nonstimulated_cell_spikes']=result['recorded_neuronal_spikes']-result['stimulated_cell_spikes']
        result['spikes_by_ei_class']=dict(Counter(lookup[root]['ei_class'] for root in spike_roots))
        result['resident_ram_bytes_after_simulation']=rss_bytes();result['all_recorded_voltages_and_final_compartment_voltages_finite']=True
        np.savez_compressed(destination/'recordings.npz',time_ms=times,soma_mv=voltages,root_ids=np.asarray(ordered),spike_times_ms=event_times,spike_root_ids=spike_roots)
        (destination/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
        print(json.dumps({k:result[k] for k in ['threads','weight_factor','solve_seconds','recorded_neuronal_spikes','nonstimulated_cell_spikes','resident_ram_bytes_after_simulation']},indent=2),flush=True)
    simulate(output,args.weight_factor,False)
    if args.paired_control:simulate(output/'zero_transmission_control',0.,True)


if __name__=='__main__':main()
