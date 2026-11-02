"""Fetch each column neuron's public reconstruction, linked by stable nucleus ID."""

import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import shutil
import time
import urllib.request

import numpy as np
from analyze import load_swc

META_BASE='https://storage.googleapis.com/mat_dbs/public/minnie65_phase3_v1/v661'
META_TABLE='allen_v1_column_types_slanted_ref_merged'
SKEL_BASE='https://bossdb-open-data.s3.amazonaws.com/iarpa_microns/minnie/minnie65/skeletons/v661/skeletons'


def repair_zero_radii(rows, parent, order):
    """Fill zero radii from the closest measured upstream radius; keep positives."""
    radii=rows[:,5].copy();repairs=[]
    if np.any(radii<0):raise ValueError('Negative reconstructed radius')
    root=int(order[0])
    if radii[root]<=0:raise ValueError('Soma radius is missing; no invented soma is permitted')
    for i in order:
        if radii[i]==0:
            ancestor=int(parent[i])
            while rows[ancestor,5]==0:ancestor=int(parent[ancestor])
            radii[i]=radii[ancestor]
            repairs.append({'node_id':int(rows[i,0]),'original_radius_um':0.,'radius_um':float(radii[i]),'measured_ancestor_node_id':int(rows[ancestor,0])})
    return radii,repairs


def collapse_identical_vertices(rows, parent, order):
    """Remove only redundant parent/child aliases at the identical cable point."""
    representatives=np.arange(len(rows));collapsed=[]
    for i in order:
        p=int(parent[i])
        if p>=0 and np.array_equal(rows[i,2:5],rows[p,2:5]):
            if rows[i,1]!=rows[p,1] or rows[i,5]!=rows[p,5]:
                raise ValueError('Coincident vertices have different compartment labels or radii')
            representatives[i]=representatives[p]
            collapsed.append({'node_id':int(rows[i,0]),'equivalent_node_id':int(rows[representatives[i],0])})
    keep=np.flatnonzero(representatives==np.arange(len(rows)))
    cleaned=rows[keep].copy()
    old_to_new={int(old):i for i,old in enumerate(keep)}
    for j,i in enumerate(keep):
        if parent[i]>=0:cleaned[j,6]=rows[representatives[parent[i]],0]
    ids={int(row[0]):i for i,row in enumerate(cleaned)}
    new_parent=np.array([-1 if row[6]==-1 else ids[int(row[6])] for row in cleaned])
    new_order=np.array([old_to_new[int(i)] for i in order if int(i) in old_to_new])
    return cleaned,new_parent,new_order,collapsed


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--workers',type=int,default=8)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    raw=args.output/'raw';ready=args.output/'swc';raw.mkdir(exist_ok=True);ready.mkdir(exist_ok=True)
    with (args.data/'nodes.csv').open() as stream:nodes=list(csv.DictReader(stream))
    sources={}
    for suffix in ['_header.csv','.csv.gz']:
        url=META_BASE+'/'+META_TABLE+suffix
        with urllib.request.urlopen(url,timeout=60) as stream:content=stream.read()
        (args.output/(META_TABLE+suffix)).write_bytes(content)
        sources[suffix]={'url':url,'sha256':hashlib.sha256(content).hexdigest()}
    columns=[r[0] for r in csv.reader((args.output/(META_TABLE+'_header.csv')).open())]
    with gzip.open(args.output/(META_TABLE+'.csv.gz'),'rt') as stream:
        old={r['target_id']:r for r in csv.DictReader(stream,fieldnames=columns)}
    def fetch(node):
        root,nucleus=node['root_id'],node['nucleus_id'];reference=old[nucleus]
        filename=f"{reference['pt_root_id']}_{nucleus}.swc";url=SKEL_BASE+'/'+filename;path=raw/filename
        if not path.exists():
            for attempt in range(3):
                try:
                    with urllib.request.urlopen(url,timeout=60) as stream:content=stream.read()
                    temporary=path.with_suffix('.part');temporary.write_bytes(content);temporary.replace(path);break
                except Exception:
                    if attempt==2:raise
                    time.sleep(1)
        content=path.read_bytes();rows,parent,order=load_swc(path,require_positive_radii=False)
        if not set(rows[:,1]).issubset({1,2,3,4}):raise ValueError('Unexpected reconstruction compartment label')
        root_row=rows[int(order[0])]
        if root_row[1]!=1:raise ValueError('Reconstruction is not rooted at its soma')
        soma=np.array([float(node[a+'_um']) for a in 'xyz'])
        distance=float(np.linalg.norm(root_row[2:5]-soma))
        if distance>.01:raise ValueError(f'Soma/nucleus coordinate mismatch: {distance:g} um')
        original=rows.copy();original_parent=parent.copy()
        rows,parent,order,collapsed=collapse_identical_vertices(rows,parent,order)
        old_length=sum(np.linalg.norm(original[i,2:5]-original[p,2:5]) for i,p in enumerate(original_parent) if p>=0)
        new_length=sum(np.linalg.norm(rows[i,2:5]-rows[p,2:5]) for i,p in enumerate(parent) if p>=0)
        if not np.isclose(old_length,new_length,atol=1e-6,rtol=1e-12):raise RuntimeError('Cable length changed')
        radii,repairs=repair_zero_radii(rows,parent,order)
        target=ready/(root+'.swc')
        if repairs or collapsed:
            values={int(rows[i,0]):repr(float(radii[i])) for i in np.flatnonzero(rows[:,5]==0)}
            parents={int(row[0]):str(int(row[6])) for row in rows}
            lines=[]
            for line in content.decode().splitlines():
                fields=line.split()
                if fields and not line.lstrip().startswith('#'):
                    node_id=int(fields[0])
                    if node_id not in parents:continue
                    if node_id in values:fields[5]=values[node_id]
                    fields[6]=parents[node_id];line=' '.join(fields)
                lines.append(line)
            target.write_text('\n'.join(lines)+'\n')
        else:shutil.copyfile(path,target)
        imported,_,_=load_swc(target)
        if not np.array_equal(imported[:,[0,1,2,3,4,6]],rows[:,[0,1,2,3,4,6]]):raise RuntimeError('Branch coordinates or topology changed')
        return {'root_id':root,'nucleus_id':nucleus,'reconstruction_root_id':reference['pt_root_id'],
                'source_url':url,'raw_sha256':hashlib.sha256(content).hexdigest(),'import_sha256':hashlib.sha256(target.read_bytes()).hexdigest(),
                'manual_cell_type':node['manual_cell_type'],'raw_nodes':len(original),'nodes':len(rows),'types':dict(Counter(int(x) for x in rows[:,1])),
                'soma_nucleus_distance_um':distance,'has_axon':bool(np.any(rows[:,1]==2)),'has_apical_labels':bool(np.any(rows[:,1]==4)),
                'minimum_positive_raw_radius_um':float(rows[rows[:,5]>0,5].min()),'repairs':repairs,'collapsed_identical_vertices':collapsed,'cable_length_um':float(new_length)}
    records=[];failures=[]
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        pending={pool.submit(fetch,node):node for node in nodes}
        for future in as_completed(pending):
            node=pending[future]
            try:records.append(future.result())
            except Exception as error:failures.append({'root_id':node['root_id'],'nucleus_id':node['nucleus_id'],'error':str(error)})
            if (len(records)+len(failures))%100==0:print(f'Reconstructions checked {len(records)+len(failures)}/{len(nodes)}; failures {len(failures)}',flush=True)
    records.sort(key=lambda row:row['root_id'])
    summary={'reconstruction_version':661,'wiring_version':943,'link_method':'Stable nucleus ID plus soma position agreement within 0.01 um',
             'metadata_sources':sources,'requested':len(nodes),'ready':len(records),'failures':failures,'nodes':sum(r['nodes'] for r in records),
             'with_recorded_axon':sum(r['has_axon'] for r in records),'with_apical_labels':sum(r['has_apical_labels'] for r in records),
             'soma_only_cells':sum(r['nodes']==1 for r in records),'cells_with_radius_repairs':sum(bool(r['repairs']) for r in records),
             'zero_radius_repairs':sum(len(r['repairs']) for r in records),'records':records,
             'collapsed_identical_vertices':sum(len(r['collapsed_identical_vertices']) for r in records),
             'repair_policy':'Only zero radii inherit a positive upstream radius. Parent/child duplicates with identical coordinates, label and radius are merged without changing cable length or branch layout. All retained coordinates and positive radii are unchanged. No branches or axons are added.'}
    (args.output/'manifest.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps({k:summary[k] for k in ['requested','ready','failures','nodes','with_recorded_axon','with_apical_labels','soma_only_cells','cells_with_radius_repairs','zero_radius_repairs']},indent=2),flush=True)
    if failures:raise SystemExit('Some reconstructions failed validation; no silent template fallback is allowed')


if __name__=='__main__':main()
