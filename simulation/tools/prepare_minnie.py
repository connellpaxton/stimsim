"""Extract a version-consistent, cell-typed Minnie column connectivity table."""

import argparse
from collections import Counter
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import time
import urllib.request

BASE = 'https://storage.googleapis.com/mat_dbs/public/minnie65_phase3_v1/v943'
TABLE = 'allen_v1_column_types_slanted_ref_merged'
EDGE_COLUMNS = ['pre_pt_root_id','post_pt_root_id','n_syn','sum_size','pre_nuc_id','post_nuc_id']


def download(url, path):
    """Cache immutable public input and retain its hash and source URL."""
    checksum = path.with_suffix(path.suffix+'.sha256')
    if path.exists() and checksum.exists():
        return {'url':url, 'bytes':path.stat().st_size, 'sha256':checksum.read_text().strip()}
    partial = path.with_suffix(path.suffix+'.part')
    digest = hashlib.sha256()
    with urllib.request.urlopen(url, timeout=120) as source, partial.open('wb') as target:
        size = int(source.headers.get('Content-Length',0))
        copied, next_report = 0, 0
        while chunk := source.read(8*1024*1024):
            target.write(chunk)
            digest.update(chunk)
            copied += len(chunk)
            if copied >= next_report:
                print(f'Downloaded {path.name}: {copied/2**20:.0f}/{size/2**20:.0f} MiB',flush=True)
                next_report = copied+256*1024*1024
    if size and copied != size:
        raise RuntimeError('Incomplete public data download')
    partial.replace(path)
    checksum.write_text(digest.hexdigest()+'\n')
    return {'url':url,'bytes':copied,'sha256':digest.hexdigest()}


def column_nodes(rows):
    """Retain uniquely identified neurons with a measured E/I classification."""
    multiplicity = Counter(row['pt_root_id'] for row in rows)
    selected, excluded = [], []
    for row in rows:
        system = row['classification_system']
        if system not in {'aibs_coarse_excitatory','aibs_coarse_inhibitory'}:
            excluded.append({'nucleus_id':row['target_id'],'reason':'No unambiguous excitatory/inhibitory classification'})
            continue
        if multiplicity[row['pt_root_id']] != 1 or int(row['pt_root_id']) <= 0:
            excluded.append({'nucleus_id':row['target_id'],'reason':'Ambiguous or invalid segmentation root ID'})
            continue
        selected.append({
            'root_id':row['pt_root_id'],'nucleus_id':row['target_id'],
            'manual_cell_type':row['cell_type'],
            'ei_class':'E' if system=='aibs_coarse_excitatory' else 'I',
            'x_um':float(row['pt_position_x'])*.004,
            'y_um':float(row['pt_position_y'])*.004,
            'z_um':float(row['pt_position_z'])*.04,
        })
    return selected, excluded


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    args.cache.mkdir(parents=True,exist_ok=True)
    args.output.mkdir(parents=True,exist_ok=True)
    sources={}
    for name in [TABLE+'_header.csv',TABLE+'.csv.gz','connections_with_nuclei_header.csv','connections_with_nuclei.csv.gz']:
        sources[name]=download(BASE+'/'+name,args.cache/name)
    columns=[r[0] for r in csv.reader((args.cache/(TABLE+'_header.csv')).open())]
    with gzip.open(args.cache/(TABLE+'.csv.gz'),'rt') as stream:
        rows=list(csv.DictReader(stream,fieldnames=columns))
    nodes,excluded=column_nodes(rows)
    with (args.output/'nodes.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(nodes[0]));writer.writeheader();writer.writerows(nodes)
    edge_header=[r[0] for r in csv.reader((args.cache/'connections_with_nuclei_header.csv').open())]
    if edge_header != EDGE_COLUMNS:
        raise RuntimeError('Unexpected archived connectivity schema')
    import pandas as pd
    ids={int(row['root_id']) for row in nodes}
    totals=Counter();seen=0;started=time.perf_counter();edge_count=synapses=0
    dtypes={name:('float64' if name=='sum_size' else 'int64') for name in EDGE_COLUMNS}
    with (args.output/'edges.csv').open('w',newline='') as stream:
        for chunk in pd.read_csv(args.cache/'connections_with_nuclei.csv.gz',header=None,names=EDGE_COLUMNS,dtype=dtypes,chunksize=2_000_000):
            selected=chunk.loc[chunk.pre_pt_root_id.isin(ids) & chunk.post_pt_root_id.isin(ids)]
            if len(selected):
                if (selected.n_syn<=0).any():raise RuntimeError('Nonpositive anatomical synapse count')
                selected.to_csv(stream,index=False,header=edge_count==0)
                edge_count+=len(selected);synapses+=int(selected.n_syn.sum())
            seen+=len(chunk)
            print(f'Scanned {seen:,} pairs; retained {edge_count:,} pairs / {synapses:,} anatomical contacts; {time.perf_counter()-started:.1f} s',flush=True)
    summary={
        'source':'MICrONS Minnie65, archived materialization v943',
        'source_documentation':'https://tutorial.microns-explorer.org/static-repositories.html',
        'sources':sources,'column_annotation_rows':len(rows),'selected_neurons':len(nodes),
        'excluded_annotations':excluded,'cell_types':dict(Counter(row['manual_cell_type'] for row in nodes)),
        'ei_classes':dict(Counter(row['ei_class'] for row in nodes)),
        'scanned_cell_pairs':seen,'internal_directed_pairs':edge_count,'internal_anatomical_synapses':synapses,
        'interpretation':'Observed, version-consistent within-column cell-pair counts; electrical parameters, precise synapse placement and missing external drive require modeling assumptions.',
    }
    (args.output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps({k:summary[k] for k in ['selected_neurons','ei_classes','internal_directed_pairs','internal_anatomical_synapses']},indent=2),flush=True)


if __name__=='__main__':main()
