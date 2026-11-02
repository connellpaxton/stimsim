"""Inspect the supplied SWC/connectivity archive without modifying it."""

import argparse
import ast
from collections import Counter
import csv
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from neuronsim.swc import load_swc

DEFAULT_DATA = Path('/home/shisoft/Downloads/neuron_10-26-2026-20261027T002512Z-1-001/neuron_10-26-2026')
SAMPLE = ['864691135640672571', '864691135993429185', '864691136025333561']




def analyze(data, output):
    output.mkdir(parents=True, exist_ok=True)
    paths = sorted(data.glob('SWC_*/*.swc'))
    if not paths:
        raise ValueError(f'No SWC morphologies found in {data}')
    by_id = {p.stem: p for p in paths}
    layers = {}
    invalid, morphology_rows, sample_cells = [], [], []
    valid_ids = set()
    for path in paths:
        layer = path.parent.name.removeprefix('SWC_')
        stats = layers.setdefault(layer, {'cells': 0, 'positive_radius_cells': 0, 'nodes': 0, 'nonpositive_radius_nodes': 0, 'cable_length_um': 0, 'zero_length_edges': 0, 'types': Counter()})
        try:
            rows, parent, order = load_swc(path, require_positive_radii=False)
        except ValueError as error:
            invalid.append(str(error))
            continue
        bad_radii = int((rows[:, 5] <= 0).sum())
        stats['nonpositive_radius_nodes'] += bad_radii
        if bad_radii:
            invalid.append(f'{path}: {bad_radii} nodes with non-positive radii')
        else:
            stats['positive_radius_cells'] += 1
            valid_ids.add(path.stem)
        edge = np.flatnonzero(parent >= 0)
        lengths = np.linalg.norm(rows[edge, 2:5] - rows[parent[edge], 2:5], axis=1)
        stats['cells'] += 1
        stats['nodes'] += len(rows)
        stats['cable_length_um'] += float(lengths.sum())
        stats['zero_length_edges'] += int((lengths == 0).sum())
        stats['types'].update(map(int, rows[:, 1]))
        morphology_rows.append({'cell_id': path.stem, 'layer': layer, 'nodes': len(rows), 'nonpositive_radius_nodes': bad_radii, 'cable_length_um': float(lengths.sum())})
        if path.stem in SAMPLE:
            sample_cells.append({'id': path.stem, 'xyz': rows[:, 2:5].round(3).tolist(), 'parent': parent.tolist(), 'types': rows[:, 1].astype(int).tolist()})
    tables, sample_contacts = {}, []
    contacts_seen, duplicate_contacts, internal_unique, validated_unique = set(), 0, set(), set()
    layer_pairs = Counter()
    for path in sorted(data.glob('*_all.csv')):
        with path.open(newline='') as stream:
            records = list(csv.DictReader(stream))
        ids = set()
        internal, within_layer = 0, 0
        pairs = Counter()
        for row in records:
            pre, post = row['pre_pt_root_id'], row['post_pt_root_id']
            ids.update((pre, post))
            key = tuple(row.values())
            if key in contacts_seen:
                duplicate_contacts += 1
            contacts_seen.add(key)
            if pre in by_id and post in by_id:
                internal += 1
                internal_unique.add(key)
                if pre in valid_ids and post in valid_ids:
                    validated_unique.add(key)
                pairs[(pre, post)] += 1
                a, b = by_id[pre].parent.name[4:], by_id[post].parent.name[4:]
                layer_pairs[(a, b)] += 1
                within_layer += a == b == path.stem.removesuffix('_all')
            if pre in SAMPLE and post in SAMPLE:
                coords = {side: [float(row[f'{side}_pt_position_{axis}']) / 1000 for axis in 'xyz'] for side in ('pre', 'post')}
                distances = {}
                for side, cell_id in [('pre', pre), ('post', post)]:
                    cell = next(c for c in sample_cells if c['id'] == cell_id)
                    distances[side] = float(np.linalg.norm(np.asarray(cell['xyz']) - coords[side], axis=1).min())
                sample_contacts.append({'pre': pre, 'post': post, 'table': path.name, 'xyz_um': coords, 'nearest_node_distance_um': distances})
        tables[path.name] = {'contacts': len(records), 'referenced_cells': len(ids), 'available_cells': len(ids & by_id.keys()), 'contacts_with_both_morphologies': internal, 'contacts_with_both_in_named_layer': within_layer, 'directed_pairs_with_both_morphologies': len(pairs)}
    source = ast.parse((data / 'cell.py').read_text())
    mechanisms = sorted({n.args[0].value for n in ast.walk(source) if isinstance(n, ast.Call) and isinstance(n.func, ast.Attribute) and n.func.attr == 'insert' and n.args and isinstance(n.args[0], ast.Constant)})
    summary = {'source': str(data.resolve()), 'morphology_count': len(paths), 'unique_cell_ids': len(by_id), 'validated_morphology_count': len(valid_ids), 'layers': layers, 'invalid_morphologies': invalid, 'tables': tables, 'unique_contacts': len(contacts_seen), 'duplicate_contact_records': duplicate_contacts, 'unique_contacts_with_both_morphologies': len(internal_unique), 'unique_contacts_with_both_validated_morphologies': len(validated_unique), 'layer_contacts': {f'{a}->{b}': n for (a, b), n in sorted(layer_pairs.items())}, 'sample_contacts': sample_contacts, 'required_custom_mechanisms': [m for m in mechanisms if m not in ('pas', 'extracellular')], 'supplied_mod_files': [str(p.relative_to(data)) for p in data.rglob('*.mod')]}
    (output / 'analysis.json').write_text(json.dumps(summary, indent=2) + '\n')
    with (output / 'morphologies.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=['cell_id', 'layer', 'nodes', 'nonpositive_radius_nodes', 'cable_length_um'])
        writer.writeheader()
        writer.writerows(morphology_rows)
    passive_path = output / 'passive_summary.json'
    payload = {'summary': summary, 'cells': sorted(sample_cells, key=lambda c: SAMPLE.index(c['id'])), 'passive': json.loads(passive_path.read_text()) if passive_path.exists() else None}
    template = Path(__file__).with_name('explorer.html').read_text()
    (output / 'explorer.html').write_text(template.replace('__DATA__', json.dumps(payload, separators=(',', ':'))))
    print(json.dumps({k: summary[k] for k in ('morphology_count', 'validated_morphology_count', 'unique_contacts', 'unique_contacts_with_both_morphologies', 'unique_contacts_with_both_validated_morphologies')}, indent=2))
    return summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', type=Path, default=DEFAULT_DATA)
    parser.add_argument('--output', type=Path, default=Path('outputs'))
    args = parser.parse_args()
    analyze(args.data, args.output)
