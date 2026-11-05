"""Validate the one-second propagation recording against the prior benchmark."""

from collections import Counter
import argparse
import csv
import json
from pathlib import Path

import numpy as np

from compare_minnie_threads import digest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workspace', type=Path, required=True)
    parser.add_argument('--run', default='aberra241165_ei96_1s')
    parser.add_argument('--reference-run', default='aberra241165_ei96')
    args = parser.parse_args(); base = args.workspace
    run = base / 'outputs' / args.run; reference = base / 'outputs' / args.reference_run
    summary = json.loads((run / 'long_protocol/summary.json').read_text())
    old = json.loads((reference / 'matched_protocol/summary.json').read_text())
    for key in ('source_commit', 'profile_sha256', 'geometry_manifest_sha256', 'input_hashes', 'neurons', 'sections', 'segments',
                'temperature_c', 'dt_ms', 'v_init_mv', 'conductance_us', 'synaptic_reversal_mv', 'synaptic_tau1_ms',
                'synaptic_tau2_ms', 'synaptic_delay_ms', 'placement_seed', 'stimulus', 'threads', 'segments_per_thread'):
        if summary[key] != old[key]: raise RuntimeError(f'Model configuration changed: {key}')
    for name in ('cell_models.csv', 'modeled_sites.csv.gz'):
        if digest(run / name) != digest(reference / name): raise RuntimeError(f'Model assignments differ: {name}')
    geometry = json.loads((run / 'branch_geometry.json').read_text())
    data = np.load(run / 'long_protocol/recordings.npz'); branch = np.load(run / 'long_protocol/branch_recordings.npz')
    previous = np.load(reference / 'matched_protocol/recordings.npz')
    if not np.array_equal(data['root_ids'], previous['root_ids']): raise RuntimeError('Recording neuron identities changed')
    prefix = previous['soma_mv'].shape[1]
    delta = data['soma_mv'][:, :prefix] - previous['soma_mv']
    same_voltages = bool(np.array_equal(data['soma_mv'][:, :prefix], previous['soma_mv']))
    mask = data['spike_times_ms'] <= old['duration_ms'] + old['dt_ms'] / 2
    same_spikes = bool(np.array_equal(data['spike_times_ms'][mask], previous['spike_times_ms']) and np.array_equal(data['spike_root_ids'][mask], previous['spike_root_ids']))
    if not np.array_equal(data['time_ms'], branch['time_ms']): raise RuntimeError('Branch and soma time axes disagree')
    if not np.isfinite(data['soma_mv']).all() or not np.isfinite(branch['segment_mv']).all(): raise RuntimeError('Nonfinite recording')
    soma_checks = {}
    for root in geometry['roots']:
        row = next(i for i, entry in enumerate(geometry['rows']) if entry['root_id'] == root and entry['section'] == 'soma_0' and entry['x'] == .5)
        index = list(data['root_ids']).index(root)
        soma_checks[root] = bool(np.array_equal(branch['segment_mv'][row], data['soma_mv'][index].astype(np.float32)))
    if not all(soma_checks.values()): raise RuntimeError('Selected branch/soma recordings are misaligned')
    histogram, bins = np.histogram(data['spike_times_ms'], bins=np.arange(0, summary['duration_ms'] + 100., 100.))
    with (base / 'data/column/nodes.csv').open() as stream: nodes = {row['root_id']: row for row in csv.DictReader(stream)}
    result = {'duration_ms': summary['duration_ms'], 'actual_final_time_ms': summary['actual_final_time_ms'],
              'same_model_assignments_contacts_and_stimulus': True,
              'first_100_ms_voltages_identical_to_prior_run': same_voltages,
              'first_100_ms_max_voltage_difference_mv': float(np.max(np.abs(delta))),
              'first_100_ms_spike_times_and_ids_identical': same_spikes,
              'branch_soma_matches_after_float32_conversion': soma_checks,
              'soma_voltage_samples': int(data['soma_mv'].size), 'branch_voltage_samples': int(branch['segment_mv'].size),
              'recording_interval_ms': .1, 'spike_counts_per_100_ms': histogram.tolist(), 'bin_edges_ms': bins.tolist(),
              'first_spike_ms': float(data['spike_times_ms'].min()) if len(data['spike_times_ms']) else None,
              'last_spike_ms': float(data['spike_times_ms'].max()) if len(data['spike_times_ms']) else None,
              'spikes_by_ei': dict(Counter(nodes[root]['ei_class'] for root in data['spike_root_ids'])),
              'selected_cells': [{**nodes[root], 'recorded_segments': sum(row['root_id'] == root for row in geometry['rows'])} for root in geometry['roots']]}
    (run / 'validation.json').write_text(json.dumps(result, indent=2) + '\n')
    if not same_voltages or not same_spikes: raise RuntimeError('Long recording diverged from the existing 100 ms benchmark; inspect validation.json')
    print(json.dumps(result, indent=2))


if __name__ == '__main__': main()
