"""Check numerical agreement and resource use between two Minnie thread counts."""

import argparse
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np


def digest(path):
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', type=Path, required=True)
    parser.add_argument('--candidate', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = {'baseline': str(args.baseline), 'candidate': str(args.candidate), 'conditions': {}}
    for filename in ('cell_models.csv', 'modeled_sites.csv.gz'):
        if digest(args.baseline / filename) != digest(args.candidate / filename):
            raise RuntimeError(f'Neuron assignments or contacts differ: {filename}')
    result['identical_assignments_and_decompressed_contact_table'] = True
    fields = ('hostname', 'neuron_version', 'source_commit', 'profile_sha256', 'geometry_manifest_sha256', 'input_hashes', 'neurons', 'sections', 'segments',
              'temperature_c', 'dt_ms', 'v_init_mv', 'conductance_us', 'inhibitory_outputs_enabled',
              'contacts_by_presynaptic_anatomical_class', 'synaptic_reversal_mv', 'synaptic_tau1_ms', 'synaptic_tau2_ms',
              'synaptic_delay_ms', 'placement_seed', 'stimulus', 'duration_ms', 'weight_factor', 'mpi_ranks',
              'zero_conductance_segment_instances_removed')
    for condition in ('original_stimulus', 'original_stimulus_control', 'matched_protocol', 'matched_protocol_control'):
        left, right = args.baseline / condition, args.candidate / condition
        a, b = (json.loads((path / 'summary.json').read_text()) for path in (left, right))
        for key in fields:
            if a[key] != b[key]: raise RuntimeError(f'{condition}: configuration differs at {key}')
        for summary in (a, b):
            loads = summary['segments_per_thread']
            if len(loads) != summary['threads'] or sum(loads) != summary['segments']:
                raise RuntimeError('Invalid native thread-load accounting')
        with np.load(left / 'recordings.npz') as x, np.load(right / 'recordings.npz') as y:
            for key in ('root_ids', 'time_ms'):
                if not np.array_equal(x[key], y[key]): raise RuntimeError(f'Recording axes differ: {key}')
            voltage_a, voltage_b = x['soma_mv'], y['soma_mv']
            if voltage_a.shape != voltage_b.shape or not np.isfinite(voltage_a).all() or not np.isfinite(voltage_b).all():
                raise RuntimeError('Invalid voltage recordings')
            delta = voltage_b - voltage_a
            spike_times_identical = bool(np.array_equal(x['spike_times_ms'], y['spike_times_ms']))
            spike_roots_identical = bool(np.array_equal(x['spike_root_ids'], y['spike_root_ids']))
            timings = {}
            for name, summary in (('baseline', a), ('candidate', b)):
                timings[name] = {key: summary[key] for key in ('threads', 'construction_seconds', 'initialization_seconds',
                                                               'solve_seconds', 'resident_ram_bytes_after_simulation',
                                                               'recorded_neuronal_spikes', 'nonstimulated_cell_spikes')}
                timings[name]['segments_per_thread_min_max'] = [min(summary['segments_per_thread']), max(summary['segments_per_thread'])]
            result['conditions'][condition] = {'same_model_data_and_protocol': True,
                                               'voltage_samples': int(voltage_a.size),
                                               'voltages_identical': bool(np.array_equal(voltage_a, voltage_b)),
                                               'maximum_absolute_voltage_difference_mv': float(np.max(np.abs(delta))),
                                               'rms_voltage_difference_mv': float(np.sqrt(np.mean(delta ** 2))),
                                               'spike_times_identical': spike_times_identical,
                                               'spike_neuron_ids_identical': spike_roots_identical,
                                               'solver_speedup': a['solve_seconds'] / b['solve_seconds'],
                                               **timings}
    result['all_recordings_identical'] = all(row['voltages_identical'] and row['spike_times_identical'] and row['spike_neuron_ids_identical']
                                              for row in result['conditions'].values())
    result['timing_note'] = 'Single trial per condition; solver timings exclude construction, initialization and file export. RAM values are resident snapshots, not peaks.'
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__': main()
