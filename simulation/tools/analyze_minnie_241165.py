"""Summarize the original-family Minnie experiment and its transmission controls."""

import argparse
from collections import Counter
import csv
import gzip
from itertools import zip_longest
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def compare_polarity(run, reference, cells):
    """Verify the paired configurations differ only in inhibitory reversal."""
    checks = {'reference_run': str(reference), 'changed_contact_reversals': 0, 'contacts_by_presynaptic_class': {}, 'protocols': {}}
    with gzip.open(run / 'modeled_sites.csv.gz', 'rt') as left, gzip.open(reference / 'modeled_sites.csv.gz', 'rt') as right:
        counts = Counter()
        for row, old in zip_longest(csv.DictReader(left), csv.DictReader(right)):
            if row is None or old is None: raise RuntimeError('Contact counts differ from the all-excitatory reference')
            if {k: v for k, v in row.items() if k != 'reversal_mv'} != {k: v for k, v in old.items() if k != 'reversal_mv'}:
                raise RuntimeError('Contact placement or a parameter other than reversal changed')
            kind = cells[row['pre_root_id']]['ei_class']; counts[kind] += 1
            if float(old['reversal_mv']) != 0.: raise RuntimeError('Reference contains inhibitory reversals')
            if kind == 'E' and float(row['reversal_mv']) != 0.: raise RuntimeError('Excitatory output reversal changed')
            if kind == 'I' and float(row['reversal_mv']) >= 0.: raise RuntimeError('Inhibitory output lacks a negative reversal')
            checks['changed_contact_reversals'] += row['reversal_mv'] != old['reversal_mv']
    checks['contacts_by_presynaptic_class'] = dict(counts)
    checks['same_contacts_positions_strengths_and_delays'] = True
    fields = ('source_commit', 'profile_sha256', 'geometry_manifest_sha256', 'input_hashes', 'neurons', 'sections', 'segments',
              'temperature_c', 'dt_ms', 'v_init_mv', 'conductance_us', 'synaptic_tau1_ms', 'synaptic_tau2_ms',
              'synaptic_delay_ms', 'placement_seed', 'stimulus', 'duration_ms', 'weight_factor', 'segments_per_thread')
    for protocol in ('original_stimulus', 'matched_protocol'):
        if not (run / protocol).exists(): continue
        summary = json.loads((run / protocol / 'summary.json').read_text())
        old_summary = json.loads((reference / protocol / 'summary.json').read_text())
        if not summary.get('inhibitory_outputs_enabled') or old_summary.get('inhibitory_outputs_enabled'):
            raise RuntimeError('Expected inhibitory and all-excitatory configurations')
        for key in fields:
            if summary[key] != old_summary[key]: raise RuntimeError(f'Configuration changed beyond polarity: {key}')
        a = np.load(run / protocol / 'recordings.npz'); b = np.load(reference / protocol / 'recordings.npz')
        for key in ('root_ids', 'time_ms'):
            if not np.array_equal(a[key], b[key]): raise RuntimeError('Reference recording axes differ')
        if not np.array_equal(a['soma_mv'][:, 0], b['soma_mv'][:, 0]): raise RuntimeError('Reference initial soma voltages differ')
        control = np.load(run / (protocol + '_control') / 'recordings.npz')
        old_control = np.load(reference / (protocol + '_control') / 'recordings.npz')
        if any(not np.array_equal(control[key], old_control[key]) for key in control.files):
            raise RuntimeError('Zero-transmission control changed with polarity')
        inhibitory_ids = {root for root, row in cells.items() if row['ei_class'] == 'I'}
        output = {'all_other_checked_parameters_identical': True, 'zero_transmission_recordings_identical': True,
                  'all_excitatory_spikes': old_summary['recorded_neuronal_spikes'], 'EI_spikes': summary['recorded_neuronal_spikes'],
                  'all_excitatory_spiking_cells': old_summary['cells_with_spikes'], 'EI_spiking_cells': summary['cells_with_spikes'],
                  'inhibitory_source_spikes': sum(root in inhibitory_ids for root in a['spike_root_ids']),
                  'all_excitatory_final_soma_above_minus40_mv': int((b['soma_mv'][:, -1] > -40).sum()),
                  'EI_final_soma_above_minus40_mv': int((a['soma_mv'][:, -1] > -40).sum())}
        checks['protocols'][protocol] = output
        selected = set(summary['stimulus']['cells'])
        difference = np.max(np.abs(a['soma_mv'] - b['soma_mv']), axis=1)
        i = int(np.argmax(np.where([root not in selected for root in a['root_ids']], difference, -np.inf)))
        fig, axes = plt.subplots(2, 1, figsize=(11, 6), layout='constrained')
        index = {root: j for j, root in enumerate(a['root_ids'])}
        for data, label, color in ((b, 'All contacts excitatory', '#ce6b28'), (a, 'Minnie E/I outputs', '#2878b5')):
            axes[0].plot(data['time_ms'], data['soma_mv'][i], color=color, label=label)
            axes[1].scatter(data['spike_times_ms'], [index[root] for root in data['spike_root_ids']], color=color, s=6, alpha=.6, label=label)
        axes[0].set_title(f'Same recorded neuron, different inhibitory reversal: {a["root_ids"][i]}')
        axes[0].set_ylabel('Soma voltage (mV)'); axes[1].set_ylabel('Neuron')
        axes[1].set_title('Population spike comparison')
        for ax in axes:
            ax.legend(loc='best'); ax.set_xlabel('Time (ms)'); ax.set_xlim(a['time_ms'][0], a['time_ms'][-1])
        fig.savefig(run / protocol / 'inhibition_comparison.png', dpi=160, bbox_inches='tight'); plt.close(fig)
    (run / 'polarity_comparison.json').write_text(json.dumps(checks, indent=2) + '\n')
    return checks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workspace', type=Path, required=True)
    parser.add_argument('--run', default='aberra241165_32')
    parser.add_argument('--reference-run', help='All-excitatory run used to verify the isolated polarity change')
    args = parser.parse_args()
    base = args.workspace; run = base / 'outputs' / args.run
    with (run / 'cell_models.csv').open() as stream:
        cells = {row['root_id']: row for row in csv.DictReader(stream)}
    comparisons = {}
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
    for protocol in ('original_stimulus', 'matched_protocol'):
        folder = run / protocol
        if not folder.exists(): continue
        summary = json.loads((folder / 'summary.json').read_text())
        control_summary = json.loads((run / (protocol + '_control') / 'summary.json').read_text())
        coupled = np.load(folder / 'recordings.npz')
        control = np.load(run / (protocol + '_control') / 'recordings.npz')
        for key in ('root_ids', 'time_ms'):
            if not np.array_equal(coupled[key], control[key]): raise RuntimeError('Control recording identities or times differ')
        if not np.array_equal(coupled['soma_mv'][:, 0], control['soma_mv'][:, 0]):
            raise RuntimeError('Initial soma voltages differ')
        if summary['input_hashes'] != control_summary['input_hashes']:
            raise RuntimeError('Coupled and control graphs differ')
        roots, time = coupled['root_ids'], coupled['time_ms']; stimulated = set(summary['stimulus']['cells'])
        mask = np.array([root not in stimulated for root in roots])
        delta = coupled['soma_mv'] - control['soma_mv']; maximum = np.max(np.abs(delta), axis=1)
        nonstim_roots = set(coupled['spike_root_ids'].tolist()) - stimulated
        comparisons[protocol] = {'same_input_graph_initial_soma_voltages_and_recording_axes': True,
                                 'recorded_voltage_samples_per_condition': int(delta.size),
                                 'uninjected_spiking_cells': len(nonstim_roots),
                                 'uninjected_spiking_cells_by_anatomical_class': dict(Counter(cells[r]['ei_class'] for r in nonstim_roots)),
                                 'uninjected_cells_delta_over_0_01_mv': int(((maximum > .01) & mask).sum()),
                                 'uninjected_cells_delta_over_0_1_mv': int(((maximum > .1) & mask).sum()),
                                 'maximum_absolute_soma_delta_mv': float(maximum.max()),
                                 'coupled_spikes': summary['recorded_neuronal_spikes'],
                                 'control_spikes': control_summary['recorded_neuronal_spikes'],
                                 'control_uninjected_spikes': control_summary['nonstimulated_cell_spikes'],
                                 'coupled_final_soma_above_minus40_mv': int((coupled['soma_mv'][:, -1] > -40).sum()),
                                 'control_final_soma_above_minus40_mv': int((control['soma_mv'][:, -1] > -40).sum())}
        selected = int(np.argmax(np.where(mask, maximum, -np.inf)))
        fig, axes = plt.subplots(3, 1, figsize=(11, 9), layout='constrained')
        scale = max(.01, float(np.quantile(np.abs(delta), .995)))
        heat = axes[0].imshow(delta, origin='lower', aspect='auto', extent=[time[0], time[-1], 0, len(roots)],
                              cmap='RdBu_r', vmin=-scale, vmax=scale)
        fig.colorbar(heat, ax=axes[0], label='Coupled − control (mV)')
        axes[0].set_ylabel('Neuron'); axes[0].set_title('Published 241165 regional model: effect of synaptic transmission')
        axes[1].plot(time, coupled['soma_mv'][selected], color='#2878b5', label='Coupled')
        axes[1].plot(time, control['soma_mv'][selected], color='#888888', linestyle='--', label='Zero transmission')
        axes[1].set_ylabel('Soma voltage (mV)'); axes[1].legend(loc='best')
        axes[1].set_title(f'Largest-response uninjected neuron: {cells[roots[selected]]["manual_type"]}, {roots[selected]}')
        index = {root: i for i, root in enumerate(roots)}
        for ei, color in [('E', '#2878b5'), ('I', '#c54d4d')]:
            event_mask = np.array([cells[r]['ei_class'] == ei for r in coupled['spike_root_ids']], dtype=bool)
            axes[2].scatter(coupled['spike_times_ms'][event_mask], [index[r] for r in coupled['spike_root_ids'][event_mask]],
                            s=8, color=color, label=f'Anatomical {ei}')
        axes[2].set_ylabel('Neuron'); axes[2].legend(loc='best')
        axes[2].set_title('Somatic spike raster — Minnie labels determine E/I outputs' if summary.get('inhibitory_outputs_enabled') else 'Somatic spike raster — all contacts modeled as excitatory')
        stimulus = summary['stimulus']
        for i, ax in enumerate(axes):
            ax.set_xlim(time[0], time[-1]); ax.set_xlabel('Time (ms)')
            if i > 0: ax.axvspan(stimulus['onset_ms'], stimulus['onset_ms'] + stimulus['duration_ms'], color='#e8c46e', alpha=.12)
        fig.savefig(folder / 'network_response.png', dpi=160, bbox_inches='tight'); plt.close(fig)
    if (run / 'matched_protocol/summary.json').exists() and (base / 'outputs/reconstructed32/summary.json').exists():
        a = json.loads((run / 'matched_protocol/summary.json').read_text())
        b = json.loads((base / 'outputs/reconstructed32/summary.json').read_text())
        comparisons['M1_comparison'] = {'same_graph': a['input_hashes'] == b['input_hashes'],
                                      'same_recorded_geometry_manifest': a['geometry_manifest_sha256'] == b['geometry']['manifest_sha256'],
                                      'same_injected_neurons': set(a['stimulus']['cells']) == set(b['stimulus']['cells']),
                                      'same_current_amplitude': a['stimulus']['amplitude_na'] == b['stimulus']['amplitude_na'],
                                      'same_injection_interval': all(a['stimulus'][key] == b['stimulus'][key] for key in ('onset_ms', 'duration_ms')),
                                      'interpretation': 'Check the matching flags before comparing. Even with matched injected currents and anatomy, channel models, synaptic conductance/kinetics/polarity/delays, initial voltage, discretization and placements differ; the experiment does not isolate ion channels alone.'}
    (run / 'comparison.json').write_text(json.dumps(comparisons, indent=2) + '\n')
    if args.reference_run:
        comparisons['polarity_comparison'] = compare_polarity(run, base / 'outputs' / args.reference_run, cells)
    print(json.dumps(comparisons, indent=2))


if __name__ == '__main__': main()
