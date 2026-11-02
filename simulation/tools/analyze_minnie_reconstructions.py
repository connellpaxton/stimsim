"""Compare reconstructed-column recordings and plot observed branch geometry."""

import argparse
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workspace', type=Path, required=True)
    parser.add_argument('--run', default='reconstructed32')
    args = parser.parse_args()
    base = args.workspace
    run = base / 'outputs' / args.run
    morphology = base / 'data/reconstructions'
    manifest = json.loads((morphology / 'manifest.json').read_text())
    with (run / 'cell_models.csv').open() as stream:
        cells = {r['root_id']: r for r in csv.DictReader(stream)}
    coupled = np.load(run / 'recordings.npz')
    control = np.load(run / 'zero_transmission_control/recordings.npz')
    for key in ('root_ids', 'time_ms'):
        if not np.array_equal(coupled[key], control[key]):
            raise RuntimeError('Recording identities or time axes disagree')
    if not np.array_equal(coupled['soma_mv'][:, 0], control['soma_mv'][:, 0]):
        raise RuntimeError('Initial voltages disagree')
    roots = coupled['root_ids']; time = coupled['time_ms']
    stimulated = np.array([cells[r]['stimulated'] == 'True' for r in roots])
    delta = coupled['soma_mv'] - control['soma_mv']
    maximum = np.max(np.abs(delta), axis=1)
    result = {'same_recording_ids_and_times': True, 'same_initial_soma_voltages': True,
              'recorded_samples_per_condition': int(delta.size),
              'unstimulated_cells': int((~stimulated).sum()),
              'unstimulated_cells_max_delta_over_0_01_mv': int(((maximum > .01) & ~stimulated).sum()),
              'unstimulated_cells_max_delta_over_0_1_mv': int(((maximum > .1) & ~stimulated).sum()),
              'maximum_absolute_soma_delta_mv': float(maximum.max()),
              'unstimulated_rms_soma_delta_mv': float(np.sqrt(np.mean(delta[~stimulated] ** 2)))}
    (run / 'comparison.json').write_text(json.dumps(result, indent=2) + '\n')

    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
    fig, axes = plt.subplots(1, 4, figsize=(14, 7), layout='constrained')
    colors = {1: '#222222', 2: '#43995d', 3: '#2878b5', 4: '#ce6b28'}
    for ax, kind in zip(axes, ('23P', '4P', '5P-PT', 'MC')):
        candidates = sorted([r for r in manifest['records'] if r['manual_cell_type'] == kind], key=lambda r: r['nodes'])
        record = candidates[len(candidates) // 2]
        rows = np.loadtxt(morphology / 'swc' / (record['root_id'] + '.swc'))
        indices = {int(r[0]): i for i, r in enumerate(rows)}
        parent = np.array([-1 if r[6] == -1 else indices[int(r[6])] for r in rows])
        soma = rows[parent == -1][0, 2:5]
        coordinates = rows[:, 2:5] - soma
        for region in (2, 3, 4):
            children = np.flatnonzero((rows[:, 1] == region) & (parent >= 0))
            lines = np.stack((coordinates[parent[children], :2], coordinates[children, :2]), axis=1)
            ax.add_collection(LineCollection(lines, colors=colors[region], linewidths=.7 if region != 2 else .5,
                                             alpha=.8 if region != 2 else .7))
        ax.scatter([0], [0], s=25, color=colors[1], zorder=3)
        ax.autoscale(); ax.set_aspect('equal'); ax.invert_yaxis()
        ax.set_title(f'{kind}\n{record["nodes"]:,} skeleton points', fontsize=11)
        ax.set_xlabel('x from soma (µm)')
        ax.set_ylabel('y from soma (µm)')
        ax.text(.02, -.17, record['root_id'], transform=ax.transAxes, fontsize=7)
    fig.suptitle('Individual Minnie reconstructions used as model geometry\nBlue: basal dendrites · orange: apical dendrites · green: recorded axon', fontsize=13)
    fig.savefig(run / 'recorded_branches.png', dpi=170)
    plt.close(fig)

    selected = int(np.argmax(np.where(~stimulated, maximum, -np.inf)))
    fig, axes = plt.subplots(3, 1, figsize=(11, 9), layout='constrained')
    order = sorted(range(len(roots)), key=lambda i: (cells[roots[i]]['manual_type'], roots[i]))
    scale = max(float(np.quantile(np.abs(delta), .995)), .01)
    heat = axes[0].imshow(delta[order], aspect='auto', origin='lower', extent=[time[0], time[-1], 0, len(roots)],
                          cmap='RdBu_r', vmin=-scale, vmax=scale)
    fig.colorbar(heat, ax=axes[0], label='Coupled − control (mV)')
    axes[0].set_ylabel('Neuron (sorted by type)')
    axes[0].set_title('Effect of transmission across the recorded network')
    axes[1].plot(time, coupled['soma_mv'][selected], label='Coupled', color='#2878b5')
    axes[1].plot(time, control['soma_mv'][selected], label='Zero transmission', color='#888888', linestyle='--')
    axes[1].set_ylabel('Soma voltage (mV)'); axes[1].legend(loc='best')
    axes[1].set_title(f'Largest-response uninjected cell: {cells[roots[selected]]["manual_type"]}, {roots[selected]}')
    index = {root: i for i, root in enumerate(roots)}
    for ei, color in [('E', '#2878b5'), ('I', '#c54d4d')]:
        mask = np.array([cells[r]['ei_class'] == ei for r in coupled['spike_root_ids']], dtype=bool)
        axes[2].scatter(coupled['spike_times_ms'][mask], [index[r] for r in coupled['spike_root_ids'][mask]],
                        s=9, color=color, label=ei)
    axes[2].set_ylabel('Neuron'); axes[2].legend(loc='best'); axes[2].set_title('Coupled somatic spike raster')
    for i, ax in enumerate(axes):
        ax.set_xlim(time[0], time[-1]); ax.set_xlabel('Time (ms)')
        if i > 0: ax.axvspan(10, 60, color='#e8c46e', alpha=.12)
    fig.savefig(run / 'network_response.png', dpi=160)
    plt.close(fig)

    probe_dir = base / 'outputs/branch_probe'
    if probe_dir.exists():
        probe = np.load(probe_dir / 'recordings.npz')
        info = json.loads((probe_dir / 'summary.json').read_text())
        fig, axes = plt.subplots(1, 2, figsize=(10, 4), layout='constrained')
        for i, (condition, color) in enumerate((('proximal', '#2878b5'), ('distal', '#ce6b28')), 1):
            response = probe[condition] - probe['baseline']
            distance = info['sites'][i - 1]['soma_path_um']
            label = f'{condition.capitalize()} ({distance:.0f} µm path)'
            axes[0].plot(probe['time_ms'], response[i], color=color, label=label)
            axes[1].plot(probe['time_ms'], response[0], color=color, label=label)
        for ax, title in zip(axes, ('At the input branch', 'At the soma')):
            ax.set_title(title); ax.set_xlim(45, 90); ax.set_xlabel('Time (ms)'); ax.set_ylabel('Change from baseline (mV)')
            ax.legend(fontsize=8)
        fig.suptitle(f'Same 1 nS AMPA input at two sites on one reconstructed {info["manual_type"]} neuron')
        fig.savefig(run / 'branch_signal_transfer.png', dpi=170)
        plt.close(fig)
        fine_path = base / 'outputs/branch_probe_fine/summary.json'
        if fine_path.exists():
            fine = json.loads(fine_path.read_text())
            if fine['root_id'] != info['root_id'] or fine['swc_sha256'] != info['swc_sha256']:
                raise RuntimeError('Probe refinement used different neuron geometry')
            convergence = {'coarse_segments': info['segments'], 'fine_segments': fine['segments'], 'sites': []}
            for coarse_site, fine_site in zip(info['sites'], fine['sites']):
                if coarse_site['section'] != fine_site['section'] or not np.isclose(coarse_site['x'], fine_site['x']):
                    raise RuntimeError('Probe refinement changed the synaptic site')
                row = {'condition': coarse_site['condition'], 'soma_path_um': coarse_site['soma_path_um']}
                for location in ('soma', 'local'):
                    a, b = (site[f'{location}_peak_delta_mv'] for site in (coarse_site, fine_site))
                    row[f'coarse_{location}_peak_mv'] = a
                    row[f'fine_{location}_peak_mv'] = b
                    row[f'{location}_relative_peak_change_percent'] = 100 * abs(b - a) / abs(b)
                convergence['sites'].append(row)
            (probe_dir / 'convergence.json').write_text(json.dumps(convergence, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
