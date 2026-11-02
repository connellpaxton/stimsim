"""Render recorded Minnie soma and cable voltages as a propagation video."""

import argparse
import csv
import html
import json
from pathlib import Path
import time

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
import numpy as np


def playback_times(duration):
    early = min(200., duration)
    return np.r_[np.arange(0., early, .2), np.arange(early, duration, 2.), duration]


def write_player(output, summary, metadata):
    """Provide a local browser player alongside the standalone MP4."""
    player = output.with_suffix('.html')
    player.write_text(f'''<!doctype html>
<html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Minnie signal propagation</title>
<style>
body{{margin:0;background:#0b1220;color:#e5edf7;font:16px system-ui,sans-serif}}
main{{max-width:1400px;margin:28px auto;padding:0 24px}}h1{{font-size:28px;margin-bottom:8px}}
p{{color:#b8c8dc;line-height:1.55}}video{{width:100%;display:block;background:#080e18;border:1px solid #33435c;border-radius:8px}}
.facts{{display:flex;gap:20px;flex-wrap:wrap;margin:14px 0 20px}}.facts span{{padding:9px 14px;background:#172236;border-radius:6px}}
</style><main><h1>Stimulated signal propagation in the Minnie column</h1>
<p>Recorded simulation voltages on the actual reconstructed branches and within-column wiring.</p>
<div class="facts"><span>{summary['duration_ms']/1000:g} second of model time</span>
<span>{summary['neurons']:,} neurons · {summary['anatomical_contact_count']:,} contacts</span>
<span>{metadata['video_duration_seconds']:.1f} seconds of slowed playback</span></div>
<video controls autoplay muted playsinline preload="metadata" poster="frame_0014.00ms.png">
<source src="{html.escape(output.name, quote=True)}" type="video/mp4"></video>
<p>Left: soma voltages throughout the column. Right: cable-segment voltages on one injected neuron and connected excitatory and inhibitory recipients. Gold rings mark injected cells; white halos mark spikes in the preceding 2 ms.</p>
<p>Input: {summary['stimulus']['amplitude_na']:g} nA into {len(summary['stimulus']['cells'])} neurons, from {summary['stimulus']['onset_ms']:g} to {summary['stimulus']['onset_ms']+summary['stimulus']['duration_ms']:g} ms. The first 200 ms plays more slowly; the on-screen clock always shows model time.</p>
<p>Geometry, connectivity and E/I labels come from Minnie. Electrical parameters use ModelDB 241165; synaptic strengths, timing and attachment locations are modeling assumptions.</p>
</main></html>''')
    return player


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workspace', type=Path, required=True)
    parser.add_argument('--run', default='aberra241165_ei96_1s')
    parser.add_argument('--output', type=Path)
    parser.add_argument('--fps', type=int, default=30)
    parser.add_argument('--max-frames', type=int, help='Render only an initial subset for a smoke test')
    args = parser.parse_args()
    base = args.workspace; run = base / 'outputs' / args.run; folder = run / 'long_protocol'
    output = args.output or run / 'video/propagation.mp4'; output.parent.mkdir(parents=True, exist_ok=True)
    with (base / 'data/column/nodes.csv').open() as stream: nodes = {row['root_id']: row for row in csv.DictReader(stream)}
    summary = json.loads((folder / 'summary.json').read_text())
    geometry = json.loads((run / 'branch_geometry.json').read_text())
    data = np.load(folder / 'recordings.npz'); branch = np.load(folder / 'branch_recordings.npz')
    roots = data['root_ids']; selected = geometry['roots']; times = data['time_ms']
    soma_mv = data['soma_mv']; segment_mv = branch['segment_mv']
    if not np.array_equal(times, branch['time_ms']): raise RuntimeError('Branch and soma time axes differ')
    if segment_mv.shape[0] != len(geometry['rows']): raise RuntimeError('Branch geometry and recording rows disagree')
    index = {root: i for i, root in enumerate(roots)}
    for root in selected:
        i = next(i for i, row in enumerate(geometry['rows']) if row['root_id'] == root and row['section'] == 'soma_0' and row['x'] == .5)
        if not np.array_equal(segment_mv[i], soma_mv[index[root]].astype(np.float32)):
            raise RuntimeError('Native branch soma recording does not match the population recording')
    xyz = np.array([[float(nodes[root][a + '_um']) for a in 'xyz'] for root in roots])
    angle = np.deg2rad(25.)
    projection = np.array([[np.cos(angle), 0., np.sin(angle)], [0., 1., 0.]])
    origin = xyz.mean(axis=0); xy = (xyz - origin) @ projection.T
    spike_times = data['spike_times_ms']; spike_ids = np.array([index[root] for root in data['spike_root_ids']], dtype=int)
    ei = np.array([nodes[root]['ei_class'] for root in roots]); stimulated = set(summary['stimulus']['cells'])
    stim_indices = [index[root] for root in roots if root in stimulated]
    frame_times = playback_times(summary['duration_ms'])
    if args.max_frames: frame_times = frame_times[:args.max_frames]
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10, 'text.color': '#e5edf7',
                         'axes.labelcolor': '#b8c8dc', 'xtick.color': '#9dafc5', 'ytick.color': '#9dafc5',
                         'axes.edgecolor': '#41506a', 'axes.facecolor': '#0b1220', 'savefig.facecolor': '#0b1220'})
    fig = plt.figure(figsize=(16, 10), facecolor='#0b1220')
    fig.text(.045, .956, 'MINNIE  |  STIMULATED SIGNAL PROPAGATION', fontsize=21, weight='bold')
    fig.text(.045, .923, f'{len(roots):,} reconstructed neurons  ·  {summary["anatomical_contact_count"]:,} anatomical contacts  ·  ModelDB 241165 + E/I outputs', fontsize=12)
    clock = fig.text(.045, .884, '', fontsize=15, color='#f7c65f')
    speed = fig.text(.96, .884, '', fontsize=10, ha='right')
    population = fig.add_axes([.055, .14, .235, .69])
    population.set_title('Whole column · soma voltage', pad=12)
    population.set_aspect('equal'); population.invert_yaxis()
    population.set_xlabel('Oblique horizontal position (µm)'); population.set_ylabel('Cortical depth projection (µm)')
    norm = Normalize(-90., 40.); cmap = plt.get_cmap('turbo')
    populations = []
    for kind, marker in [('E', 'o'), ('I', '^')]:
        mask = ei == kind
        artist = population.scatter(xy[mask, 0], xy[mask, 1], c=soma_mv[mask, 0], s=10 if kind == 'E' else 18,
                                    marker=marker, cmap=cmap, norm=norm, edgecolors='none', alpha=.9)
        populations.append((artist, mask))
    population.scatter(xy[stim_indices, 0], xy[stim_indices, 1], s=38, facecolors='none', edgecolors='#f7c65f', linewidths=.8)
    halo = population.scatter([], [], s=45, facecolors='none', edgecolors='white', linewidths=.8)
    population.legend(handles=[Line2D([], [], marker='o', linestyle='', color='#b7c8dc', label='Excitatory'),
                               Line2D([], [], marker='^', linestyle='', color='#b7c8dc', label='Inhibitory'),
                               Line2D([], [], marker='o', linestyle='', markerfacecolor='none', color='#f7c65f', label='Injected')],
                      loc='lower left', fontsize=8, framealpha=.85, facecolor='#172033', edgecolor='#41506a')
    colors = ['#f7c65f', '#7dd3fc', '#d0a5ff']; tree_artists = []
    roles = ['Injected neuron', 'Connected excitatory cell', 'Connected inhibitory cell']
    for k, root in enumerate(selected):
        ax = fig.add_axes([.335 + .207 * k, .42, .184, .40])
        rows = [i for i, row in enumerate(geometry['rows']) if row['root_id'] == root]
        soma = np.array([float(nodes[root][a + '_um']) for a in 'xyz'])
        lines = [(np.array(geometry['rows'][i]['points_um']) - soma) @ projection.T for i in rows]
        ax.add_collection(LineCollection(lines, colors='#5e6a7c', linewidths=1.5, alpha=.3))
        collection = LineCollection(lines, cmap=cmap, norm=norm, linewidths=.85)
        collection.set_array(segment_mv[rows, 0]); ax.add_collection(collection); tree_artists.append((collection, rows))
        ax.scatter([0], [0], color=colors[k], s=16, zorder=4)
        ax.autoscale(); ax.set_aspect('equal'); ax.invert_yaxis()
        ax.set_title(f'{k + 1}. {roles[k]}\n{nodes[root]["manual_cell_type"]}', fontsize=11, color=colors[k], pad=9)
        ax.set_xlabel('µm from soma', fontsize=8); ax.tick_params(labelsize=8)
        population.annotate(str(k + 1), xy=xy[index[root]], xytext=(12 if k != 2 else -22, 14 if k != 2 else -18),
                            textcoords='offset points', color=colors[k], fontsize=13, weight='bold',
                            arrowprops={'arrowstyle': '-', 'color': colors[k], 'linewidth': .7})
    color_axis = fig.add_axes([.965, .44, .008, .34])
    fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=color_axis, ticks=[-90, -65, -40, 0, 40])
    color_axis.set_title('mV', fontsize=9)
    traces = fig.add_axes([.335, .14, .615, .205])
    traces.set_title('Recorded soma responses · selected cells', loc='left', fontsize=11)
    traces.set_ylim(-95, 50); traces.set_xlim(0, min(200., summary['duration_ms']))
    traces.set_xlabel('Model time (ms)'); traces.set_ylabel('mV')
    traces.axvspan(summary['stimulus']['onset_ms'], summary['stimulus']['onset_ms'] + summary['stimulus']['duration_ms'], color='#f7c65f', alpha=.12)
    trace_lines = [traces.plot([], [], color=colors[k], linewidth=1., label=str(k + 1))[0] for k in range(len(selected))]
    cursor = traces.axvline(0, color='white', linewidth=.7, alpha=.7)
    traces.legend(loc='upper right', ncol=3, framealpha=.8, facecolor='#172033', edgecolor='#41506a', fontsize=8)
    activity = fig.text(.055, .085, '', fontsize=10, color='#b8c8dc')
    fig.text(.335, .077, 'Lines: imported 3D cable centerlines · colors: actual simulated segment voltages', fontsize=10)
    fig.text(.045, .035, 'Geometry, wiring and E/I labels: Minnie. Electrical parameters, synaptic strengths/locations and current stimulus are modeled.', fontsize=9, color='#8fa3bd')
    writer = FFMpegWriter(fps=args.fps, codec='libx264', extra_args=['-crf', '20', '-preset', 'fast', '-pix_fmt', 'yuv420p', '-movflags', '+faststart'])
    started = time.perf_counter(); preview_indices = set(np.clip(np.searchsorted(frame_times, [0., 14., 35., 80., 200., summary['duration_ms']]), 0, len(frame_times) - 1).tolist())
    with writer.saving(fig, str(output), dpi=100):
        for frame, model_time in enumerate(frame_times):
            sample = min(int(round(model_time / .1)), len(times) - 1)
            for artist, mask in populations: artist.set_array(soma_mv[mask, sample])
            for artist, rows in tree_artists: artist.set_array(segment_mv[rows, sample])
            lo = np.searchsorted(spike_times, model_time - 2., side='right'); hi = np.searchsorted(spike_times, model_time, side='right')
            active = np.unique(spike_ids[lo:hi]); halo.set_offsets(xy[active] if len(active) else np.empty((0, 2)))
            on = summary['stimulus']['onset_ms'] <= model_time < summary['stimulus']['onset_ms'] + summary['stimulus']['duration_ms']
            clock.set_text(f'Model time: {model_time:7.2f} / {summary["duration_ms"]:g} ms   |   Stimulus {"ON" if on else "OFF"} · {len(stimulated)} cells, {summary["stimulus"]["amplitude_na"]:g} nA')
            speed.set_text(f'Playback: 1 s video = {(.2 if model_time < 200. else 2.) * args.fps:g} ms model time')
            activity.set_text(f'Spikes in preceding 2 ms: {hi - lo}  ·  white halos = recent spikes')
            if model_time >= 200.: traces.set_xlim(0, summary['duration_ms'])
            for k, root in enumerate(selected): trace_lines[k].set_data(times[:sample + 1], soma_mv[index[root], :sample + 1])
            cursor.set_xdata([model_time, model_time])
            writer.grab_frame()
            if frame in preview_indices: fig.savefig(output.parent / f'frame_{model_time:07.2f}ms.png', dpi=100)
            if frame % 100 == 0: print(f'Video {frame}/{len(frame_times)} frames; model {model_time:g} ms; render {time.perf_counter() - started:.1f} s', flush=True)
    plt.close(fig)
    metadata = {'source_run': str(run), 'video': str(output), 'model_duration_ms': summary['duration_ms'], 'frames': len(frame_times),
                'fps': args.fps, 'video_duration_seconds': len(frame_times) / args.fps, 'selected_roots': selected,
                'playback': '0–200 ms at 0.2 ms/frame; remainder at 2 ms/frame. Actual model time is labeled.',
                'encoding': 'H.264, yuv420p', 'render_seconds': time.perf_counter() - started,
                'data': 'Population colors are recorded soma voltage. Branch colors are recorded computational-segment voltage. White halos mark recorded spikes in the preceding 2 ms. No activity or traveling signal arrows are synthesized.'}
    metadata['player'] = str(write_player(output, summary, metadata))
    output.with_suffix('.json').write_text(json.dumps(metadata, indent=2) + '\n')
    print(json.dumps(metadata, indent=2))


if __name__ == '__main__': main()
