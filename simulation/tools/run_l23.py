"""Run a documented active L23 network demonstration using supplied contacts."""

import argparse
from collections import Counter
import csv
import hashlib
import json
import os
from pathlib import Path
import platform
import time

ROOT = Path(__file__).resolve().parents[1]
PUBLISHED = ROOT / 'vendor' / 'aberra2018'
DEFAULT_DATA = Path('/home/shisoft/Downloads/neuron_10-26-2026-20261027T002512Z-1-001/neuron_10-26-2026')


def choose_cells(valid_ids, contacts, limit):
    if limit < 0:
        raise ValueError('Cell limit must be nonnegative')
    outgoing = Counter(row['pre_pt_root_id'] for row in contacts)
    seed = min(valid_ids, key=lambda name: (-outgoing[name], name))
    if limit == 0 or limit >= len(valid_ids):
        return [seed] + sorted(valid_ids - {seed})
    chosen = {seed}
    while len(chosen) < limit:
        scores = Counter()
        for row in contacts:
            a, b = row['pre_pt_root_id'], row['post_pt_root_id']
            if a in chosen and b not in chosen:
                scores[b] += 1
            if b in chosen and a not in chosen:
                scores[a] += 1
        remaining = valid_ids - chosen
        chosen.add(min(remaining, key=lambda name: (-scores[name], name)))
    return [seed] + sorted(chosen - {seed})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', type=Path, default=DEFAULT_DATA)
    parser.add_argument('--output', type=Path, default=ROOT / 'outputs/l23_pilot')
    parser.add_argument('--max-cells', type=int, default=4, help='Connected subset size; 0 selects every valid L23 morphology')
    parser.add_argument('--duration', type=float, default=25., help='Duration in ms')
    parser.add_argument('--amplitude', type=float, default=1., help='Soma pulse amplitude in nA')
    parser.add_argument('--synaptic-weight', type=float, default=.01, help='Peak synaptic conductance in uS; 0 provides a no-transmission control')
    parser.add_argument('--threads', type=int, default=1, help='NEURON solver threads within one process')
    parser.add_argument('--original-family', action='store_true', help='Use the current cACint209 profile and original-wrapper discretization instead of the earlier L23 profile')
    args = parser.parse_args()
    if args.duration <= 0 or args.max_cells < 0 or args.synaptic_weight < 0:
        parser.error('Duration must be positive; cell limit and synaptic weight must be nonnegative')
    if args.threads < 1:
        parser.error('Thread count must be positive')
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=True)
    os.environ['NEURON_MODULE_OPTIONS'] = os.environ.get('NEURON_MODULE_OPTIONS', '') + ' -nogui'
    os.environ['MPLBACKEND'] = 'Agg'
    os.environ['MPLCONFIGDIR'] = str(args.output / 'matplotlib-cache')
    os.environ['OMP_NUM_THREADS'] = '1'
    os.environ['OPENBLAS_NUM_THREADS'] = '1'
    import numpy as np
    import neuron
    from neuron import h
    from analyze import load_swc
    from published_swc import PublishedSWCCell

    beginning = time.perf_counter()
    files = sorted((args.data / 'SWC_L23').glob('*.swc'))
    if not files:
        parser.error('No L23 morphologies found')
    valid, skipped = {}, []
    for path in files:
        try:
            load_swc(path)
        except ValueError as error:
            skipped.append({'cell_id': path.stem, 'reason': str(error)})
        else:
            valid[path.stem] = path
    if not valid:
        parser.error('No valid L23 morphologies')
    connectivity = args.data / 'L23_all.csv'
    with connectivity.open(newline='') as stream:
        all_contacts = list(csv.DictReader(stream))
    supplied_ids = {path.stem for path in files}
    within_layer = [row for row in all_contacts if row['pre_pt_root_id'] in supplied_ids and row['post_pt_root_id'] in supplied_ids]
    valid_contacts = [row for row in within_layer if row['pre_pt_root_id'] in valid and row['post_pt_root_id'] in valid]
    selected = choose_cells(set(valid), valid_contacts, args.max_cells)
    selected_set = set(selected)
    contacts = [row for row in valid_contacts if row['pre_pt_root_id'] in selected_set and row['post_pt_root_id'] in selected_set]
    print(json.dumps({'selected_cells': len(selected), 'selected_contacts': len(contacts), 'valid_cells': len(valid), 'excluded_invalid_cells': len(skipped), 'stimulated_cell': selected[0]}), flush=True)
    if args.original_family:
        from neuronsim.network import load_published_model
        manifest, profile = load_published_model(h, neuron)
        profile_name = 'L4_LBC_cACint209_1'
        profile_url = 'https://modeldb.science/241165'
    else:
        neuron.load_mechanisms(str(PUBLISHED / 'mechanisms'))
        h.load_file('stdrun.hoc')
        h.load_file('import3d.hoc')
        profile_path = PUBLISHED / 'cells/L23_PC_cADpyr229_1/biophysics.hoc'
        manifest = json.loads((PUBLISHED / 'L23_SOURCE_MANIFEST.json').read_text())
        for entry in manifest['files']:
            if hashlib.sha256((ROOT / entry['path']).read_bytes()).hexdigest() != entry['sha256']:
                raise RuntimeError('Published profile source has changed')
        h.load_file(str(profile_path))
        profile = h.cADpyr229_biophys()
        profile_name = 'L23_PC_cADpyr229_1'
        profile_url = manifest['files'][0]['display_url']
    pc = h.ParallelContext()
    # Construct sequentially, then distribute complete cells as in neuronsim.
    pc.nthread(1)
    h.CVode().cache_efficient(1)
    h.cvode_active(0)
    h.celsius = 34
    h.dt = .025
    cells = {}
    for index, cell_id in enumerate(selected, 1):
        cells[cell_id] = PublishedSWCCell(valid[cell_id], profile, original_discretization=args.original_family)
        if index == len(selected) or index % 25 == 0:
            print(f'Loaded {index}/{len(selected)} cells; elapsed {time.perf_counter() - beginning:.1f} s', flush=True)
    # Largest-first placement balances compartment counts without splitting cells.
    partitions = [h.SectionList() for _ in range(args.threads)]
    partition_loads = [0] * args.threads
    for cell in sorted(cells.values(), key=lambda cell: (-len(cell.segments), cell.cell_id)):
        thread = min(range(args.threads), key=lambda index: partition_loads[index])
        for section in cell.sections:
            if section.parentseg() is None:
                partitions[thread].append(sec=section)
        partition_loads[thread] += len(cell.segments)
    pc.nthread(args.threads)
    for thread, roots in enumerate(partitions):
        pc.partition(thread, roots)
    # NEURON 9's Python binding coerces sec_in_thread's integer ID to bool.
    # Returning it through a HOC function preserves thread IDs above one.
    h('objref l23_thread_context')
    h.l23_thread_context = pc
    h('func l23_thread_id() { return l23_thread_context.sec_in_thread() }')
    sources, source_monitors = {}, []

    def record_source(segment):
        key = (segment.sec.name(), float(segment.x))
        if key not in sources:
            monitor = h.NetCon(segment._ref_v, None, sec=segment.sec)
            monitor.threshold = 0
            vector = h.Vector()
            monitor.record(vector)
            sources[key] = vector
            source_monitors.append(monitor)
        return sources[key]

    soma_events = {name: record_source(cell.soma[0](.5)) for name, cell in cells.items()}
    synapses, connections, connection_events, mappings = [], [], [], []
    for row in contacts:
        pre, post = row['pre_pt_root_id'], row['post_pt_root_id']
        pre_segment, pre_distance = cells[pre].nearest_segment([float(row[f'pre_pt_position_{axis}']) for axis in 'xyz'])
        post_segment, post_distance = cells[post].nearest_segment([float(row[f'post_pt_position_{axis}']) for axis in 'xyz'])
        synapse = h.Exp2Syn(post_segment)
        synapse.tau1, synapse.tau2, synapse.e = .5, 3., 0.
        connection = h.NetCon(pre_segment._ref_v, synapse, sec=pre_segment.sec)
        connection.threshold, connection.delay, connection.weight[0] = 0., 1., args.synaptic_weight
        connection_events.append(record_source(pre_segment))
        synapses.append(synapse)
        connections.append(connection)
        mappings.append({'pre_cell': pre, 'post_cell': post, 'pre_segment': pre_segment.sec.name(), 'pre_x': float(pre_segment.x), 'post_segment': post_segment.sec.name(), 'post_x': float(post_segment.x), 'pre_distance_um': pre_distance, 'post_distance_um': post_distance})
    # Tie every native recorder to its owning thread. Record each fixed step
    # and downsample afterward, matching the validated current package.
    recording_owners = {}
    soma = {}
    for name, cell in cells.items():
        owner = h.IClamp(cell.soma[0](.5))
        owner.amp = owner.dur = 0.
        recording_owners[name] = owner
        soma[name] = h.Vector().record(owner, cell.soma[0](.5)._ref_v)
    stimulus = h.IClamp(cells[selected[0]].soma[0](.5))
    stimulus.amp, stimulus.delay, stimulus.dur = args.amplitude, 0., 1.
    # Global t has no section: use the clamp to identify its recording thread.
    t = h.Vector().record(stimulus, h._ref_t)
    pc.set_maxstep(1.)
    setup_seconds = time.perf_counter() - beginning
    print(f'Starting fixed-step simulation: {len(cells)} cells, {len(connections)} synapses, {args.duration:g} ms, {args.threads} thread(s)', flush=True)
    integration_start = time.perf_counter()
    integration_cpu_start = time.process_time()
    h.finitialize(-65)
    initialization_seconds = time.perf_counter() - integration_start
    stop = 0.
    while stop < args.duration:
        stop = min(args.duration, stop + 5.)
        pc.psolve(stop)
        print(f'Simulated {stop:g}/{args.duration:g} ms; integration {time.perf_counter() - integration_start:.1f} s', flush=True)
    integration_seconds = time.perf_counter() - integration_start
    integration_cpu_seconds = time.process_time() - integration_cpu_start
    thread_segments = [0] * int(pc.nthread())
    for cell in cells.values():
        for section in cell.sections:
            thread_segments[int(h.l23_thread_id(sec=section))] += section.nseg
    times = np.asarray(t)[::4].copy()
    voltages = np.asarray([np.asarray(soma[name])[::4] for name in selected])
    if not np.isfinite(voltages).all() or not np.isfinite(times).all():
        raise RuntimeError('Simulation returned non-finite recordings')
    if abs(float(h.t) - args.duration) > max(1e-7, args.duration * 1e-10):
        raise RuntimeError('Simulation did not reach its requested stop time')
    if thread_segments != partition_loads:
        raise RuntimeError('Actual native partition differs from planned loads')
    if not all(np.isfinite(segment.v) for cell in cells.values() for segment in cell.segments):
        raise RuntimeError('Non-finite final compartment voltage')
    np.savez_compressed(args.output / 'recordings.npz', time_ms=times, soma_mv=voltages, cell_ids=np.asarray(selected))
    with (args.output / 'soma_voltage.csv').open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['time_ms'] + selected)
        writer.writerows(zip(times, *voltages))
    if mappings:
        with (args.output / 'synapse_locations.csv').open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(mappings[0]))
            writer.writeheader()
            writer.writerows(mappings)
    spikes = {name: np.asarray(vector).tolist() for name, vector in soma_events.items()}
    (args.output / 'soma_spikes.json').write_text(json.dumps(spikes, indent=2) + '\n')
    distances = [value for row in mappings for value in [row['pre_distance_um'], row['post_distance_um']]]
    rss = next(int(line.split()[1]) * 1024 for line in Path('/proc/self/status').read_text().splitlines() if line.startswith('VmRSS:'))
    summary = {
        'scope': 'Active-channel demonstration on downloaded SWC geometry, with original contact coordinates; not a calibrated reconstruction of the missing custom model',
        'layer': 'L23', 'available_morphologies': len(files), 'valid_morphologies': len(valid), 'excluded_invalid_morphologies': skipped,
        'all_table_contacts': len(all_contacts), 'contacts_between_supplied_l23_cells': len(within_layer), 'contacts_between_valid_l23_cells': len(valid_contacts),
        'simulated_cells': len(cells), 'selected_cell_ids': selected, 'simulated_contacts': len(connections),
        'directed_pairs': len({(r['pre_pt_root_id'], r['post_pt_root_id']) for r in contacts}),
        'sections': sum(len(c.sections) for c in cells.values()), 'segments': sum(len(c.segments) for c in cells.values()),
        'biophysics_profile': 'ModelDB 241165 / ' + profile_name, 'profile_source': profile_url,
        'original_wrapper_discretization': args.original_family,
        'source_hashes': {str(path.relative_to(args.data)): hashlib.sha256(path.read_bytes()).hexdigest() for path in [connectivity, *[valid[root] for root in selected]]},
        'recording_method': 'Owning-thread point-process anchors, fixed-step native vectors, downsampled by four to 0.1 ms',
        'temperature_c': 34, 'dt_ms': .025, 'record_dt_ms': .1, 'duration_ms': args.duration, 'actual_final_time_ms': float(h.t), 'threads': int(pc.nthread()),
        'hostname': platform.node(), 'neuron_version': neuron.__version__, 'segments_per_thread': thread_segments,
        'zero_conductance_pruning_enabled': True, 'removed_segment_instances': dict(sum((Counter(c.removed_zero_channels) for c in cells.values()), Counter())),
        'stimulus': {'cell_id': selected[0], 'amplitude_na': args.amplitude, 'delay_ms': 0, 'duration_ms': 1},
        'synapse_parameters': {'mechanism': 'Exp2Syn', 'tau1_ms': .5, 'tau2_ms': 3, 'reversal_mv': 0, 'weight_us': args.synaptic_weight, 'delay_ms': 1, 'source_threshold_mv': 0},
        'soma_spike_count': sum(len(v) for v in spikes.values()), 'cells_with_soma_spikes': sum(bool(v) for v in spikes.values()),
        'estimated_delivered_synaptic_events': sum(int(np.count_nonzero(np.asarray(v) + 1 <= args.duration + 1e-8)) for v in connection_events),
        'contact_mapping_distance_um': {'median': float(np.median(distances)), 'p95': float(np.percentile(distances, 95)), 'maximum': float(max(distances)), 'endpoints_over_20um': sum(d > 20 for d in distances)} if distances else None,
        'setup_seconds': setup_seconds, 'integration_seconds': integration_seconds, 'initialization_seconds': initialization_seconds, 'solve_seconds': integration_seconds - initialization_seconds, 'integration_cpu_seconds': integration_cpu_seconds, 'resident_ram_bytes_before_plotting': rss, 'recorded_samples': len(times), 'all_recordings_finite': True,
        'assumptions': ['The same published ' + profile_name + ' profile is assigned to every selected morphology; cell-type-specific fits are not provided.', 'SWC type 3 is basal dendrite; type 4, if present, is apical. No missing apical labels are inferred.', 'Contacts are modeled as excitatory using the original script defaults; the table provides no excitatory/inhibitory labels or calibrated weights.', 'Synapses attach to the nearest computational segment; mapping distances are reported.', 'No extracellular stimulation or myelination is added.'],
    }
    (args.output / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(12, 5))
    if len(selected) <= 8:
        for index, name in enumerate(selected):
            ax.plot(times, voltages[index], label=name, linewidth=1)
        ax.legend(fontsize=8)
    else:
        ax.fill_between(times, voltages.min(axis=0), voltages.max(axis=0), alpha=.2, label='Population min–max')
        ax.plot(times, np.mean(voltages, axis=0), label='Population mean')
        ax.plot(times, voltages[0], label='Stimulated cell', linewidth=1)
        ax.legend()
    ax.set(xlabel='Time (ms)', ylabel='Soma voltage (mV)', title=f'L23 demonstration: {len(cells)} cells, {len(connections)} contacts; {args.amplitude:g} nA pulse')
    fig.tight_layout()
    fig.savefig(args.output / 'soma_voltage.png', dpi=150)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(12, 5))
    for index, name in enumerate(selected):
        if spikes[name]:
            ax.vlines(spikes[name], index - .4, index + .4, linewidth=.8)
    ax.set(xlim=(0, args.duration), ylim=(-1, len(selected)), xlabel='Time (ms)', ylabel='Cell index (stimulated cell = 0)', title=f'L23 soma spikes: {summary["soma_spike_count"]} events; weight {args.synaptic_weight:g} uS')
    fig.tight_layout()
    fig.savefig(args.output / 'spike_raster.png', dpi=150)
    plt.close(fig)
    print(json.dumps({k: summary[k] for k in ['simulated_cells', 'simulated_contacts', 'segments', 'soma_spike_count', 'estimated_delivered_synaptic_events', 'setup_seconds', 'integration_seconds', 'resident_ram_bytes_before_plotting']}, indent=2), flush=True)


if __name__ == '__main__':
    main()
