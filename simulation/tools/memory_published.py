"""Measure four published cells and calculate voltage-recording storage requirements."""

import argparse
import json
import math
import os
from pathlib import Path
import resource
import runpy


ROOT = Path(__file__).resolve().parents[1]
PUBLISHED = ROOT / 'vendor' / 'aberra2018'


def memory_snapshot():
    values = {}
    for line in Path('/proc/self/status').read_text().splitlines():
        key, _, rest = line.partition(':')
        if key in ('VmRSS', 'VmHWM', 'VmSize'):
            values[key + '_bytes'] = int(rest.split()[0]) * 1024
    values['process_peak_rss_bytes'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024
    return values


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cells', type=int, default=4)
    parser.add_argument('--duration', type=float, default=3000., help='Target recording duration in ms')
    args = parser.parse_args()
    if args.cells < 1 or args.duration <= 0:
        parser.error('Cell count and target duration must be positive')
    output = ROOT / 'outputs' / 'memory'
    output.mkdir(parents=True, exist_ok=True)
    os.environ['NEURON_MODULE_OPTIONS'] = os.environ.get('NEURON_MODULE_OPTIONS', '') + ' -nogui'
    import neuron
    from neuron import h
    if __package__:
        from .simulation_optimizations import prune_zero_conductance_channels
    else:
        from simulation_optimizations import prune_zero_conductance_channels
    neuron.load_mechanisms(str(PUBLISHED / 'mechanisms'))
    h.ParallelContext().nthread(1)
    os.chdir(PUBLISHED / 'cells' / 'L4_LBC_cACint209_1')
    namespace = runpy.run_path('run.py', run_name='published_memory_probe')
    namespace['init_simulation']()
    snapshots = {'runtime_loaded': memory_snapshot()}
    cells, stimuli, pruned = [], [], []
    for index in range(args.cells):
        cell = namespace['create_cell'](add_synapses=False)
        cells.append(cell)
        pruned.append(prune_zero_conductance_channels(cell.all))
        stimuli.extend(namespace['create_stimuli'](cell, 1))
        assert len(list(h.allsec())) == sum(len(list(c.all)) for c in cells)
        snapshots[f'after_cell_{index + 1}'] = memory_snapshot()
    segments = sum(s.nseg for c in cells for s in c.all)
    samples = math.ceil(args.duration / .1) + 1
    # Preallocate only shared time and one soma vector per cell, not the much
    # larger all-compartment animation buffers. Initialization retains capacity.
    time_vector = h.Vector(samples).record(h._ref_t, .1)
    soma_vectors = [h.Vector(samples).record(c.soma[0](.5)._ref_v, .1) for c in cells]
    snapshots['soma_buffers_allocated'] = memory_snapshot()
    h.cvode_active(0)
    h.finitialize(-65)
    h.continuerun(min(25., args.duration))
    snapshots['after_short_probe'] = memory_snapshot()
    capacities = [int(v.buffer_size()) for v in [time_vector, *soma_vectors]]
    layouts = []
    for duration in sorted({25., args.duration}):
        for interval in [.1, .025]:
            count = math.ceil(duration / interval) + 1
            raw_all = (segments + 1) * count * 8
            # The supplied animator creates per-vector float32 copies and then
            # stacks them, temporarily holding two float32 matrices as well as
            # the original double-precision NEURON recordings.
            peak_arrays = raw_all + 2 * segments * count * 4
            layouts.append({'duration_ms': duration, 'record_dt_ms': interval, 'samples': count, 'soma_plus_shared_time_bytes': (args.cells + 1) * count * 8, 'all_segments_plus_time_raw_bytes': raw_all, 'animator_peak_array_payload_bytes': peak_arrays, 'excludes': 'Vector spare capacity, object overhead, plotting geometry and runtime RSS'})
    report = {
        'scope': 'Four copies of published L4_LBC_cACint209_1 by default; independent cells, no synapses or animation',
        'cells': len(cells),
        'sections': sum(len(list(c.all)) for c in cells),
        'segments': segments,
        'neuron_version': neuron.__version__,
        'threads': int(h.ParallelContext().nthread()),
        'cvode_active': bool(h.cvode_active()),
        'target_duration_ms': args.duration,
        'probe_duration_ms': min(25., args.duration),
        'soma_record_dt_ms': .1,
        'preallocated_soma_recording_capacity_bytes': sum(capacities) * 8,
        'pruned_segment_instances_by_cell': pruned,
        'snapshots': snapshots,
        'recording_storage_calculations': layouts,
        'method': 'Linux process RSS measured during setup and a short fixed-step probe. Large all-compartment recording payloads are calculated, not allocated or measured. Target duration is not the duration simulated by this probe.',
    }
    (output / f'published_{args.cells}_cells.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))
    assert len(stimuli) == args.cells * 2


if __name__ == '__main__':
    main()
