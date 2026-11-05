"""Run a scaled, published mouse M1 network from ModelDB 260015."""

import argparse
from collections import Counter
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import runpy
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'vendor/m1_260015'


def connection_totals(cells):
    """Count actual receptor objects and separately count connection-site keys."""
    receptors = Counter()
    netcons = 0
    sites = set()
    pairs = set()
    digest = hashlib.sha256()
    for cell in cells:
        for section in getattr(cell, 'secs', {}).values():
            for synapse in section.get('synMechs', []):
                if synapse.get('hObj') is not None:
                    receptors[synapse['label']] += 1
        for connection in cell.conns:
            if connection.get('hObj') is None:
                continue
            netcons += 1
            pre = connection['preGid']
            pairs.add((pre, cell.gid))
            sites.add((pre, cell.gid, connection.get('sec'), connection.get('loc')))
            digest.update((repr((pre, cell.gid, connection.get('sec'), connection.get('loc'), connection.get('synMech'), connection.get('weight'), connection.get('delay')))+'\n').encode())
    return {
        'synaptic_point_processes': sum(receptors.values()),
        'synaptic_netcons': netcons,
        'unique_connection_site_keys': len(sites),
        'directed_cell_pairs': len(pairs),
        'receptors': dict(receptors),
        'connectivity_sha256': digest.hexdigest(),
    }


def rss_bytes():
    return next(int(line.split()[1]) * 1024 for line in Path('/proc/self/status').read_text().splitlines() if line.startswith('VmRSS:'))


def non_thread_safe_mechanisms(source_dir):
    """Read the compiler's vectorization flags rather than changing MOD files."""
    unsafe = set()
    for path in (source_dir / 'x86_64').glob('*.cpp'):
        text = path.read_text()
        if '#define NRN_VECTORIZED 0' in text:
            match = re.search(r'_mechanism\[\]\s*=\s*\{\s*"[^"]*",\s*"([^"]+)"', text)
            if match:
                unsafe.add(match.group(1))
    return unsafe


def scale_populations(populations, scale):
    """Apply the same count multiplier to density and fixed-count populations."""
    for population in populations.values():
        if 'density' in population:
            population['density'] *= scale
        elif 'numCells' in population:
            population['numCells'] = int(population['numCells'] * scale)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--scale', type=float, default=.02, help='Published population-count multiplier')
    parser.add_argument('--threads', type=int, default=1)
    parser.add_argument('--duration', type=float, default=25., help='Simulated milliseconds')
    parser.add_argument('--output', type=Path, default=ROOT / 'outputs/mouse_m1_pilot')
    parser.add_argument('--build-only', action='store_true')
    parser.add_argument('--minimum-synapses', type=int, default=0)
    parser.add_argument('--minimum-contact-sites', type=int, default=0)
    parser.add_argument('--record-cells', type=int, default=32)
    parser.add_argument('--evoked', action='store_true', help="Enable the publication's optional 0.5 nA IT5B current input")
    parser.add_argument('--save-connections', action='store_true', help='Export the generated receptor-level connection table')
    args = parser.parse_args()
    if args.scale <= 0 or args.threads < 1 or args.duration <= 0 or args.record_cells < 1:
        parser.error('Scale, threads, duration, and recording count must be positive')
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    os.environ['NEURON_MODULE_OPTIONS'] = os.environ.get('NEURON_MODULE_OPTIONS', '') + ' -nogui'
    os.environ['MPLBACKEND'] = 'Agg'
    os.environ['MPLCONFIGDIR'] = str(output / 'matplotlib-cache')
    os.environ['OMP_NUM_THREADS'] = '1'
    os.environ['OPENBLAS_NUM_THREADS'] = '1'
    import numpy as np
    import neuron
    from neuron import h
    import netpyne
    from netpyne import sim
    from netpyne.cell.compartCell import CompartCell

    started = time.perf_counter()
    manifest = json.loads((SOURCE / 'SOURCE_MANIFEST.json').read_text())
    for item in manifest['files']:
        if hashlib.sha256((SOURCE / item['path']).read_bytes()).hexdigest() != item['sha256']:
            raise RuntimeError(f"Published source changed: {item['path']}")
    source_dir = SOURCE / 'sim'
    neuron.load_mechanisms(str(source_dir))
    os.chdir(source_dir)
    sys.path.insert(0, str(source_dir))
    # The unchanged netParams.py intentionally reads cfg from __main__.
    global cfg
    cfg = runpy.run_path(str(source_dir / 'cfg.py'))['cfg']
    cfg.scale = 1.
    cfg.duration = args.duration
    cfg.recordLFP = []
    cfg.recordDipole = False
    cfg.recordDipoles = False
    cfg.recordTraces = {}
    cfg.recordCells = []
    cfg.recordTime = False
    cfg.analysis = {}
    cfg.savePickle = False
    cfg.saveJson = False
    cfg.printRunTime = False
    cfg.verbose = False
    if args.evoked:
        cfg.addIClamp = 1
    # The published parameter modification starts at 500 ms; short benchmarks
    # end before it, and do not invoke the publication's interval callback.
    if args.duration >= cfg.modifyMechs['startTime']:
        parser.error('This benchmark must end before the published 500 ms mechanism-change protocol')
    net_params = runpy.run_path(str(source_dir / 'netParams.py'))['netParams']
    # NetPyNE 1.1.1 applies params.scale to fixed counts but not density-based
    # populations. Scale both forms explicitly and keep its scale at one.
    net_params.scale = 1.
    scale_populations(net_params.popParams, args.scale)
    sim.initialize(simConfig=cfg, netParams=net_params)
    # Build on thread zero, including the published non-thread-safe SOM channel.
    sim.pc.nthread(1)
    sim.net.createPops()
    sim.net.createCells()
    print(f'Cells created in {time.perf_counter()-started:.2f} s', flush=True)
    sim.net.connectCells()
    sim.net.addStims()
    cells = [cell for cell in sim.net.cells if isinstance(cell, CompartCell)]
    totals = connection_totals(cells)
    hoc_receptors = {name: int(h.List(name).count()) for name in ['MyExp2SynBB', 'MyExp2SynNMDABB']}
    if sum(hoc_receptors.values()) != totals['synaptic_point_processes']:
        raise RuntimeError('Python and NEURON synaptic object counts disagree')
    sections = [section['hObj'] for cell in cells for section in cell.secs.values()]
    segments = sum(section.nseg for section in sections)
    partitions = [h.SectionList() for _ in range(args.threads)]
    loads = [0] * args.threads
    unsafe = non_thread_safe_mechanisms(source_dir)
    pinned = []
    ordered = sorted(cells, key=lambda cell: (-sum(s['hObj'].nseg for s in cell.secs.values()), cell.gid))
    for cell in ordered:
        restricted = {mechanism for section in cell.secs.values() for mechanism in section.get('mechs', {})} & unsafe
        thread = 0 if restricted else min(range(args.threads), key=lambda i: loads[i])
        if restricted:
            pinned.append({'gid': cell.gid, 'mechanisms': sorted(restricted)})
        for section in cell.secs.values():
            sec = section['hObj']
            if sec.parentseg() is None:
                partitions[thread].append(sec=sec)
        loads[thread] += sum(s['hObj'].nseg for s in cell.secs.values())
    sim.pc.nthread(args.threads)
    for thread, roots in enumerate(partitions):
        sim.pc.partition(thread, roots)
    h('objref m1_thread_context')
    h.m1_thread_context = sim.pc
    h('func m1_thread_id() { return m1_thread_context.sec_in_thread() }')
    selected = sorted(cells, key=lambda cell: cell.gid)
    if len(selected) > args.record_cells:
        indices = np.linspace(0, len(selected)-1, args.record_cells, dtype=int)
        selected = [selected[i] for i in indices]
    records, anchors = [], []
    for cell in selected:
        sec = cell.secs['soma']['hObj']
        anchor = next((synapse['hObj'] for section in cell.secs.values() for synapse in section.get('synMechs', []) if synapse.get('hObj') is not None), None)
        if anchor is None:
            anchor = h.IClamp(sec(.5))
            anchor.amp = anchor.dur = 0
        anchors.append(anchor)
        # Record each fixed step, then decimate by integer indices on export.
        # No timed interpolation callbacks or shared global-time recorder.
        records.append(h.Vector().record(anchor, sec(.5)._ref_v))
    spike_times, spike_gids = h.Vector(), h.Vector()
    sim.pc.spike_record(-1, spike_times, spike_gids)
    summary = {
        'scope': 'Scaled published mouse M1 model; generated connectivity, not a measured synapse list',
        'modeldb_url': manifest['modeldb_url'], 'source_commit': manifest['commit'],
        'hostname': platform.node(), 'neuron_version': neuron.__version__, 'netpyne_version': netpyne.__version__,
        'scale': args.scale, 'threads': args.threads, 'mpi_ranks': int(sim.pc.nhost()),
        'optional_current_input_enabled': args.evoked,
        'optional_current_input_parameters': dict(cfg.IClamp1) if args.evoked else None,
        'biophysical_cells': len(cells), 'artificial_input_cells': len(sim.net.cells)-len(cells),
        'populations': dict(Counter(cell.tags['pop'] for cell in sim.net.cells)),
        'sections': len(sections), 'segments': segments, **totals,
        'cells_restricted_to_thread_zero': pinned,
        'hoc_synaptic_object_counts': hoc_receptors,
        'dt_ms': float(cfg.dt), 'duration_ms': args.duration,
        'recorded_cell_ids': [cell.gid for cell in selected],
        'recording_method': 'Every fixed step, with explicit cell point-process anchors; exported every four steps',
        'construction_seconds': time.perf_counter()-started, 'resident_ram_bytes_after_construction': rss_bytes(),
        'benchmark_changes': ['Population counts scaled; published cell models, synaptic parameters, connection rules, and seeds retained.', 'LFP, dipole, analysis plotting and publication saving disabled.', 'Run ends before the 500 ms conductance-modification protocol.'],
    }
    (output / 'construction.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps({k: summary[k] for k in ['biophysical_cells','artificial_input_cells','synaptic_point_processes','synaptic_netcons','unique_connection_site_keys','segments','construction_seconds']}, indent=2), flush=True)
    if totals['synaptic_point_processes'] < args.minimum_synapses:
        raise RuntimeError('Constructed fewer synaptic point processes than required')
    if totals['unique_connection_site_keys'] < args.minimum_contact_sites:
        raise RuntimeError('Constructed fewer distinct modeled contact-site keys than required')
    if args.save_connections:
        with gzip.open(output/'connections.csv.gz', 'wt', compresslevel=1, newline='') as stream:
            writer = csv.writer(stream)
            writer.writerow(['pre_gid','post_gid','synaptic_mechanism','post_section','post_location','weight','delay_ms'])
            for cell in cells:
                for connection in cell.conns:
                    if connection.get('hObj') is not None:
                        writer.writerow([connection['preGid'],cell.gid,connection.get('synMech'),connection.get('sec'),connection.get('loc'),connection.get('weight'),connection.get('delay')])
        with (output/'cells.csv').open('w', newline='') as stream:
            writer = csv.writer(stream)
            writer.writerow(['gid','population','cell_type','cell_model','x_um','y_um','z_um'])
            for cell in sim.net.cells:
                writer.writerow([cell.gid]+[cell.tags.get(key) for key in ['pop','cellType','cellModel','x','y','z']])
    if args.build_only:
        return
    init_start = time.perf_counter()
    sim.preRun()
    h.finitialize(float(cfg.hParams['v_init']))
    summary['initialization_seconds'] = time.perf_counter()-init_start
    actual_loads = [0] * args.threads
    for section in sections:
        actual_loads[int(h.m1_thread_id(sec=section))] += section.nseg
    if actual_loads != loads:
        raise RuntimeError('Actual thread partitions disagree with assigned cells')
    summary['segments_per_thread'] = actual_loads
    solve_start, cpu_start = time.perf_counter(), time.process_time()
    stop = 0.
    while stop < args.duration:
        stop = min(args.duration, stop+5.)
        sim.pc.psolve(stop)
        print(f'Simulated {stop:g}/{args.duration:g} ms; solve {time.perf_counter()-solve_start:.2f} s', flush=True)
    summary['solve_seconds'] = time.perf_counter()-solve_start
    summary['solve_cpu_seconds'] = time.process_time()-cpu_start
    summary['actual_final_time_ms'] = float(h.t)
    arrays = [np.asarray(vector).copy() for vector in records]
    if not arrays or len({len(a) for a in arrays}) != 1:
        raise RuntimeError('Inconsistent fixed-step recording lengths')
    voltages = np.asarray(arrays)[:, ::4]
    times = np.arange(voltages.shape[1]) * float(cfg.dt)*4
    if not np.isfinite(voltages).all() or abs(float(h.t)-args.duration) > 1e-7:
        raise RuntimeError('Nonfinite voltage or incorrect final time')
    event_times = np.asarray(spike_times).copy()
    event_gids = np.asarray(spike_gids, dtype=np.int64).copy()
    order = np.lexsort((event_gids, event_times))
    event_times, event_gids = event_times[order], event_gids[order]
    real_ids = {cell.gid for cell in cells}
    real_mask = np.asarray([int(gid) in real_ids for gid in event_gids], dtype=bool)
    summary['biophysical_cell_spikes'] = int(real_mask.sum())
    summary['artificial_input_spikes'] = int((~real_mask).sum())
    summary['recorded_samples'] = len(times)
    summary['all_recordings_finite'] = True
    summary['resident_ram_bytes_after_simulation'] = rss_bytes()
    np.savez_compressed(output/'recordings.npz',time_ms=times,soma_mv=voltages,cell_ids=np.asarray(summary['recorded_cell_ids']),spike_times_ms=event_times,spike_cell_ids=event_gids)
    (output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps({k:summary[k] for k in ['threads','synaptic_point_processes','initialization_seconds','solve_seconds','biophysical_cell_spikes','artificial_input_spikes','resident_ram_bytes_after_simulation']},indent=2),flush=True)


if __name__ == '__main__':
    main()
