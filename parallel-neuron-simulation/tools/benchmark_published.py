"""Profile the published example and compare explicitly selected single-thread solvers."""

import argparse
import cProfile
import io
import json
import os
from pathlib import Path
import pstats
import runpy
import time


ROOT = Path(__file__).resolve().parents[1]
PUBLISHED = ROOT / 'vendor' / 'aberra2018'


def compare(candidate_t, candidate_v, reference_t, reference_v):
    import numpy as np
    common = (reference_t >= candidate_t[0]) & (reference_t <= candidate_t[-1] + 1e-7)
    times = reference_t[common]
    reference = reference_v[common]
    difference = np.interp(times, candidate_t, candidate_v) - reference
    index = int(np.argmax(np.abs(difference)))
    return {
        'compared_samples': len(times),
        'max_abs_soma_difference_mv': float(np.abs(difference).max()),
        'rms_soma_difference_mv': float(np.sqrt(np.mean(difference ** 2))),
        'time_of_max_difference_ms': float(times[index]),
        'final_soma_difference_mv': float(difference[-1]),
        'candidate_soma_upward_zero_crossings': int(((candidate_v[:-1] < 0) & (candidate_v[1:] >= 0)).sum()),
        'reference_soma_upward_zero_crossings': int(((reference_v[:-1] < 0) & (reference_v[1:] >= 0)).sum()),
        'interpretation': 'Trace differences against the selected reference; not a certified error bound or a validation of spike/threshold protocols',
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--solver', choices=('baseline', 'prune-zero', 'cvode', 'cvode-tight'), default='baseline')
    args = parser.parse_args()
    output = ROOT / 'outputs' / 'benchmarks'
    output.mkdir(parents=True, exist_ok=True)
    os.environ['NEURON_MODULE_OPTIONS'] = os.environ.get('NEURON_MODULE_OPTIONS', '') + ' -nogui'
    os.environ['OMP_NUM_THREADS'] = '1'
    os.environ['OPENBLAS_NUM_THREADS'] = '1'
    import numpy as np
    import neuron
    from neuron import h

    total_start = time.perf_counter()
    neuron.load_mechanisms(str(PUBLISHED / 'mechanisms'))
    h.ParallelContext().nthread(1)
    os.chdir(PUBLISHED / 'cells' / 'L4_LBC_cACint209_1')
    namespace = runpy.run_path('run.py', run_name='published_benchmark')
    namespace['init_simulation']()
    setup_start = time.perf_counter()
    cell = namespace['create_cell'](add_synapses=False)
    removed = {}
    if args.solver == 'prune-zero':
        if __package__:
            from .simulation_optimizations import prune_zero_conductance_channels
        else:
            from simulation_optimizations import prune_zero_conductance_channels
        removed = prune_zero_conductance_channels(cell.all)
    stimuli = namespace['create_stimuli'](cell, 1)
    recordings = namespace['create_recordings'](cell)
    setup_seconds = time.perf_counter() - setup_start
    cvode = h.CVode()
    h.tstop = 3000
    if args.solver in ('baseline', 'prune-zero'):
        cvode.active(0)
    else:
        tolerance = 1e-6 if args.solver == 'cvode' else 1e-8
        cvode.atol(tolerance)
        cvode.rtol(tolerance)
        cvode.atolscale('cai', 1e-3)
        cvode.active(1)
    sections = list(cell.all)
    metadata = {
        'solver': args.solver,
        'neuron_version': neuron.__version__,
        'threads': int(h.ParallelContext().nthread()),
        'mpi_ranks': int(h.ParallelContext().nhost()),
        'cache_efficient': bool(cvode.cache_efficient()),
        'sections': len(sections),
        'segments': sum(s.nseg for s in sections),
        'duration_ms': float(h.tstop),
        'fixed_dt_ms': float(h.dt) if args.solver in ('baseline', 'prune-zero') else None,
        'atol': float(cvode.atol()) if cvode.active() else None,
        'rtol': float(cvode.rtol()) if cvode.active() else None,
        'cai_atolscale': float(cvode.atolscale('cai')) if cvode.active() else None,
        'setup_seconds': setup_seconds,
        'removed_zero_conductance_segment_instances': removed,
        'source_helpers': 'Unchanged published helpers, step 1; identical spatial discretization. prune-zero removes only static zero-conductance Nap_Et2/K_Tst instances; CVode modes change the integration method.',
    }
    print(json.dumps({'starting': metadata}), flush=True)
    profiler = cProfile.Profile()
    start, cpu_start = time.perf_counter(), time.process_time()
    def integrate():
        # HocObject calls bypass cProfile's builtin-call hooks. A Python
        # boundary makes their inclusive elapsed time visible in the profile.
        h.run()

    profiler.runcall(integrate)
    metadata['integration_seconds'] = time.perf_counter() - start
    metadata['integration_cpu_seconds'] = time.process_time() - cpu_start
    metadata['actual_final_time_ms'] = float(h.t)
    post_start = time.perf_counter()
    t = np.asarray(recordings['time']).copy()
    v = np.asarray(recordings['soma(0.5)']).copy()
    if len(t) < 30000 or not np.isfinite(v).all():
        raise RuntimeError('Invalid benchmark recording')
    np.savez(output / f'{args.solver}.npz', time_ms=t, soma_mv=v)
    metadata['samples'] = len(t)
    metadata['array_copy_and_binary_export_seconds'] = time.perf_counter() - post_start
    metadata['total_seconds_excluding_imports'] = time.perf_counter() - total_start
    stream = io.StringIO()
    pstats.Stats(profiler, stream=stream).sort_stats('cumulative').print_stats(25)
    (output / f'{args.solver}_profile.txt').write_text(stream.getvalue())
    reference_path = output / 'baseline.npz'
    if args.solver != 'baseline' and reference_path.exists():
        with np.load(reference_path) as reference:
            metadata['comparison_to_fixed_step'] = compare(t, v, reference['time_ms'], reference['soma_mv'])
        baseline = json.loads((output / 'baseline.json').read_text())
        metadata['integration_speedup_vs_baseline'] = baseline['integration_seconds'] / metadata['integration_seconds']
    if args.solver == 'baseline':
        reference = np.loadtxt(ROOT / 'outputs/published_model/soma_voltage_step1.dat')
        metadata['comparison_to_original_run'] = compare(t, v, reference[:, 0], reference[:, 1])
    (output / f'{args.solver}.json').write_text(json.dumps(metadata, indent=2) + '\n')
    print(json.dumps(metadata, indent=2), flush=True)
    # Keep the cell and both current clamps alive until recordings have been exported.
    assert cell is not None and len(stimuli) == 2


if __name__ == '__main__':
    main()
