"""Run the published current-clamp example with validated fixed-step optimizations."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import runpy


ROOT = Path(__file__).resolve().parents[1]
PUBLISHED = ROOT / 'vendor' / 'aberra2018'
CELL_NAME = 'L4_LBC_cACint209_1'


def summarize(step, output, optimization=None):
    import numpy as np
    import matplotlib.pyplot as plt

    original = output / f'soma_voltage_step{step}.dat'
    values = np.loadtxt(original)
    if values.shape not in ((30000, 2), (30001, 2)) or not np.isfinite(values).all():
        raise RuntimeError('Published run did not produce the expected finite recording')
    times, voltage = values.T
    # Fixed-step floating-point timing can omit a recording exactly at tstop.
    if not np.isclose(times[0], 0) or not 2999.89999 <= times[-1] <= 3000.00001:
        raise RuntimeError('Published recording has unexpected simulation times')
    if not np.allclose(np.diff(times), .1, atol=1e-8):
        raise RuntimeError('Published recording has unexpected sample intervals')
    np.savetxt(output / f'soma_voltage_step{step}.csv', values, delimiter=',', header='time_ms,soma_voltage_mv', comments='')
    fig, ax = plt.subplots(figsize=(12, 4))
    ax.axvspan(700, 2700, color='#dddddd', alpha=.6, label='Step current interval')
    ax.plot(times, voltage, linewidth=.8, color='#225fa8', label='Soma voltage')
    ax.set(xlabel='Time (ms)', ylabel='Membrane potential (mV)', title=f'ModelDB 241165 bundled example · {CELL_NAME} · step {step}')
    ax.legend(loc='lower right')
    fig.tight_layout()
    fig.savefig(output / f'soma_voltage_step{step}.png', dpi=150)
    plt.close(fig)
    summary = {
        'source': 'https://modeldb.science/241165',
        'commit': '93e01d4a927ab0cc5ca137329b86925777e63f02',
        'cell': CELL_NAME,
        'entry_point': 'Published run.py helpers with selected runtime optimizations: init_simulation(); run_step(step, plot_traces=False)',
        'optimization': optimization or {'prune_zero_conductance': False},
        'solver': 'fixed-step',
        'scope': 'Bundled Blue Brain baseline example; not the adult myelinated model, a paper figure reproduction, or the downloaded three-cell SWC network',
        'step': step,
        'duration_ms': 3000,
        'dt_ms': .025,
        'record_dt_ms': .1,
        'temperature_c': 34,
        'samples': len(times),
        'last_recorded_time_ms': float(times[-1]),
        'step_current_na': [.0630936, .0683514, .0736092][step - 1],
        'step_delay_ms': 700,
        'step_duration_ms': 2000,
        'holding_current_na': -.031738,
        'minimum_soma_mv': float(voltage.min()),
        'maximum_soma_mv': float(voltage.max()),
        'final_soma_mv': float(voltage[-1]),
        'soma_upward_zero_crossings': int(((voltage[:-1] < 0) & (voltage[1:] >= 0)).sum()),
        'all_values_finite': True,
    }
    (output / f'summary_step{step}.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--step', type=int, choices=(1, 2, 3), default=1)
    parser.add_argument('--no-prune-zero', action='store_true', help='Retain zero-conductance mechanisms to reproduce the unoptimized baseline')
    args = parser.parse_args()
    output = ROOT / 'outputs' / 'published_model'
    output.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((PUBLISHED / 'SOURCE_MANIFEST.json').read_text())
    for entry in manifest['files']:
        if hashlib.sha256((PUBLISHED / entry['path']).read_bytes()).hexdigest() != entry['sha256']:
            raise RuntimeError(f'Published file modified: {entry["path"]}')
    os.environ['NEURON_MODULE_OPTIONS'] = os.environ.get('NEURON_MODULE_OPTIONS', '') + ' -nogui'
    os.environ['MPLBACKEND'] = 'Agg'
    os.environ.setdefault('MPLCONFIGDIR', str(output / 'matplotlib-cache'))
    import neuron
    from neuron import h
    h.ParallelContext().nthread(1)
    if not neuron.load_mechanisms(str(PUBLISHED / 'mechanisms')):
        raise RuntimeError('Compile the published mechanisms before running this example')
    directory = PUBLISHED / 'cells' / CELL_NAME
    recordings = directory / 'python_recordings'
    if not recordings.exists():
        recordings.symlink_to(output, target_is_directory=True)
    if recordings.resolve() != output.resolve():
        raise RuntimeError('Published example recording directory points to another location')
    os.chdir(directory)
    namespace = runpy.run_path(str(directory / 'run.py'), run_name='published_cell_example')
    optimization = {'prune_zero_conductance': not args.no_prune_zero, 'removed_segment_instances': {}}
    if not args.no_prune_zero:
        if __package__:
            from .simulation_optimizations import prune_zero_conductance_channels
        else:
            from simulation_optimizations import prune_zero_conductance_channels
        original_create_cell = namespace['create_cell']

        def create_optimized_cell(add_synapses=True):
            cell = original_create_cell(add_synapses=add_synapses)
            optimization['removed_segment_instances'] = prune_zero_conductance_channels(cell.all)
            return cell

        # runpy returns a copy of the namespace; update the function's actual globals.
        namespace['run_step'].__globals__['create_cell'] = create_optimized_cell
    namespace['init_simulation']()
    namespace['run_step'](args.step, plot_traces=False)
    summarize(args.step, output, optimization)


if __name__ == '__main__':
    main()
