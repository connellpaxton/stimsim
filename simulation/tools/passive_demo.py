"""Passive cable demonstration on real SWC trees; not the original active NEURON model."""

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from analyze import DEFAULT_DATA, SAMPLE, load_swc


def factor_tree(diagonal, coupling, parent, order):
    """Hines elimination for a symmetric tree matrix with off-diagonal -coupling."""
    pivots = diagonal.copy()
    for i in order[:0:-1]:
        p = parent[i]
        pivots[p] -= coupling[i] ** 2 / pivots[i]
    if (pivots <= 0).any():
        raise ValueError('Cable matrix is not positive definite')
    return pivots


def solve_tree(rhs, pivots, coupling, parent, order):
    values = rhs.copy()
    for i in order[:0:-1]:
        values[parent[i]] += coupling[i] * values[i] / pivots[i]
    root = order[0]
    values[root] /= pivots[root]
    for i in order[1:]:
        values[i] = (values[i] + coupling[i] * values[parent[i]]) / pivots[i]
    return values


def simulate(path, duration=25., dt=.025, amplitude=1., delay=1., pulse_duration=1.):
    """Backward Euler: C/dt*(Vnew-Vold) + G*Vnew = I.

    Units: um geometry, nF capacitance, uS conductance, mV voltage, nA current.
    Each non-root node represents its parent edge as a cylindrical compartment.
    The root is a sphere. This is a deliberately simple, uncalibrated passive model.
    """
    if dt <= 0 or duration <= 0 or delay < 0 or pulse_duration < 0:
        raise ValueError('Invalid simulation timing')
    rows, parent, order = load_swc(path)
    root = order[0]
    if int(rows[root, 1]) != 1:
        raise ValueError('Passive demo requires a soma root')
    n = len(rows)
    edge = np.flatnonzero(parent >= 0)
    radius = rows[:, 5]
    lengths = np.zeros(n)
    lengths[edge] = np.linalg.norm(rows[edge, 2:5] - rows[parent[edge], 2:5], axis=1)
    # Preserve recorded radii. Only degenerate edge lengths receive a 0.01 um floor.
    clipped = int((lengths[edge] < .01).sum())
    lengths[edge] = np.maximum(lengths[edge], .01)
    area = 2 * np.pi * radius * lengths
    area[root] = 4 * np.pi * radius[root] ** 2
    capacitance = area * 1e-5  # cm=1 uF/cm2 -> nF/um2
    leak = area * 6.3e-7  # g_pas=6.3e-5 S/cm2 -> uS/um2
    coupling = np.zeros(n)
    # Series resistance between adjacent compartment centers; Ra=100 ohm cm.
    resistance = np.zeros(n)
    resistance[edge] = .5 * lengths[edge] / radius[edge] ** 2
    resistance[edge] += np.where(parent[edge] == root, 0, .5 * lengths[parent[edge]] / radius[parent[edge]] ** 2)
    coupling[edge] = np.pi / (100 * .01 * resistance[edge])
    diagonal = capacitance / dt + leak + coupling
    np.add.at(diagonal, parent[edge], coupling[edge])
    pivots = factor_tree(diagonal, coupling, parent, order)
    voltage = np.zeros(n)  # Relative to -65 mV resting potential.
    steps = round(duration / dt)
    times = np.arange(steps + 1) * dt
    traces = np.empty((steps + 1, 4))
    distal = int(np.argmax(np.linalg.norm(rows[:, 2:5] - rows[root, 2:5], axis=1)))
    traces[0] = (-65, -65, -65, -65)
    peak = np.zeros(n)
    for step in range(steps):
        t = step * dt
        rhs = capacitance / dt * voltage
        if delay <= t < delay + pulse_duration:
            rhs[root] += amplitude
        voltage = solve_tree(rhs, pivots, coupling, parent, order)
        peak = np.maximum(peak, voltage)
        traces[step + 1] = (voltage[root] - 65, voltage[distal] - 65, voltage.min() - 65, voltage.max() - 65)
    if not np.isfinite(traces).all():
        raise ValueError('Non-finite simulation result')
    return times, traces, {'cell_id': path.stem, 'compartments': n, 'area_um2': float(area.sum()), 'short_edges_floored': clipped, 'peak_soma_mv': float(traces[:, 0].max()), 'peak_distal_mv': float(traces[:, 1].max()), 'final_soma_mv': float(traces[-1, 0]), 'max_relative_voltage_mv': float(peak.max()), 'dt_ms': dt, 'duration_ms': float(times[-1]), 'amplitude_na': amplitude, 'delay_ms': delay, 'pulse_duration_ms': pulse_duration, 'passive_time_constant_ms': float(capacitance[root] / leak[root])}


def trace_svg(times, traces, title):
    low, high = min(-66., float(traces[:, :2].min()) - 1), max(-60., float(traces[:, :2].max()) + 1)
    def path(column):
        return ' '.join(f'{"M" if i == 0 else "L"}{65 + t / times[-1] * 680:.2f},{270 - (v - low) / (high - low) * 225:.2f}' for i, (t, v) in enumerate(zip(times, traces[:, column])))
    return f'''<svg xmlns="http://www.w3.org/2000/svg" width="800" height="350" viewBox="0 0 800 350" role="img" aria-label="Passive cable voltage response"><rect width="800" height="350" fill="white"/><g font-family="sans-serif" font-size="13" fill="#17212b"><text x="65" y="24">{title}: passive cable demonstration</text><path d="M65 45 V270 H745" fill="none" stroke="#aaa"/><text x="65" y="294">0</text><text x="720" y="294">{times[-1]:g} ms</text><text x="65" y="325">Soma (blue) · farthest node (orange) · no active channels or synapses</text><text x="10" y="49">{high:.1f}</text><text x="10" y="270">{low:.1f}</text><text x="10" y="150">mV</text></g><path d="{path(0)}" fill="none" stroke="#1764c0" stroke-width="2"/><path d="{path(1)}" fill="none" stroke="#bd5900" stroke-width="2"/></svg>'''


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', type=Path, default=DEFAULT_DATA)
    parser.add_argument('--output', type=Path, default=Path('outputs'))
    parser.add_argument('--cells', nargs='+', default=SAMPLE)
    parser.add_argument('--dt', type=float, default=.025)
    parser.add_argument('--duration', type=float, default=25.)
    parser.add_argument('--amplitude', type=float, default=1.)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    summaries = []
    for cell_id in args.cells:
        matches = list(args.data.glob(f'SWC_*/{cell_id}.swc'))
        if len(matches) != 1:
            parser.error(f'Expected one morphology for {cell_id}, found {len(matches)}')
        times, traces, summary = simulate(matches[0], duration=args.duration, dt=args.dt, amplitude=args.amplitude)
        with (args.output / f'passive_{cell_id}.csv').open('w', newline='') as stream:
            writer = csv.writer(stream)
            writer.writerow(['time_ms', 'soma_mv', 'farthest_node_mv', 'minimum_node_mv', 'maximum_node_mv'])
            writer.writerows(zip(times, *traces.T))
        (args.output / f'passive_{cell_id}.svg').write_text(trace_svg(times, traces, cell_id))
        summaries.append(summary)
        print(json.dumps(summary), flush=True)
    (args.output / 'passive_summary.json').write_text(json.dumps({'model': 'uncalibrated passive cable; no active channels, synapses, or extracellular stimulation', 'cells': summaries}, indent=2) + '\n')
