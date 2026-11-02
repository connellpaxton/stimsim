"""Preflight and run the supplied active NEURON model without editing the download."""

import argparse
import ast
import csv
import importlib.util
import os
from pathlib import Path
import runpy
import sys

from analyze import DEFAULT_DATA


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', type=Path, default=DEFAULT_DATA)
    parser.add_argument('--model', type=Path, default=Path(__file__).resolve().parents[1] / 'model', help='Working copy of the supplied scripts')
    parser.add_argument('--mechanisms', type=Path, help='Directory containing compiled custom mechanisms')
    parser.add_argument('--output', type=Path, default=Path('outputs'))
    parser.add_argument('--check', action='store_true', help='Only check prerequisites')
    args = parser.parse_args()
    missing = [name for name in ['numpy', 'matplotlib', 'neuron'] if importlib.util.find_spec(name) is None]
    if missing:
        print('Missing Python packages: ' + ', '.join(missing))
    tree = ast.parse((args.model / 'cell.py').read_text())
    required = {n.args[0].value for n in ast.walk(tree) if isinstance(n, ast.Call) and isinstance(n.func, ast.Attribute) and n.func.attr == 'insert' and n.args and isinstance(n.args[0], ast.Constant)}
    if 'neuron' in missing:
        print(f'Custom mechanisms must also be supplied ({len(required - {"pas", "extracellular"})} referenced; no MOD files included in download).')
        return 1
    os.environ['NEURON_MODULE_OPTIONS'] = os.environ.get('NEURON_MODULE_OPTIONS', '') + ' -nogui'
    os.environ['MPLBACKEND'] = 'Agg'
    os.environ.setdefault('MPLCONFIGDIR', str(args.output.resolve() / 'matplotlib-cache'))
    import neuron
    from neuron import h
    if args.mechanisms:
        neuron.load_mechanisms(str(args.mechanisms.resolve()))
    available = set()
    mechanisms = h.MechanismType(0)
    name = h.ref('')
    for i in range(int(mechanisms.count())):
        mechanisms.select(i)
        mechanisms.selected(name)
        available.add(str(name[0]))
    absent = sorted(required - available)
    if absent:
        print('Missing NEURON mechanisms: ' + ', '.join(absent))
    if missing or absent:
        return 1
    if args.check:
        print('Prerequisites available. Original model was not run.')
        return 0
    sys.dont_write_bytecode = True
    sys.path.insert(0, str(args.model.resolve()))
    os.environ['NEURONSIM_DATA'] = str(args.data.resolve())
    os.environ.pop('NEURONSIM_ANIMATION', None)
    import matplotlib.pyplot as plt
    plt.show = lambda: None
    namespace = runpy.run_path(str(args.model / 'wrapper-synaptic-connections.py'), run_name='__main__')
    cells, synapses = namespace['cells'], namespace['sc']
    time, voltage = namespace['t'], namespace['v_vecs']
    args.output.mkdir(parents=True, exist_ok=True)
    with (args.output / 'original_voltage.csv').open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['time_ms'] + [f'{c}_axon_mv' for c in cells])
        writer.writerows(zip(time, *voltage))
    for index, number in enumerate(plt.get_fignums(), start=1):
        plt.figure(number).savefig(args.output / f'original_voltage_{index}.png', dpi=160, bbox_inches='tight')
    print(f'Ran original active model: {len(cells)} cells, {len(synapses.post_syns)} contacts, 25 ms. Saved original_voltage.csv.')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
