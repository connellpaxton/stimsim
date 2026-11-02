"""Compare identical synaptic inputs at two sites on one recorded Minnie cell."""

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import runpy
import sys

from minnie_network import PROFILES
from reconstructed_cell import ElectricalProfile, ReconstructedCell

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', type=Path, required=True)
    parser.add_argument('--morphologies', type=Path, required=True)
    parser.add_argument('--root-id', required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--d-lambda', type=float, default=.1)
    parser.add_argument('--dt', type=float, default=.025)
    args = parser.parse_args()
    if args.d_lambda <= 0 or args.dt <= 0:
        parser.error('Discretization parameters must be positive')
    data, morph, output = (p.resolve() for p in (args.data, args.morphologies, args.output))
    output.mkdir(parents=True, exist_ok=True)
    os.environ['NEURON_MODULE_OPTIONS'] = '-nogui'
    os.environ['MPLBACKEND'] = 'Agg'
    os.environ['OPENBLAS_NUM_THREADS'] = '1'
    os.environ['OMP_NUM_THREADS'] = '1'
    import numpy as np
    import neuron
    from neuron import h
    from netpyne import sim

    source = ROOT / 'vendor/m1_260015'
    manifest = json.loads((source / 'SOURCE_MANIFEST.json').read_text())
    for item in manifest['files']:
        if hashlib.sha256((source / item['path']).read_bytes()).hexdigest() != item['sha256']:
            raise RuntimeError('Published source changed')
    profiles = json.loads((morph / 'electrical_profiles.json').read_text())
    if profiles['source_commit'] != manifest['commit']:
        raise RuntimeError('Electrical prior source mismatch')
    with (data / 'nodes.csv').open() as stream:
        node = next(n for n in csv.DictReader(stream) if n['root_id'] == args.root_id)
    record = next(r for r in json.loads((morph / 'manifest.json').read_text())['records'] if r['root_id'] == args.root_id)
    swc = morph / 'swc' / (args.root_id + '.swc')
    if record['nucleus_id'] != node['nucleus_id'] or hashlib.sha256(swc.read_bytes()).hexdigest() != record['import_sha256']:
        raise RuntimeError('Reconstruction identity or checksum mismatch')
    source_dir = source / 'sim'
    neuron.load_mechanisms(str(source_dir))
    os.chdir(source_dir)
    sys.path.insert(0, str(source_dir))
    global cfg
    cfg = runpy.run_path(str(source_dir / 'cfg.py'))['cfg']
    cfg.addConn = cfg.addLongConn = cfg.addSubConn = cfg.addIClamp = 0
    cfg.recordLFP = []; cfg.recordDipole = False; cfg.recordDipoles = False
    cfg.analysis = {}; cfg.verbose = False; cfg.recordTraces = {}; cfg.recordCells = []
    parameters = runpy.run_path(str(source_dir / 'netParams.py'))['netParams']
    sim.initialize(simConfig=cfg, netParams=parameters)
    sim.pc.nthread(1)
    h.load_file('stdrun.hoc'); h.load_file('import3d.hoc')
    label, kind, model = PROFILES[node['manual_cell_type']]
    cell = ReconstructedCell(swc, {**node, 'profile': label, 'cell_type': kind, 'cell_model': model},
                             0, ElectricalProfile(profiles['profiles'][label]), sim.pc, d_lambda=args.d_lambda)
    soma = cell.secs['soma']['hObj']
    h.distance(0, .5, sec=soma)
    branches = [(name, info['hObj']) for name, info in cell.secs.items() if cell.section_regions[name] == 4]
    if not branches:
        raise RuntimeError('Probe requires an apical reconstruction')
    sites = []
    for target in (60., 250.):
        candidates = []
        for name, section in branches:
            start = float(h.distance(0., sec=section))
            x = float(np.clip((target - start) / section.L, .05, .95))
            candidates.append((start + x * section.L, name, x))
        sites.append(min(candidates, key=lambda s: abs(s[0] - target)))
    if sites[0] == sites[1]:
        raise RuntimeError('Distinct proximal and distal sites are required')
    synapses, sources, owners = [], [], []
    for distance, name, x in sites:
        syn = h.MyExp2SynBB(cell.secs[name]['hObj'](x))
        for key, value in parameters.synMechParams['AMPA'].items():
            if key != 'mod': setattr(syn, key, value)
        synapses.append(syn); sources.append(h.NetCon(None, syn)); owners.append(syn)
    anchor = h.IClamp(soma(.5)); anchor.amp = anchor.dur = 0.
    vectors = [h.Vector().record(anchor, soma(.5)._ref_v)]
    for owner, (_, name, x) in zip(owners, sites):
        vectors.append(h.Vector().record(owner, cell.secs[name]['hObj'](x)._ref_v))
    h.celsius = 34.; h.dt = args.dt; h.steps_per_ms = 1. / args.dt; h.CVode().active(0)
    initializer = h.FInitializeHandler(0, cell.initV)
    sim.pc.set_maxstep(10.)
    recordings = {}
    for condition in ('baseline', 'proximal', 'distal'):
        for source in sources: source.weight[0] = .001
        h.finitialize(float(cfg.hParams['v_init']))
        if condition != 'baseline': sources[0 if condition == 'proximal' else 1].event(50.)
        sim.pc.psolve(100.)
        recordings[condition] = np.array([np.asarray(v).copy() for v in vectors])
    if not all(np.isfinite(v).all() for v in recordings.values()):
        raise RuntimeError('Nonfinite branch-probe voltage')
    time_ms = np.arange(recordings['baseline'].shape[1]) * args.dt
    result = {'root_id': args.root_id, 'manual_type': node['manual_cell_type'], 'M1_profile': label,
              'swc_sha256': record['import_sha256'], 'conductance_us': .001, 'event_ms': 50.,
              'dt_ms': args.dt, 'd_lambda': args.d_lambda,
              'segments': sum(info['hObj'].nseg for info in cell.secs.values()),
              'interpretation': 'Two modeled AMPA sites on the same recorded active tree. Responses include geometry and transferred channel-density priors; this does not isolate geometry experimentally.',
              'sites': []}
    for i, (condition, site) in enumerate(zip(('proximal', 'distal'), sites), 1):
        delta = recordings[condition] - recordings['baseline']
        window = time_ms >= 50.
        soma_peak = float(delta[0, window].max()); local_peak = float(delta[i, window].max())
        result['sites'].append({'condition': condition, 'section': site[1], 'x': float(site[2]),
                                'soma_path_um': site[0], 'soma_peak_delta_mv': soma_peak,
                                'local_peak_delta_mv': local_peak,
                                'soma_to_local_peak_ratio': soma_peak / local_peak})
    np.savez_compressed(output / 'recordings.npz', time_ms=time_ms, **recordings)
    (output / 'summary.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2), flush=True)


if __name__ == '__main__':
    main()
