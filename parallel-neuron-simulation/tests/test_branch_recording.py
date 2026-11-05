"""Check recorded cable polylines and voltage timing used by the video."""

import importlib.util
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


@unittest.skipUnless(importlib.util.find_spec('neuron'), 'Requires NEURON')
class BranchRecordingTests(unittest.TestCase):
    def test_curved_cable_is_preserved_and_voltage_samples_align(self):
        code = '''
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from neuron import h
from branch_recording import BranchRecording
h.load_file('stdrun.hoc'); h.dt=.025
section=h.Section(); h.pt3dclear(sec=section)
for x,y in [(0,0),(100,0),(100,100)]:h.pt3dadd(x,y,0,2,sec=section)
section.nseg=3; section.insert('pas'); section.e_pas=-65.; section.g_pas=.0001
cell=SimpleNamespace(soma=[section])
record=BranchRecording({'cell':cell},['cell'],{'cell':{'dend_0':section}},Path('geometry.json'))
stim=h.IClamp(section(.5)); stim.amp=.1; stim.delay=1.; stim.dur=1.
soma=h.Vector().record(stim,section(.5)._ref_v)
h.finitialize(-65.); h.continuerun(5.)
times=np.arange(len(np.asarray(soma)[::4]))*.1
record.save(Path('.'),times)
data=np.load('branch_recordings.npz')
assert np.array_equal(data['segment_mv'][1],np.asarray(soma)[::4].astype(np.float32))
points=[np.array(row['points_um']) for row in record.geometry['rows']]
assert any(np.any(np.all(p==[100,0,0],axis=1)) for p in points)
assert np.allclose(sum(np.linalg.norm(np.diff(p,axis=0),axis=1).sum() for p in points),200.)
assert np.allclose(points[0][0],[0,0,0]) and np.allclose(points[-1][-1],[100,100,0])
assert data['segment_mv'].max()>-65.
'''
        with tempfile.TemporaryDirectory() as directory:
            env = dict(os.environ, NEURON_MODULE_OPTIONS='-nogui', PYTHONPATH=str(ROOT / 'tools'))
            result = subprocess.run([sys.executable, '-c', code], cwd=directory, env=env, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)


if __name__ == '__main__': unittest.main()
