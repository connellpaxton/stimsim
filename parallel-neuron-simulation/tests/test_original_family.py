"""Check the original-family regional profile on a recorded-style branched tree."""

import importlib.util
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'tools'))
from run_minnie_241165 import synaptic_reversal


class PolaritySelectionTests(unittest.TestCase):
    def test_unknown_transmitter_class_is_rejected(self):
        with self.assertRaises(ValueError): synaptic_reversal({'ei_class': 'unknown'}, inhibition=True)


@unittest.skipUnless(importlib.util.find_spec('neuron'), 'Requires NEURON')
class SynapticPolarityTests(unittest.TestCase):
    def test_inhibitory_event_hyperpolarizes_with_positive_conductance(self):
        code = '''
from neuron import h
from run_minnie_241165 import synaptic_reversal
h.load_file('stdrun.hoc'); h.dt=.025
sections=[]; synapses=[]; sources=[]; vectors=[]
for kind in ['E','I']:
    sec=h.Section(); sec.L=sec.diam=20.; sec.insert('pas'); sec.g_pas=.0001; sec.e_pas=-65.
    syn=h.Exp2Syn(sec(.5)); syn.tau1=.5; syn.tau2=3.
    syn.e=synaptic_reversal({'ei_class':kind}, inhibition=True)
    nc=h.NetCon(None,syn); nc.weight[0]=.01
    sections.append(sec); synapses.append(syn); sources.append(nc)
    vectors.append(h.Vector().record(syn,sec(.5)._ref_v))
h.finitialize(-65.)
for nc in sources:nc.event(1.)
h.continuerun(5.)
assert max(vectors[0])>-64.
assert min(vectors[1])<-66.
assert all(nc.weight[0]>0 for nc in sources)
assert synaptic_reversal({'ei_class':'I'}, inhibition=False)==0.
'''
        env = dict(os.environ, NEURON_MODULE_OPTIONS='-nogui', PYTHONPATH=str(ROOT / 'tools'))
        result = subprocess.run([sys.executable, '-c', code], env=env, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)


@unittest.skipUnless(importlib.util.find_spec('neuron') and (ROOT / 'vendor/aberra2018/mechanisms/x86_64/libnrnmech.so').exists(), 'Requires published compiled mechanisms')
class OriginalFamilyTests(unittest.TestCase):
    def test_published_regional_values_gradient_and_original_discretization(self):
        code = '''
import math, os
from pathlib import Path
import neuron
from neuron import h
from published_swc import PublishedSWCCell
root = Path(os.environ['NEURONSIM_TEST_ROOT'])
neuron.load_mechanisms(str(root/'vendor/aberra2018/mechanisms'))
h.load_file('stdrun.hoc'); h.load_file('import3d.hoc')
h.load_file(str(root/'vendor/aberra2018/cells/L4_LBC_cACint209_1/biophysics.hoc'))
cell = PublishedSWCCell(Path(os.environ['NEURONSIM_TEST_SWC']), h.cACint209_biophys(), original_discretization=True, coordinate_cache=False)
assert cell.segment_xyz_um is None
assert abs(cell.soma[0](.5).gNaTs2_tbar_NaTs2_t-.395881)<1e-12
assert abs(cell.axon[0](.5).gNaTa_tbar_NaTa_t-3.993125)<1e-12
assert abs(cell.soma[0](.5).gamma_CaDynamics_E2-.000814)<1e-12
values = []
for section in cell.sections:
    full = int((section.L/(.01*h.lambda_f(100,sec=section))+.999)/2)*2+1
    assert section.nseg == int(full/9/2)*2+1
for section in cell.apic:
    for seg in section:
        # Upstream for(x) also visits x=1, overwriting the last segment's
        # density with the endpoint value. Preserve this source behavior.
        x = 1. if seg.x > 1.-1./section.nseg else seg.x
        expected = (-.869600+2.087000*math.exp(h.distance(x,sec=section)*.003000))*.000023
        # The upstream distribute() formats physical distance with C's %g.
        assert math.isclose(seg.gIhbar_Ih, expected, rel_tol=2e-6, abs_tol=1e-12), (seg.x, h.distance(seg.x,sec=section), seg.gIhbar_Ih, expected)
        values.append(seg.gIhbar_Ih)
assert max(values)>min(values)
assert not h.ismembrane('Nap_Et2',sec=cell.axon[0])
pc=h.ParallelContext(); pc.set_maxstep(1); h.finitialize(-65); pc.psolve(.1)
assert all(math.isfinite(seg.v) for seg in cell.segments)
'''
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / '864691134941559651.swc'
            path.write_text('0 1 0 0 0 5 -1\n1 2 0 -5 0 .5 0\n2 2 0 -100 0 .3 1\n3 3 -5 0 0 1 0\n4 3 -150 0 0 .5 3\n5 4 0 5 0 1 0\n6 4 0 500 0 .4 5\n')
            env = dict(os.environ, NEURON_MODULE_OPTIONS='-nogui', PYTHONPATH=str(ROOT / 'tools'),
                       NEURONSIM_TEST_ROOT=str(ROOT), NEURONSIM_TEST_SWC=str(path))
            result = subprocess.run([sys.executable, '-c', code], env=env, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)


if __name__ == '__main__': unittest.main()
