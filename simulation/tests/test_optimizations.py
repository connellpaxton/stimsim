"""Check that zero-channel removal preserves every nonzero section."""

import importlib.util
import os
from pathlib import Path
import sys
import unittest

os.environ['NEURON_MODULE_OPTIONS'] = os.environ.get('NEURON_MODULE_OPTIONS', '') + ' -nogui'
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'tools'))
from simulation_optimizations import prune_zero_conductance_channels


@unittest.skipUnless(importlib.util.find_spec('neuron') and (ROOT / 'mechanisms/x86_64/libnrnmech.so').exists(), 'Requires compiled NEURON mechanisms')
class OptimizationTests(unittest.TestCase):
    def test_only_entirely_zero_mechanisms_in_selected_sections_are_removed(self):
        import neuron
        from neuron import h
        if not hasattr(h, 'stim_xtra'):
            neuron.load_mechanisms(str(ROOT / 'mechanisms'))
        section = h.Section(name='pruning_test')
        other = h.Section(name='unselected_test')
        section.nseg = 3
        for sec in [section, other]:
            sec.insert('Nap_Et2')
            sec.insert('K_Tst')
            sec.gNap_Et2bar_Nap_Et2 = 0
            sec.gK_Tstbar_K_Tst = 0
        section(.5).gK_Tstbar_K_Tst = .01
        removed = prune_zero_conductance_channels([section])
        self.assertEqual(removed, {'Nap_Et2': 3})
        self.assertFalse(h.ismembrane('Nap_Et2', sec=section))
        self.assertTrue(h.ismembrane('K_Tst', sec=section))
        self.assertAlmostEqual(section(.5).gK_Tstbar_K_Tst, .01)
        self.assertTrue(h.ismembrane('Nap_Et2', sec=other))
        self.assertEqual(prune_zero_conductance_channels([section]), {})
        section.gK_Tstbar_K_Tst = 0
        self.assertEqual(prune_zero_conductance_channels([section]), {'K_Tst': 3})


if __name__ == '__main__':
    unittest.main()
