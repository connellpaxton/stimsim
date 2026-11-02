"""Check the downloaded xtra mechanism against the supplied Python interface."""

import importlib.util
import os
from pathlib import Path
import unittest

os.environ['NEURON_MODULE_OPTIONS'] = os.environ.get('NEURON_MODULE_OPTIONS', '') + ' -nogui'
MECHANISMS = Path(__file__).resolve().parents[1] / 'mechanisms'


@unittest.skipUnless(importlib.util.find_spec('neuron') and (MECHANISMS / 'x86_64' / 'libnrnmech.so').exists(), 'Requires NEURON and compiled xtra')
class XtraTests(unittest.TestCase):
    def test_coordinate_interface_and_external_voltage_pointer(self):
        import neuron
        from neuron import h
        if not hasattr(h, 'stim_xtra'):
            self.assertTrue(neuron.load_mechanisms(str(MECHANISMS)))
        section = h.Section(name='xtra_test')
        section.nseg = 3
        section.insert('pas')
        section.insert('extracellular')
        section.insert('xtra')
        for segment in section:
            h.setpointer(segment._ref_e_extracellular, 'ex', segment.xtra)
            segment.xtra.x, segment.xtra.y, segment.xtra.z = 10., 20., 30.
            segment.xtra.es = 2.
        h.stim_xtra = .5
        h.load_file('stdrun.hoc')
        h.finitialize(-65)
        for segment in section:
            self.assertEqual((segment.x_xtra, segment.y_xtra, segment.z_xtra), (10., 20., 30.))
            self.assertAlmostEqual(segment.e_extracellular, 1.)
        h.stim_xtra = 0.
        h.fadvance()
        for segment in section:
            self.assertAlmostEqual(segment.e_extracellular, 0.)


if __name__ == '__main__':
    unittest.main()
