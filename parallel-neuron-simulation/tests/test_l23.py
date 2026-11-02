"""Check layer selection and the published-profile SWC adapter."""

import importlib.util
import os
from pathlib import Path
import sys
import tempfile
import unittest

os.environ['NEURON_MODULE_OPTIONS'] = os.environ.get('NEURON_MODULE_OPTIONS', '') + ' -nogui'
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'tools'))
from run_l23 import choose_cells


class SelectionTests(unittest.TestCase):
    def test_contacts_drive_selection_without_changing_ids(self):
        ids = {'864691134941559651', '864691134941559652', '864691134941559653', '864691134941559654'}
        a, b, c, d = sorted(ids)
        contacts = [{'pre_pt_root_id': pre, 'post_pt_root_id': post} for pre, post in [(a, b), (a, b), (b, c), (c, d)]]
        self.assertEqual(choose_cells(ids, contacts, 3), [a, b, c])
        self.assertEqual(choose_cells(ids, contacts, 0), [a, b, c, d])


@unittest.skipUnless(importlib.util.find_spec('neuron') and (ROOT / 'vendor/aberra2018/mechanisms/x86_64/libnrnmech.so').exists(), 'Requires published compiled mechanisms')
class SWCProfileTests(unittest.TestCase):
    def test_region_parameters_coordinates_and_multiple_cells(self):
        import neuron
        from neuron import h
        if not hasattr(h, 'stim_xtra'):
            neuron.load_mechanisms(str(ROOT / 'vendor/aberra2018/mechanisms'))
        h.load_file('stdrun.hoc')
        h.load_file('import3d.hoc')
        h.load_file(str(ROOT / 'vendor/aberra2018/cells/L23_PC_cADpyr229_1/biophysics.hoc'))
        from published_swc import PublishedSWCCell
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / '864691134941559651.swc'
            path.write_text('0 1 0 0 0 5 -1\n1 2 10 0 0 .5 0\n2 2 20 0 0 .5 1\n3 3 0 10 0 1 0\n4 3 0 20 0 1 3\n5 4 0 0 10 1 0\n6 4 0 0 20 1 5\n')
            profile = h.cADpyr229_biophys()
            first = PublishedSWCCell(path, profile)
            second = PublishedSWCCell(path, profile)
            self.assertNotEqual(str(first), str(second))
            self.assertGreater(first.soma[0].L, 0)
            self.assertAlmostEqual(first.soma[0](.5).gNaTs2_tbar_NaTs2_t, .926705)
            self.assertAlmostEqual(first.axon[0](.5).gNaTa_tbar_NaTa_t, 3.429725)
            self.assertAlmostEqual(first.dend[0](.5).gIhbar_Ih, .00008)
            self.assertFalse(h.ismembrane('NaTs2_t', sec=first.dend[0]))
            self.assertAlmostEqual(first.apic[0](.5).gNaTs2_tbar_NaTs2_t, .012009)
            segment, distance = first.nearest_segment(first.segment_xyz_um[-1] * 1000)
            self.assertAlmostEqual(distance, 0)
            self.assertEqual(segment.sec, first.segments[-1].sec)
            self.assertAlmostEqual(segment.x, first.segments[-1].x)


if __name__ == '__main__':
    unittest.main()
