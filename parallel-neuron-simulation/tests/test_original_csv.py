"""Verify CSV parsing preserves the supplied model's synaptic connections."""

import csv
import importlib.util
import os
from pathlib import Path
import sys
import tempfile
import unittest

os.environ['NEURON_MODULE_OPTIONS'] = os.environ.get('NEURON_MODULE_OPTIONS', '') + ' -nogui'
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'model'))


@unittest.skipUnless(importlib.util.find_spec('neuron'), 'Requires the NEURON environment')
class OriginalCsvTests(unittest.TestCase):
    def test_precise_ids_coordinates_and_synapse_parameters(self):
        from neuron import h
        from cell import SynapseConstructor

        class TestCell:
            def __init__(self, name):
                self.name = name
                self.section = h.Section(name='csv_test', cell=self)

            def __str__(self):
                return self.name

        ids = ['864691135640672571', '864691135993429185', '864691136025333561']
        cells = [TestCell(name) for name in ids]
        fields = ['pre_pt_root_id', 'post_pt_root_id'] + [f'{side}_pt_position_{axis}' for side in ('pre', 'post') for axis in 'xyz']
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'contacts.csv'
            with path.open('w', newline='') as stream:
                writer = csv.writer(stream)
                writer.writerow(fields)
                for i in range(3):
                    writer.writerow([ids[i], ids[(i + 1) % 3], 1000., 2000., 3000., 4000., 5000., 6000.])
                writer.writerow(['864691135640672572', ids[0], 1, 2, 3, 4, 5, 6])
            constructor = SynapseConstructor()
            constructor.get_all_cells = lambda: (cells, ids)
            looked_up = []

            def lookup(cell, xyz, units):
                looked_up.append((str(cell), xyz, units))
                return cell.section(.5)

            constructor.find_synapse_segment = lookup
            constructor.insert_synapses(path)
            self.assertEqual(len(constructor.post_syns), 3)
            self.assertEqual(len(constructor.pre_syns), 3)
            self.assertEqual(looked_up[0], (ids[1], [4000., 5000., 6000.], 'nm'))
            self.assertEqual(looked_up[1], (ids[0], [1000., 2000., 3000.], 'nm'))
            for i, connection in enumerate(constructor.pre_syns):
                self.assertEqual(connection.preseg().sec, cells[i].section)
                self.assertEqual(connection.postseg().sec, cells[(i + 1) % 3].section)
                self.assertAlmostEqual(connection.weight[0], .01)
                self.assertAlmostEqual(connection.delay, 1.)
                self.assertAlmostEqual(connection.threshold, 0.)
            for synapse in constructor.post_syns:
                self.assertAlmostEqual(synapse.tau1, .5)
                self.assertAlmostEqual(synapse.tau2, 3.)
                self.assertAlmostEqual(synapse.e, 0.)


if __name__ == '__main__':
    unittest.main()
