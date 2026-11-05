"""Validate synapse accounting and compiler-based thread restrictions."""

from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'tools'))
from run_mouse_m1 import connection_totals, non_thread_safe_mechanisms, scale_populations


class MouseM1Tests(unittest.TestCase):
    def test_density_and_fixed_counts_are_both_scaled(self):
        populations = {'cortex': {'density': 1000., 'ynormRange': [.1, .3]}, 'inputs': {'numCells': 1000}}
        scale_populations(populations, .1)
        self.assertEqual(populations['cortex']['density'], 100.)
        self.assertEqual(populations['inputs']['numCells'], 100)
        self.assertEqual(populations['cortex']['ynormRange'], [.1, .3])

    def test_receptors_are_separate_from_connection_sites(self):
        cell = SimpleNamespace(gid=2, secs={'soma': {'synMechs': [
            {'label': 'AMPA', 'hObj': object()},
            {'label': 'NMDA', 'hObj': object()},
            {'label': 'GABAA', 'hObj': None},
        ]}}, conns=[
            {'preGid': 1, 'sec': 'soma', 'loc': .5, 'hObj': object()},
            {'preGid': 1, 'sec': 'soma', 'loc': .5, 'hObj': object()},
        ])
        result = connection_totals([cell])
        self.assertEqual(result['synaptic_point_processes'], 2)
        self.assertEqual(result['synaptic_netcons'], 2)
        self.assertEqual(result['unique_connection_site_keys'], 1)
        self.assertEqual(result['directed_cell_pairs'], 1)

    def test_compiler_flags_identify_actual_mechanism_suffix(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'x86_64'
            path.mkdir()
            (path / 'legacy_filename.cpp').write_text('#define NRN_VECTORIZED 0\nstatic const char *_mechanism[] = {"7.7.0", "catcb", 0};')
            (path / 'safe.cpp').write_text('#define NRN_VECTORIZED 1\nstatic const char *_mechanism[] = {"7.7.0", "safe", 0};')
            self.assertEqual(non_thread_safe_mechanisms(Path(directory)), {'catcb'})


if __name__ == '__main__':
    unittest.main()
