"""Regression checks for threshold events exactly at a replay cutoff."""

from copy import deepcopy
from pathlib import Path
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
from analyze_nas_trial import validate_replay


class ReplayValidationTests(unittest.TestCase):
    def recordings(self):
        main = {"root_ids": np.array(["864691135999999999"]),
                "time_ms": np.array([0, .025, .05, .075]),
                "soma_mv": np.array([[-65., -50., 10., -30.]]),
                "spike_times_ms": np.array([.025 + 1e-10, .05 + 1e-10]),
                "spike_root_ids": np.array(["864691135999999999"] * 2)}
        replay = deepcopy(main)
        replay["time_ms"] = replay["time_ms"][:3]
        replay["soma_mv"] = replay["soma_mv"][:, :3]
        replay["spike_times_ms"] = replay["spike_times_ms"][:1]
        replay["spike_root_ids"] = replay["spike_root_ids"][:1]
        return main, replay

    def test_cutoff_event_is_reported_separately_without_changing_raw_data(self):
        main, replay = self.recordings()
        original = main["spike_times_ms"].copy()
        check = validate_replay(main, replay, .05)
        self.assertEqual(check["interior_spikes_compared"], 1)
        self.assertEqual(check["main_spikes_at_window_boundary"], 1)
        self.assertEqual(check["replay_spikes_at_window_boundary"], 0)
        self.assertFalse(check["all_saved_arrays_identical_to_main_prefix"])
        self.assertTrue(check["all_interior_spike_times_and_ids_identical"])
        np.testing.assert_array_equal(main["spike_times_ms"], original)

    def test_missing_interior_spike_is_rejected(self):
        main, replay = self.recordings()
        replay["spike_times_ms"] = replay["spike_times_ms"][:0]
        replay["spike_root_ids"] = replay["spike_root_ids"][:0]
        with self.assertRaisesRegex(RuntimeError, "before cutoff"):
            validate_replay(main, replay, .05)

    def test_changed_voltage_at_cutoff_is_rejected(self):
        main, replay = self.recordings()
        replay["soma_mv"][0, -1] += .01
        with self.assertRaisesRegex(RuntimeError, "soma_mv"):
            validate_replay(main, replay, .05)


if __name__ == "__main__":
    unittest.main()
