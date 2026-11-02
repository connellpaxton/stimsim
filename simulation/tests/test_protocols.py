"""Check that transmission controls preserve the experimental stimulus."""

from dataclasses import asdict
import unittest

from neuronsim.cli import parse_args
from neuronsim.config import make_protocols


class ProtocolTests(unittest.TestCase):
    def options(self, *extra):
        return parse_args([
            "--data", "data", "--morphologies", "morphologies", "--output", "outputs/test",
            *extra,
        ])

    def test_one_second_protocol_preserves_exact_neuron_ids(self):
        roots = ["864691135395143157", "864691136444943747"]
        protocols = make_protocols(self.options("--duration-ms", "1000"), roots)
        self.assertEqual(len(protocols), 1)
        self.assertEqual(protocols[0].cells, roots)
        self.assertEqual(protocols[0].duration_ms, 1000.0)
        self.assertEqual(protocols[0].stimulus_metadata(), {
            "cells": roots, "amplitude_na": 0.5, "onset_ms": 10.0, "duration_ms": 50.0,
        })

    def test_paired_controls_change_transmission_but_preserve_stimulus(self):
        roots = ["864691135395143157", "864691136444943747"]
        protocols = make_protocols(self.options("--paired-control", "--matched-protocol"), roots)
        self.assertEqual([p.label for p in protocols], [
            "original_stimulus", "original_stimulus_control",
            "matched_protocol", "matched_protocol_control",
        ])
        for experiment, control in zip(protocols[::2], protocols[1::2]):
            left, right = asdict(experiment), asdict(control)
            self.assertEqual(left.pop("weight_factor"), 1.0)
            self.assertEqual(right.pop("weight_factor"), 0.0)
            left.pop("label")
            right.pop("label")
            self.assertEqual(left, right)

    def test_validation_replay_restores_transmission_and_retains_stimulus(self):
        roots = ["864691135395143157", "864691136444943747"]
        protocols = make_protocols(self.options("--duration-ms", "1000", "--validation-window-ms", "100"), roots)
        self.assertEqual([p.label for p in protocols], ["long_protocol", "validation_zero_transmission", "validation_replay"])
        self.assertEqual([p.duration_ms for p in protocols], [1000.0, 100.0, 100.0])
        self.assertEqual([p.weight_factor for p in protocols], [1.0, 0.0, 1.0])
        self.assertTrue(all(p.stimulus_metadata() == protocols[0].stimulus_metadata() for p in protocols))


if __name__ == "__main__":
    unittest.main()
