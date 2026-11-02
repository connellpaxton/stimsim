"""Minnie column: 1,352 reconstructed neurons; v661 geometry linked to v943 wiring.

Run with --describe for the exact dataset identity, source URLs and assumptions.
See docs/DATASET_ENTRY_POINTS.md for preparation and execution commands.
"""

from neuronsim.datasets import COLUMN, launch


if __name__ == "__main__":
    # The shared engine owns construction, synapses, recording and integration.
    launch(COLUMN)
