"""Full eligible Minnie65 v1300 cohort: 70,995 neurons; predominantly dendrite-only geometry.

Run with --describe for the exact dataset identity, source URLs and assumptions.
See docs/DATASET_ENTRY_POINTS.md for preparation and execution commands.
"""

from neuronsim.datasets import WHOLE_MINNIE, launch


if __name__ == "__main__":
    # The shared engine owns construction, synapses, recording and integration.
    launch(WHOLE_MINNIE)
