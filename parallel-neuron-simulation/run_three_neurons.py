"""Original archive: exactly three L5 wrapper neurons and three coordinate-mapped contacts.

Run with --describe for the exact dataset identity, source URLs and assumptions.
See docs/DATASET_ENTRY_POINTS.md for preparation and execution commands.
"""

from neuronsim.datasets import THREE_NEURONS, launch


if __name__ == "__main__":
    # The shared engine owns construction, synapses, recording and integration.
    launch(THREE_NEURONS)
