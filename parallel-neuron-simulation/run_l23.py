"""Original archive: SWC_L23/ plus L23_all.csv; all 313 valid neurons by default.

Run with --describe for the exact dataset identity, source URLs and assumptions.
See docs/DATASET_ENTRY_POINTS.md for preparation and execution commands.
"""

from neuronsim.datasets import L23, launch


if __name__ == "__main__":
    # The shared engine owns construction, synapses, recording and integration.
    launch(L23)
