"""Compatibility launcher for `python -m neuronsim`; see docs/CODE_WALKTHROUGH.md."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from neuronsim.cli import main
from neuronsim.network import synaptic_reversal

__all__ = ["main", "synaptic_reversal"]

if __name__ == "__main__":
    main()
