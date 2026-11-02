"""Compatibility import; the cell implementation is in neuronsim/cell.py."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from neuronsim.cell import PublishedSWCCell

__all__ = ['PublishedSWCCell']
