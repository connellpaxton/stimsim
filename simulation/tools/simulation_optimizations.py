"""Compatibility import for the validated static zero-channel optimization."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from neuronsim.optimizations import prune_zero_conductance_channels

__all__ = ['prune_zero_conductance_channels']
