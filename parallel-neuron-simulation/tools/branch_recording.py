"""Compatibility import; recorders live in neuronsim/recording.py."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from neuronsim.recording import BranchRecording

__all__ = ['BranchRecording']
