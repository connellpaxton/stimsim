"""Simulation assumptions, with units; intrinsic channel densities live in HOC."""

from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "vendor/aberra2018"
PROFILE = SOURCE / "cells/L4_LBC_cACint209_1/biophysics.hoc"


# SIMULATION CORE: numerical/physical constants shared by both dataset engines.
@dataclass(frozen=True)
class ModelParameters:
    """Defaults used by the validated runs, rather than Minnie measurements."""

    dt_ms: float = 0.025
    temperature_c: float = 34.0
    initial_voltage_mv: float = -65.0
    contact_conductance_us: float = 0.01
    synapse_rise_ms: float = 0.5
    synapse_decay_ms: float = 3.0
    delay_ms: float = 1.0
    placement_seed: int = 20261027
    record_stride: int = 4
    progress_interval_ms: float = 10.0


PARAMETERS = ModelParameters()


@dataclass(frozen=True)
class Protocol:
    """One current-injection condition, including a transmission control."""

    label: str
    cells: list[str]
    amplitude_na: float
    onset_ms: float
    pulse_ms: float
    duration_ms: float
    weight_factor: float = 1.0

    def stimulus_metadata(self):
        return {
            "cells": self.cells,
            "amplitude_na": self.amplitude_na,
            "onset_ms": self.onset_ms,
            "duration_ms": self.pulse_ms,
        }


def make_protocols(args, stimulated):
    """Preserve the existing CLI's condition order and stimulus settings."""
    factor = 0.0 if args.control_only else 1.0
    protocols = []
    if args.duration_ms is not None:
        label = "long_protocol_control" if args.control_only else "long_protocol"
        protocols.append(Protocol(label, stimulated, 0.5, 10.0, 50.0, args.duration_ms, factor))
        if args.paired_control:
            protocols.append(Protocol("long_protocol_control", stimulated, 0.5, 10.0, 50.0, args.duration_ms, 0.0))
        window = getattr(args, "validation_window_ms", None)
        if window is not None:
            protocols.append(Protocol("validation_zero_transmission", stimulated, 0.5, 10.0, 50.0, window, 0.0))
            protocols.append(Protocol("validation_replay", stimulated, 0.5, 10.0, 50.0, window, factor))
        return protocols

    label = "original_stimulus_control" if args.control_only else "original_stimulus"
    protocols.append(Protocol(label, stimulated[:1], 1.0, 0.0, 1.0, 25.0, factor))
    if args.paired_control:
        protocols.append(Protocol("original_stimulus_control", stimulated[:1], 1.0, 0.0, 1.0, 25.0, 0.0))
    if args.matched_protocol:
        label = "matched_protocol_control" if args.control_only else "matched_protocol"
        protocols.append(Protocol(label, stimulated, 0.5, 10.0, 50.0, 100.0, factor))
        if args.paired_control:
            protocols.append(Protocol("matched_protocol_control", stimulated, 0.5, 10.0, 50.0, 100.0, 0.0))
    return protocols
