"""Command-line options for the Minnie / ModelDB 241165 runner."""

import argparse
import math
import os
from pathlib import Path


def parse_args(argv=None, *, description=__doc__, defaults=None):
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--data", type=Path, required=True, help="Prepared Minnie nodes.csv and edges.csv directory")
    parser.add_argument("--morphologies", type=Path, required=True, help="Validated manifest.json and swc/ directory")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=32, help="Native NEURON workers within one process")
    parser.add_argument("--paired-control", action="store_true", help="Repeat the condition with transmission weights zero")
    parser.add_argument("--control-only", action="store_true", help="Run zero-transmission protocols in a fresh process")
    parser.add_argument("--matched-protocol", action="store_true", help="Also run the earlier 100 ms, 32-cell current protocol")
    parser.add_argument("--stimulus-reference", type=Path, help="Reuse the injected IDs from a previous run summary")
    polarity = parser.add_mutually_exclusive_group()
    polarity.add_argument("--inhibition", action="store_true", help="Make outputs from Minnie inhibitory cells inhibitory")
    polarity.add_argument("--no-inhibition", dest="inhibition", action="store_false", help="Explicitly request all-excitatory outputs")
    parser.set_defaults(inhibition=False)
    parser.add_argument("--inhibitory-reversal", type=float, default=-80.0, help="Assumed inhibitory Exp2Syn reversal potential, mV")
    parser.add_argument("--duration-ms", type=float, help="Model time in milliseconds, with the near-start stimulus; 1000 means one model second")
    parser.add_argument("--validation-window-ms", type=float, help="After the long protocol, run a zero-transmission control and coupled replay for this many milliseconds")
    parser.add_argument("--branch-roots", type=Path, help="JSON list of neuron IDs to record at cable-segment resolution")
    parser.add_argument("--mpi-shards", type=Path, help="Destination-owned input shards prepared by python -m neuronsim.sharding")
    parser.set_defaults(**(defaults or {}))
    args = parser.parse_args(argv)
    if args.threads < 1:
        parser.error("Thread count must be positive")
    if args.control_only and args.paired_control:
        parser.error("Choose paired controls or a control-only run")
    if not -200.0 < args.inhibitory_reversal < 0.0:
        parser.error("Inhibitory reversal must be finite and between -200 and 0 mV")
    if args.duration_ms is not None and (not math.isfinite(args.duration_ms) or args.duration_ms <= 0):
        parser.error("Duration must be finite and positive")
    if args.validation_window_ms is not None:
        if (not math.isfinite(args.validation_window_ms) or args.validation_window_ms <= 0
                or args.duration_ms is None or args.validation_window_ms > args.duration_ms
                or args.control_only):
            parser.error("Validation window requires a coupled long protocol and must be positive and no longer than its duration")
    args.data = args.data.resolve()
    args.morphologies = args.morphologies.resolve()
    args.output = args.output.resolve()
    if args.mpi_shards:
        args.mpi_shards = args.mpi_shards.resolve()
    return args


def run_args(args):
    """Shared dispatch for the generic CLI and dataset-specific entry points."""
    args.output.mkdir(parents=True, exist_ok=True)
    # Set before NumPy/NEURON imports so BLAS does not create nested workers.
    os.environ["NEURON_MODULE_OPTIONS"] = "-nogui"
    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    if args.mpi_shards:
        from .mpi import run
    else:
        from .simulation import run

    run(args)


def main(argv=None):
    run_args(parse_args(argv))
