"""Compare merged MPI recordings and every modeled contact with a reference."""

import argparse
import csv
import gzip
from itertools import zip_longest
import json
from pathlib import Path

import numpy as np


def compare(run, reference, tolerance_mv=1e-8):
    result = {"run": str(run), "reference": str(reference), "voltage_tolerance_mv": tolerance_mv, "protocols": {}}
    with gzip.open(run / "modeled_sites.csv.gz", "rt") as left, gzip.open(reference / "modeled_sites.csv.gz", "rt") as right:
        count = 0
        for new, old in zip_longest(csv.DictReader(left), csv.DictReader(right)):
            if new != old:
                raise RuntimeError(f"Modeled contact differs at row {count}")
            count += 1
    result["identical_modeled_contacts"] = count
    for directory in sorted(run.iterdir()):
        if not (directory / "recordings.npz").exists():
            continue
        old_path = reference / directory.name / "recordings.npz"
        if not old_path.exists():
            continue
        with np.load(directory / "recordings.npz") as new, np.load(old_path) as old:
            stop = float(new["time_ms"][-1])
            mask = old["time_ms"] <= stop + 1e-7
            if not np.array_equal(new["time_ms"], old["time_ms"][mask]) or not np.array_equal(new["root_ids"], old["root_ids"]):
                raise RuntimeError("Recording axes differ")
            difference = float(np.max(np.abs(new["soma_mv"] - old["soma_mv"][:, mask])))
            spikes = old["spike_times_ms"] <= stop + 1e-7
            same_spikes = (np.array_equal(new["spike_root_ids"], old["spike_root_ids"][spikes])
                           and np.array_equal(new["spike_times_ms"], old["spike_times_ms"][spikes]))
            result["protocols"][directory.name] = {"duration_ms": stop, "soma_voltage_samples": int(new["soma_mv"].size),
                "maximum_voltage_difference_mv": difference, "identical_spike_times_and_ids": same_spikes,
                "spikes": len(new["spike_times_ms"])}
            if not np.isfinite(difference) or difference > tolerance_mv or not same_spikes:
                raise RuntimeError("Distributed simulation differs: " + json.dumps(result["protocols"][directory.name]))
    if not result["protocols"]:
        raise RuntimeError("No shared protocols were compared")
    (run / "reference_comparison.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--tolerance-mv", type=float, default=1e-8)
    args = parser.parse_args()
    print(json.dumps(compare(args.run.resolve(), args.reference.resolve(), args.tolerance_mv), indent=2))


if __name__ == "__main__":
    main()
