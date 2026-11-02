"""Validate and plot the NAS proofread-cohort run and its short controls."""

import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def validate_replay(main_run, replay, window):
    """Compare inclusive voltage samples and spikes strictly before the stop time.

    A longer run can publish a threshold event at the cutoff on its next step;
    a run stopped at that cutoff has no next step. Preserve/report those events
    separately instead of treating them as a different interior trajectory.
    """
    tolerance = 1e-7  # Far below the fixed 0.025-ms integration step.
    samples = main_run["time_ms"] <= window + tolerance
    for key, expected in (("root_ids", main_run["root_ids"]),
                          ("time_ms", main_run["time_ms"][samples]),
                          ("soma_mv", main_run["soma_mv"][:, samples])):
        require(np.array_equal(expected, replay[key]), f"Replay differs: {key}")
    main_spikes = main_run["spike_times_ms"] < window - tolerance
    replay_spikes = replay["spike_times_ms"] < window - tolerance
    for key in ("spike_times_ms", "spike_root_ids"):
        require(np.array_equal(main_run[key][main_spikes], replay[key][replay_spikes]),
                f"Replay differs before cutoff: {key}")
    closed_prefix = main_run["spike_times_ms"] <= window + tolerance
    return {
        "duration_ms": window,
        "all_saved_voltage_arrays_identical_to_main_prefix": True,
        "maximum_soma_voltage_difference_mv": 0.0,
        "spike_comparison_interval": f"[0, {window:g}) ms",
        "interior_spikes_compared": int(main_spikes.sum()),
        "all_interior_spike_times_and_ids_identical": True,
        "main_spikes_at_window_boundary": int((np.abs(main_run["spike_times_ms"] - window) <= tolerance).sum()),
        "replay_spikes_at_window_boundary": int((np.abs(replay["spike_times_ms"] - window) <= tolerance).sum()),
        "boundary_policy": "Cutoff events are reported separately because a longer run can emit them on its next step. Raw recordings are unchanged.",
        "all_saved_arrays_identical_to_main_prefix": all(
            np.array_equal(main_run[key][closed_prefix], replay[key]) for key in ("spike_times_ms", "spike_root_ids")
        ),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", type=Path, required=True)
    args = parser.parse_args()
    base = args.workspace.resolve()
    output = base / "results"
    prepared = json.loads((base / "prepared/data/summary.json").read_text())
    geometry = json.loads((base / "prepared/morphologies/manifest.json").read_text())
    labels = ("long_protocol", "validation_zero_transmission", "validation_replay")
    summaries = {label: json.loads((output / label / "summary.json").read_text()) for label in labels}
    arrays = {}
    for label in labels:
        with np.load(output / label / "recordings.npz") as archive:
            arrays[label] = {key: archive[key] for key in archive.files}
    main_run, control, replay = (arrays[label] for label in labels)
    summary = summaries[labels[0]]
    require(summary["neurons"] == prepared["selected_neurons"] == geometry["ready"], "Cohort size differs")
    require(summary["directed_pairs"] == prepared["internal_directed_pairs"], "Pair count differs")
    for name, expected in summary["input_hashes"].items():
        require(hashlib.sha256((base / "prepared/data" / name).read_bytes()).hexdigest() == expected, "Prepared table changed")
    require(hashlib.sha256((base / "prepared/morphologies/manifest.json").read_bytes()).hexdigest() == summary["geometry_manifest_sha256"], "Geometry manifest changed")
    source_hashes = json.loads((base / "source_hashes.json").read_text())
    for name, expected in source_hashes.items():
        require(hashlib.sha256((base / name).read_bytes()).hexdigest() == expected, f"Code changed: {name}")
    window = summaries[labels[2]]["duration_ms"]
    replay_check = validate_replay(main_run, replay, window)
    for label, data in arrays.items():
        info = summaries[label]
        require(np.array_equal(main_run["root_ids"], data["root_ids"]), "Cell ordering changed")
        require(data["soma_mv"].shape == (info["neurons"], len(data["time_ms"])), "Recording shape differs")
        require(np.isfinite(data["soma_mv"]).all(), "Nonfinite voltage")
        require(abs(data["time_ms"][-1] - info["duration_ms"]) < 1e-7, "Incomplete recording")
        require(info["all_recorded_voltages_and_final_compartment_voltages_finite"], "Invalid compartment output")
        require(len(info["segments_per_thread"]) == info["threads"], "Thread count differs")
        require(sum(info["segments_per_thread"]) == info["segments"], "Partition size differs")
        require(info["anatomical_contact_count"] == prepared["internal_anatomical_synapses"] == info["synaptic_receptor_objects"], "Contact count differs")
        require(info["reconstruction_materialization"] == info["connectivity_materialization"] == 1300, "Materialization differs")
        require(sum(info["contacts_by_presynaptic_anatomical_class"].values()) == info["anatomical_contact_count"], "Polarity counts differ")
        for key in ("input_hashes", "geometry_manifest_sha256", "source_commit", "profile_sha256", "stimulus", "dt_ms", "conductance_us", "synaptic_reversal_mv", "segments_per_thread"):
            require(info[key] == summary[key], f"Control setup changed: {key}")
    require(summaries[labels[1]]["weight_factor"] == 0, "Control transmits")
    require(summary["weight_factor"] == summaries[labels[2]]["weight_factor"] == 1, "Coupled weights differ")
    require(np.array_equal(control["time_ms"], replay["time_ms"]), "Control axes differ")
    for record in geometry["records"]:
        require(hashlib.sha256(Path(record["source_path"]).read_bytes()).hexdigest() == record["raw_sha256"], "NAS original changed")
        path = base / "prepared/morphologies/swc" / (record["root_id"] + ".swc")
        require(hashlib.sha256(path.read_bytes()).hexdigest() == record["import_sha256"], "Prepared SWC changed")
    with (output / "cell_models.csv").open() as stream:
        cells = {row["root_id"]: row for row in csv.DictReader(stream)}
    stimulated = set(summary["stimulus"]["cells"])
    index = {root: i for i, root in enumerate(main_run["root_ids"])}
    difference = np.max(np.abs(replay["soma_mv"] - control["soma_mv"]), axis=1)
    nonstim = np.array([root not in stimulated for root in main_run["root_ids"]])
    example = int(np.argmax(np.where(nonstim, difference, -np.inf)))
    result = {
        "cohort": {key: prepared[key] for key in ("selected_neurons", "ei_classes", "internal_directed_pairs", "internal_anatomical_synapses", "full_metadata_eligible_catalog")},
        "raw_and_prepared_swc_hashes_verified_after_simulation": len(geometry["records"]),
        "simulation_source_files_verified_after_simulation": len(source_hashes),
        "replay": replay_check,
        "control": {"duration_ms": window, "nonstimulated_cells_with_voltage_difference_above_1uv": int(((difference > .001) & nonstim).sum()), "largest_nonstimulated_voltage_difference_mv": float(difference[example]), "example_root_id": str(main_run["root_ids"][example])},
        "protocols": {},
    }
    distributed = summary["mpi_ranks"] > 1
    if distributed:
        result["memory_measurement_note"] = summary["resident_ram_measurement"]
        result["launcher_resource_note"] = "GNU time for mpiexec does not measure the aggregate peak memory of all MPI ranks. Its elapsed time covers MPI execution, before the separate merge and analysis."
    for line in (output / "process_resources.txt").read_text().splitlines():
        if "Maximum resident set size (kbytes):" in line:
            key = "mpi_launcher_reported_max_rss_bytes" if distributed else "whole_process_peak_resident_bytes"
            result[key] = int(line.rsplit(":", 1)[1]) * 1024
        if "Elapsed (wall clock) time" in line:
            fields = line.split("): ", 1)[1].split(":")
            key = "mpi_execution_wall_seconds" if distributed else "whole_process_wall_seconds"
            result[key] = sum(float(value) * 60 ** power for power, value in enumerate(reversed(fields)))
    for label, data in arrays.items():
        info = summaries[label]
        times, roots = data["spike_times_ms"], data["spike_root_ids"]
        result["protocols"][label] = {
            **{key: info[key] for key in ("duration_ms", "initialization_seconds", "solve_seconds", "recorded_neuronal_spikes", "cells_with_spikes", "nonstimulated_cell_spikes", "resident_ram_bytes_after_simulation")},
            "nonstimulated_spiking_cells": len(set(roots.tolist()) - stimulated),
            "first_spike_ms": float(times.min()) if len(times) else None,
            "last_spike_ms": float(times.max()) if len(times) else None,
            "spikes_after_200ms": int((times > 200).sum()),
            "final_somas_above_minus40_mv": int((data["soma_mv"][:, -1] > -40).sum()),
        }
    fig, axes = plt.subplots(4, 1, figsize=(12, 12), layout="constrained")
    for kind, color in (("E", "#2568a0"), ("I", "#bd4b31")):
        select = np.array([cells[root]["ei_class"] == kind for root in main_run["spike_root_ids"]])
        axes[0].scatter(main_run["spike_times_ms"][select], [index[root] for root in main_run["spike_root_ids"][select]], s=6, color=color, label=kind)
    execution = (f"{summary['mpi_ranks']} MPI ranks × {summary['threads_per_rank']} native thread per rank"
                 if distributed else f"{summary['threads']} native threads")
    axes[0].set(ylabel="Neuron index", title=f"Minnie v1300: {summary['neurons']:,} cells, {summary['anatomical_contact_count']:,} contacts; {execution}")
    axes[0].legend()
    low, median, high = np.percentile(main_run["soma_mv"], [10, 50, 90], axis=0)
    axes[1].fill_between(main_run["time_ms"], low, high, alpha=.2, color="#2568a0", label="10th–90th percentile")
    axes[1].plot(main_run["time_ms"], median, color="#2568a0", label="Population median")
    axes[1].set(ylabel="Soma voltage (mV)")
    axes[1].legend()
    bins = np.arange(0, window + 2.5, 2.5)
    for data, label, color in ((replay, "Coupled", "#2568a0"), (control, "Zero transmission", "#a66b28")):
        counts, _ = np.histogram(data["spike_times_ms"], bins)
        axes[2].stairs(counts, bins, label=label, color=color)
        axes[3].plot(data["time_ms"], data["soma_mv"][example], label=label, color=color)
    axes[2].set(ylabel="Spikes per 2.5 ms", title="Short control: same current stimulus, synaptic weights disabled")
    axes[3].set(ylabel="Soma voltage (mV)", title=f"Unstimulated cell with largest control difference: {main_run['root_ids'][example]}")
    for i, ax in enumerate(axes):
        ax.axvspan(10, 60, color="#878787", alpha=.12)
        ax.set(xlabel="Model time (ms)", xlim=(0, summary["duration_ms"] if i < 2 else window))
        if i > 1:
            ax.legend()
    fig.savefig(output / "network_response.png", dpi=160)
    plt.close(fig)
    (output / "validation.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
