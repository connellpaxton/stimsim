"""Merge completed MPI rank outputs into the existing population file format."""

import argparse
from collections import Counter
import csv
import gzip
import heapq
import json
from pathlib import Path

import numpy as np


def merge(run, output):
    complete = json.loads((run / "MPI_COMPLETE.json").read_text())
    shards = Path(complete["shards"])
    manifest = json.loads((shards / "manifest.json").read_text())
    ranks = [run / f"rank{rank:04d}" for rank in range(complete["ranks"])]
    output.mkdir(parents=True, exist_ok=False)
    with (shards / "global_nodes.csv").open() as stream:
        roots = sorted(row["root_id"] for row in csv.DictReader(stream))
    index = {root: i for i, root in enumerate(roots)}
    if len(index) != len(roots) or len(roots) != manifest["neurons"]:
        raise ValueError("Global neuron identities differ")
    protocols = sorted(p.name for p in ranks[0].iterdir() if (p / "summary.json").is_file())
    if not protocols:
        raise ValueError("No completed protocols")
    for protocol in protocols:
        destination = output / protocol
        destination.mkdir()
        summaries = [json.loads((rank / protocol / "summary.json").read_text()) for rank in ranks]
        first = summaries[0]
        shared_fields = ("duration_ms", "weight_factor", "source_commit", "profile_sha256", "global_input_hashes",
                         "global_geometry_manifest_sha256", "stimulus", "dt_ms", "conductance_us", "synaptic_reversal_mv",
                         "synaptic_tau1_ms", "synaptic_tau2_ms", "synaptic_delay_ms", "placement_seed", "threads", "mpi_ranks")
        for summary in summaries:
            if any(summary[key] != first[key] for key in shared_fields):
                raise ValueError("MPI rank configurations differ")
            if summary.get("dataset") != first.get("dataset"):
                raise ValueError("MPI rank dataset identities differ")
        with np.load(ranks[0] / protocol / "recordings.npz") as archive:
            times = archive["time_ms"]
        voltages = np.empty((len(roots), len(times)), dtype=np.float64)
        occupied = set()
        spike_times, spike_roots = [], []
        for rank in ranks:
            with np.load(rank / protocol / "recordings.npz") as archive:
                local_roots = archive["root_ids"]
                if occupied.intersection(local_roots) or not np.array_equal(times, archive["time_ms"]):
                    raise ValueError("Duplicate cells or inconsistent MPI recording axes")
                occupied.update(local_roots)
                voltages[[index[root] for root in local_roots]] = archive["soma_mv"]
                spike_times.append(archive["spike_times_ms"])
                spike_roots.append(archive["spike_root_ids"])
        if occupied != set(roots) or not np.isfinite(voltages).all():
            raise ValueError("Incomplete or nonfinite distributed recording")
        event_times, event_roots = np.concatenate(spike_times), np.concatenate(spike_roots)
        order = np.lexsort((event_roots, event_times))
        event_times, event_roots = event_times[order], event_roots[order]
        np.savez_compressed(destination / "recordings.npz", time_ms=times, soma_mv=voltages,
                            root_ids=np.asarray(roots), spike_times_ms=event_times, spike_root_ids=event_roots)
        summary = {**first, "neurons": len(roots), "threads": complete["ranks"] * first["threads"],
                   "threads_per_rank": first["threads"], "mpi_ranks": complete["ranks"],
                   "input_hashes": manifest["input_hashes"], "geometry_manifest_sha256": manifest["geometry_manifest_sha256"],
                   "resident_ram_measurement": "Sum of rank RSS snapshots; shared pages may be counted on multiple ranks; not a peak or unique physical-memory measurement",
                   "shared_construction": str(run),
                   "segments_per_thread": [load for info in summaries for load in info["segments_per_thread"]],
                   "construction_seconds": max(info["construction_seconds"] for info in summaries),
                   "cell_construction_seconds": max(info["cell_construction_seconds"] for info in summaries),
                   "initialization_seconds": max(info["initialization_seconds"] for info in summaries),
                   "solve_seconds": max(info["solve_seconds"] for info in summaries),
                   "all_recorded_voltages_and_final_compartment_voltages_finite": all(info["all_recorded_voltages_and_final_compartment_voltages_finite"] for info in summaries),
                   "cells_with_spikes": len(set(event_roots.tolist())), "directed_pairs": manifest["directed_pairs"]}
        for key in ("sections", "segments", "anatomical_contact_count", "synaptic_receptor_objects", "recorded_neuronal_spikes",
                    "nonstimulated_cell_spikes", "resident_ram_bytes_after_simulation", "trees_with_recorded_axon", "root_compartment_normalizations"):
            summary[key] = sum(info.get(key) or 0 for info in summaries)
        for key in ("contacts_by_presynaptic_anatomical_class", "zero_conductance_segment_instances_removed", "spikes_by_anatomical_ei_class"):
            counter = Counter()
            for info in summaries:
                counter.update(info[key])
            summary[key] = dict(counter)
        summary.pop("mpi_rank", None)
        if summary["anatomical_contact_count"] != manifest["contacts"] or summary["recorded_neuronal_spikes"] != len(event_times):
            raise ValueError("Merged contact or spike counts differ")
        (destination / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
        if first.get("dataset"):
            (output / "dataset.json").write_text(json.dumps(first["dataset"], indent=2) + "\n")
    streams = [gzip.open(rank / "modeled_sites.csv.gz", "rt") for rank in ranks]
    try:
        readers = [csv.DictReader(stream) for stream in streams]
        fields = [key for key in readers[0].fieldnames if key != "global_contact_index"]
        with gzip.open(output / "modeled_sites.csv.gz", "wt", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            count = 0
            for row in heapq.merge(*readers, key=lambda item: int(item["global_contact_index"])):
                if int(row.pop("global_contact_index")) != count:
                    raise ValueError("MPI contact ordering has a gap or duplicate")
                writer.writerow(row)
                count += 1
            if count != manifest["contacts"]:
                raise ValueError("Incomplete modeled-site table")
    finally:
        for stream in streams:
            stream.close()
    assignments = []
    for rank in ranks:
        with (rank / "cell_models.csv").open() as stream:
            assignments.extend(csv.DictReader(stream))
    with (output / "cell_models.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(assignments[0]))
        writer.writeheader()
        writer.writerows(sorted(assignments, key=lambda row: row["root_id"]))
    (output / "MPI_MERGE_COMPLETE.json").write_text(json.dumps({**complete, "protocols": protocols, "run": str(run)}) + "\n")
    return {"neurons": len(roots), "contacts": manifest["contacts"], "protocols": protocols}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(merge(args.run.resolve(), args.output.resolve()), indent=2))


if __name__ == "__main__":
    main()
