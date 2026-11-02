"""Prepare destination-owned MPI shards without changing contact order or RNG draws."""

import argparse
from collections import Counter
from contextlib import ExitStack
import csv
import heapq
import json
from pathlib import Path
import shutil

from .data import sha256


def prepare(data, morphologies, output, ranks, stimulus_reference=None):
    with (data / "nodes.csv").open() as stream:
        nodes = sorted(csv.DictReader(stream), key=lambda row: row["root_id"])
    if not 1 <= ranks <= len(nodes):
        raise ValueError("MPI rank count must be positive and no greater than the cell count")
    gids = {node["root_id"]: gid for gid, node in enumerate(nodes)}
    if len(gids) != len(nodes):
        raise ValueError("Duplicate root ID")
    geometry = json.loads((morphologies / "manifest.json").read_text())
    records = {row["root_id"]: row for row in geometry["records"]}
    heap = [(0, rank) for rank in range(ranks)]
    owners = {}
    for node in sorted(nodes, key=lambda row: (-int(records[row["root_id"]].get("nodes") or 1), row["root_id"])):
        load, rank = heapq.heappop(heap)
        root = node["root_id"]
        owners[root] = rank
        heapq.heappush(heap, (load + int(records[root].get("nodes") or 1), rank))
    output.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(data / "nodes.csv", output / "global_nodes.csv")
    outgoing = Counter()
    contacts = pairs = 0
    per_rank_contacts = [0] * ranks
    per_rank_pairs = [0] * ranks
    rank_dirs = [output / f"rank{rank:04d}" for rank in range(ranks)]
    for directory in rank_dirs:
        (directory / "data").mkdir(parents=True)
        (directory / "morphologies").mkdir()
        (directory / "morphologies/swc").symlink_to((morphologies / "swc").resolve(), target_is_directory=True)
    with ExitStack() as stack:
        stream = stack.enter_context((data / "edges.csv").open())
        reader = csv.DictReader(stream)
        fields = reader.fieldnames + ["global_contact_offset"]
        writers = []
        for directory in rank_dirs:
            stream = stack.enter_context((directory / "data/edges.csv").open("w", newline=""))
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            writers.append(writer)
        for edge in reader:
            pre, post, count = edge["pre_pt_root_id"], edge["post_pt_root_id"], int(edge["n_syn"])
            if pre not in gids or post not in gids or count <= 0:
                raise ValueError("Invalid edge endpoint or contact count")
            owner = owners[post]
            edge["global_contact_offset"] = contacts
            writers[owner].writerow(edge)
            per_rank_contacts[owner] += count
            per_rank_pairs[owner] += 1
            outgoing[pre] += count
            contacts += count
            pairs += 1
    source_summary = json.loads((data / "summary.json").read_text())
    if "internal_anatomical_synapses" in source_summary and source_summary["internal_anatomical_synapses"] != contacts:
        raise ValueError("Source contact count differs")
    if stimulus_reference:
        stimulated = json.loads(stimulus_reference.read_text())["stimulus"]["cells"]
    else:
        eligible = [node["root_id"] for node in nodes if node["ei_class"] == "E" and node["manual_cell_type"] == "5P-IT"]
        stimulated = sorted(eligible, key=lambda root: (-outgoing[root], root))[:32]
    if not stimulated or not set(stimulated).issubset(gids):
        raise ValueError("Missing stimulus IDs; supply a reference for cohorts without 5P-IT cells")
    for rank, directory in enumerate(rank_dirs):
        local_nodes = [node for node in nodes if owners[node["root_id"]] == rank]
        local_records = [records[node["root_id"]] for node in local_nodes]
        with (directory / "data/nodes.csv").open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(nodes[0]))
            writer.writeheader()
            writer.writerows(local_nodes)
        local_geometry = {
            **geometry, "records": local_records, "ready": len(local_nodes), "requested": len(local_nodes), "failures": [],
            "nodes": sum(row["nodes"] for row in local_records) if all("nodes" in row for row in local_records) else None,
            "with_recorded_axon": sum(row.get("has_axon", False) for row in local_records),
            "root_compartment_normalizations": sum(bool(row.get("root_compartment_normalization")) for row in local_records),
            "zero_radius_repairs": sum(len(row.get("repairs", [])) for row in local_records),
            "collapsed_identical_vertices": sum(len(row.get("collapsed_identical_vertices", [])) for row in local_records),
        }
        (directory / "morphologies/manifest.json").write_text(json.dumps(local_geometry) + "\n")
        (directory / "data/summary.json").write_text(json.dumps({"rank": rank, "ranks": ranks, "neurons": len(local_nodes), "contacts": per_rank_contacts[rank], "pairs": per_rank_pairs[rank], "source_data": str(data.resolve())}) + "\n")
    manifest = {
        "ranks": ranks, "neurons": len(nodes), "contacts": contacts, "directed_pairs": pairs,
        "ownership": "Greedily balance recorded SWC point counts across ranks; global GIDs remain sorted by root ID; each contact belongs to its postsynaptic cell",
        "random_stream": "Each contact consumes one probability-choice draw and one uniform-position draw. Global offsets advance the identical PCG64 stream across skipped contacts.",
        "stimulus": {"cells": stimulated}, "source_data": str(data.resolve()), "source_morphologies": str(morphologies.resolve()),
        "input_hashes": {name: sha256(data / name) for name in ("nodes.csv", "edges.csv", "summary.json")},
        "geometry_manifest_sha256": sha256(morphologies / "manifest.json"),
        "contacts_per_rank": per_rank_contacts, "pairs_per_rank": per_rank_pairs,
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--morphologies", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--ranks", type=int, required=True)
    parser.add_argument("--stimulus-reference", type=Path)
    args = parser.parse_args()
    print(json.dumps(prepare(args.data, args.morphologies, args.output, args.ranks, args.stimulus_reference), indent=2))


if __name__ == "__main__":
    main()
