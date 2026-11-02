"""Prepare a version-consistent neuron cohort from the downloaded NAS data."""

import argparse
import base64
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
import csv
import hashlib
import json
from pathlib import Path
import shutil
import sqlite3
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))
import numpy as np
from neuronsim.swc import load_swc
from prepare_minnie_skeletons import collapse_identical_vertices, repair_zero_radii


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_csv(path, rows):
    if not rows:
        raise ValueError("Refusing to prepare an empty table")
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def prepare_cell(node, source, destination, source_record):
    """Keep measured coordinates/radii; use the column's documented cleanup."""
    raw = source.read_bytes()
    if hashlib.md5(raw).digest() != base64.b64decode(source_record["md5Hash"]):
        raise ValueError("NAS source no longer matches its published MD5")
    original, parent, order = load_swc(source, require_positive_radii=False)
    if not set(original[:, 1]).issubset({1, 2, 3, 4}):
        raise ValueError("Unsupported SWC compartment label")
    soma = original[order[0]]
    distance = float(np.linalg.norm(soma[2:5] - [float(node[a + "_um"]) for a in "xyz"]))
    if distance > .01:
        raise ValueError(f"Soma identity mismatch: {distance} um")
    normalization = None
    working = original
    if soma[1] != 1:
        if (source_record["source_kind"] != "dendrite" or soma[1] != 3
                or np.any(original[:, 1] == 1) or soma[5] <= 0):
            raise ValueError("SWC is not rooted at a recorded soma")
        # The dendrite export labels its nucleus-positioned root as dendrite.
        # Use the verified root and its supplied radius as the soma compartment.
        working = original.copy()
        working[order[0], 1] = 1
        normalization = {"node_id": int(soma[0]), "original_type": 3, "imported_type": 1,
                         "evidence": "Dendrite-export root matches this unique nucleus's recorded soma position within 0.01 um; source root radius retained"}
    original_edges = np.flatnonzero(parent >= 0)
    original_length = np.linalg.norm(original[original_edges, 2:5] - original[parent[original_edges], 2:5], axis=1).sum()
    rows, parent, order, collapsed = collapse_identical_vertices(working, parent, order)
    radii, repairs = repair_zero_radii(rows, parent, order)
    edges = np.flatnonzero(parent >= 0)
    length = np.linalg.norm(rows[edges, 2:5] - rows[parent[edges], 2:5], axis=1).sum()
    if not np.isclose(original_length, length, atol=1e-6, rtol=1e-12):
        raise RuntimeError("Cleanup changed cable length")
    if repairs or collapsed or normalization:
        parents = {int(row[0]): str(int(row[6])) for row in rows}
        changed = {int(rows[i, 0]): repr(float(radii[i])) for i in np.flatnonzero(rows[:, 5] == 0)}
        lines = []
        for line in raw.decode().splitlines():
            fields = line.split()
            if fields and not line.lstrip().startswith("#"):
                ident = int(fields[0])
                if ident not in parents:
                    continue
                fields[6] = parents[ident]
                if normalization and ident == normalization["node_id"]:
                    fields[1] = "1"
                if ident in changed:
                    fields[5] = changed[ident]
                line = " ".join(fields)
            lines.append(line)
        destination.write_text("\n".join(lines) + "\n")
    else:
        shutil.copyfile(source, destination)
    imported, _, _ = load_swc(destination)
    if not np.array_equal(imported[:, [0, 1, 2, 3, 4, 6]], rows[:, [0, 1, 2, 3, 4, 6]]):
        raise RuntimeError("Cleanup changed retained geometry or topology")
    return {"root_id": node["root_id"], "nucleus_id": node["nucleus_id"], "source_path": str(source), "source_url": source_record["url"], "source_generation": source_record["generation"], "source_kind": source_record["source_kind"], "raw_sha256": hashlib.sha256(raw).hexdigest(), "import_sha256": sha(destination), "raw_nodes": len(original), "nodes": len(rows), "types": dict(Counter(int(t) for t in rows[:, 1])), "has_axon": bool(np.any(rows[:, 1] == 2)), "repairs": repairs, "collapsed_identical_vertices": collapsed, "soma_nucleus_distance_um": distance, "root_compartment_normalization": normalization, "cable_length_um": float(length)}


def prepare_candidate(task):
    """Validate one independent SWC; report exclusions without changing its source."""
    node, source, destination, source_record = task
    try:
        record = prepare_cell(node, source, destination, source_record)
    except Exception as exc:
        return node, None, {"root_id": node["root_id"], "nucleus_id": node["nucleus_id"], "reason": f"{type(exc).__name__}: {exc}"}
    return node, record, None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=Path("/mnt/nas-storage/NeuronSim/datasets/minnie65/v1300"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--cohort", choices=("proofread", "all"), default="proofread", help="Proofread export only, or every uniquely identified E/I cell with an SWC")
    parser.add_argument("--workers", type=int, default=1, help="Independent morphology preparation processes; simulation threading is separate")
    parser.add_argument("--connections-only", action="store_true", help="Reuse checksum-verified prepared morphologies after an incomplete wiring pass")
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("Worker count must be positive")
    source, output = args.source.resolve(), args.output.resolve()
    if json.loads((source / "DOWNLOAD_STATUS.json").read_text())["status"] != "complete":
        raise RuntimeError("NAS source download has not completed")
    if (output / "data/summary.json").exists():
        raise RuntimeError("Prepared output already exists; use a new experiment directory")
    data, morphology = output / "data", output / "morphologies"
    data.mkdir(parents=True, exist_ok=True)
    (morphology / "swc").mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    connection = sqlite3.connect("file:" + str(source / "metadata/cells.sqlite") + "?mode=ro", uri=True)
    connection.row_factory = sqlite3.Row
    eligible = [dict(row) for row in connection.execute("SELECT * FROM cells WHERE swc_available_in_catalog=1 AND root_nucleus_multiplicity=1 AND broad_type IN ('excitatory','inhibitory') ORDER BY root_id")]
    connection.close()
    candidates = eligible if args.cohort == "all" else [row for row in eligible if row["swc_source_kind"] == "proofread"]
    catalog = {row["root_id"]: row for row in (json.loads(line) for line in (source / "SOURCE_CATALOG.jsonl").open()) if row.get("root_id")}
    if args.connections_only:
        geometry = json.loads((morphology / "manifest.json").read_text())
        if geometry["requested"] != len(candidates):
            raise RuntimeError("Cached morphology cohort differs")
        with (data / "nodes.csv").open() as stream:
            nodes = list(csv.DictReader(stream))
        records, skipped = geometry["records"], geometry["failures"]
        expected = {row["root_id"]: row["nucleus_id"] for row in candidates}
        if len(nodes) != len(records) or len(nodes) != geometry["ready"]:
            raise RuntimeError("Cached morphology count differs")
        for record in records:
            path = morphology / "swc" / (record["root_id"] + ".swc")
            if expected.get(record["root_id"]) != record["nucleus_id"] or sha(path) != record["import_sha256"]:
                raise RuntimeError("Cached morphology identity or checksum differs")
    else:
        records, nodes, skipped = [], [], []
        tasks = []
        for row in candidates:
            node = {"root_id": row["root_id"], "nucleus_id": row["nucleus_id"], "manual_cell_type": row["cell_type"] or row["mtype"], "ei_class": "E" if row["broad_type"] == "excitatory" else "I", **{a + "_um": row[a + "_um"] for a in "xyz"}, "cell_type_source": row["cell_type_source"], "visual_area": row["visual_area"]}
            path = source / row["swc_path"]
            tasks.append((node, path, morphology / "swc" / path.name, catalog[node["root_id"]]))
        pool = ProcessPoolExecutor(max_workers=args.workers) if args.workers > 1 else None
        results = pool.map(prepare_candidate, tasks, chunksize=16) if pool else map(prepare_candidate, tasks)
        try:
            for index, (node, record, error) in enumerate(results, 1):
                if error:
                    skipped.append(error)
                else:
                    records.append(record)
                    nodes.append(node)
                if index % 100 == 0 or index == len(candidates):
                    print(f"Prepared {index}/{len(candidates)} {args.cohort} cohort cells; retained {len(nodes)}, excluded {len(skipped)}; {time.perf_counter()-started:.1f} s", flush=True)
                    (output / "STATUS.json").write_text(json.dumps({"status": "preparing_morphologies", "cohort": args.cohort, "requested": len(candidates), "processed": index, "retained": len(nodes), "excluded": len(skipped), "elapsed_seconds": time.perf_counter()-started}) + "\n")
        finally:
            if pool:
                pool.shutdown()
        if not nodes:
            raise RuntimeError("No valid neurons remain")
        write_csv(data / "nodes.csv", nodes)
        geometry = {"reconstruction_version": 1300, "wiring_version": 1300, "link_method": "identical v1300 root ID, unique stable nucleus ID and soma position agreement within 0.01 um", "requested": len(candidates), "ready": len(records), "failures": skipped, "records": records, "nodes": sum(r["nodes"] for r in records), "with_recorded_axon": sum(r["has_axon"] for r in records), "root_compartment_normalizations": sum(r["root_compartment_normalization"] is not None for r in records), "zero_radius_repairs": sum(len(r["repairs"]) for r in records), "collapsed_identical_vertices": sum(len(r["collapsed_identical_vertices"]) for r in records), "repair_policy": "Only zero radii inherit the nearest measured positive upstream radius. Identical parent/child aliases are collapsed without cable-length changes. A dendrite-export root may be relabeled from type 3 to soma type 1 only after unique-nucleus position agreement and a positive supplied root radius. Retained coordinates and positive radii are unchanged. NAS originals remain untouched."}
        (morphology / "manifest.json").write_text(json.dumps(geometry, indent=2) + "\n")
    (output / "STATUS.json").write_text(json.dumps({"status": "filtering_connections", "cohort": args.cohort, "retained": len(nodes), "excluded": len(skipped)}) + "\n")
    import pandas as pd
    table = source / "connectivity/connections_with_nuclei.csv.gz"
    header = [row[0] for row in csv.reader((source / "connectivity/connections_with_nuclei_header.csv").open())]
    expected = ["pre_pt_root_id", "post_pt_root_id", "n_syn", "sum_size", "pre_nuc_id", "post_nuc_id"]
    if header != expected:
        raise ValueError("Unexpected connection table schema")
    ids = {int(row["root_id"]) for row in nodes}
    all_ids = {int(row["root_id"]) for row in eligible}
    nuclei = {int(row["root_id"]): int(row["nucleus_id"]) for row in nodes}
    total_rows = pairs = contacts = broad_pairs = broad_contacts = 0
    missing_nucleus_pairs = missing_nucleus_contacts = 0
    with (data / "edges.csv").open("w", newline="") as stream:
        for chunk in pd.read_csv(table, header=None, names=header, dtype={name: "float64" if name == "sum_size" else "int64" for name in header}, chunksize=2_000_000):
            broad = chunk.pre_pt_root_id.isin(all_ids) & chunk.post_pt_root_id.isin(all_ids)
            broad_pairs += int(broad.sum())
            broad_contacts += int(chunk.loc[broad, "n_syn"].sum())
            selected = chunk.loc[chunk.pre_pt_root_id.isin(ids) & chunk.post_pt_root_id.isin(ids)]
            if len(selected):
                if (selected.n_syn <= 0).any():
                    raise RuntimeError("Nonpositive contact count")
                # -1 has no nucleus annotation in this aggregate export. Root
                # identity is independently verified against metadata and SWC.
                pre_match = (selected.pre_nuc_id == -1) | (selected.pre_nuc_id == selected.pre_pt_root_id.map(nuclei))
                post_match = (selected.post_nuc_id == -1) | (selected.post_nuc_id == selected.post_pt_root_id.map(nuclei))
                if not pre_match.all() or not post_match.all():
                    selected.loc[~(pre_match & post_match)].to_csv(data / "nucleus_identity_conflicts.csv", index=False)
                    raise RuntimeError("Connection endpoint has a conflicting annotated nucleus identity")
                missing = (selected.pre_nuc_id == -1) | (selected.post_nuc_id == -1)
                missing_nucleus_pairs += int(missing.sum())
                missing_nucleus_contacts += int(selected.loc[missing, "n_syn"].sum())
                selected.to_csv(stream, index=False, header=pairs == 0)
                pairs += len(selected)
                contacts += int(selected.n_syn.sum())
            total_rows += len(chunk)
            print(f"Scanned {total_rows:,} pairs; {args.cohort} cohort {pairs:,} pairs / {contacts:,} contacts; {time.perf_counter()-started:.1f} s", flush=True)
    if pairs == 0:
        raise RuntimeError("Selected cohort has no observed internal connections")
    summary = {"source": "MICrONS Minnie65 v1300 public static NAS exports", "materialization_version": 1300, "source_root": str(source), "selection": f"Every uniquely identified E/I neuron in {'the selected SWC catalog' if args.cohort == 'all' else 'the proofread SWC export'}, followed by documented morphology/identity checks", "cohort": args.cohort, "preparation_workers": args.workers, "swc_source_kinds": dict(Counter(r["source_kind"] for r in records)), "candidate_neurons": len(candidates), "selected_neurons": len(nodes), "excluded_neurons": skipped, "ei_classes": dict(Counter(n["ei_class"] for n in nodes)), "cell_types": dict(Counter(n["manual_cell_type"] for n in nodes)), "scanned_cell_pairs": total_rows, "connection_identity_validation": {"pairs_with_missing_nucleus_annotation": missing_nucleus_pairs, "contacts_with_missing_nucleus_annotation": missing_nucleus_contacts, "policy": "Only -1 is treated as unannotated. Endpoints match the selected v1300 root IDs with independent unique-nucleus/soma validation. Every nonmissing exported nucleus ID must match metadata. Original edge annotations are retained."}, "internal_directed_pairs": pairs, "internal_anatomical_synapses": contacts, "full_metadata_eligible_catalog": {"neurons": len(eligible), "internal_directed_pairs": broad_pairs, "anatomical_contacts": broad_contacts, "note": "Metadata eligibility before morphology checks; see selected cohort counts and subsequent simulation results."}, "source_hashes": {"source_catalog": sha(source / "SOURCE_CATALOG.jsonl"), "cell_index": sha(source / "metadata/cells.sqlite"), "connection_table": sha(table)}, "preparation_seconds": time.perf_counter() - started, "synapse_positions": "Same as the column model: somatic GID source and area-sampled postsynaptic target; the individual-contact coordinate table is retained on NAS for subsequent experiments."}
    (data / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    (output / "STATUS.json").write_text(json.dumps({"status": "complete", "cohort": args.cohort, "retained": len(nodes), "excluded": len(skipped), "pairs": pairs, "contacts": contacts, "elapsed_seconds": time.perf_counter()-started}) + "\n")
    print("PREPARATION COMPLETE", json.dumps({k: summary[k] for k in ["selected_neurons", "ei_classes", "internal_directed_pairs", "internal_anatomical_synapses", "full_metadata_eligible_catalog", "preparation_seconds"]}), flush=True)


if __name__ == "__main__":
    main()
