"""Load prepared Minnie identities, contacts, and reconstruction provenance."""

import csv
from dataclasses import dataclass
import hashlib
import json

from .graph import graph_summary


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


@dataclass
class ColumnData:
    nodes: list
    edges: list
    lookup: dict
    graph: dict
    skeletons: dict
    geometry_manifest: dict


def load_column(data, morphology, branch_roots):
    """Keep long root IDs as strings and check the selected recording cells."""
    with (data / "nodes.csv").open() as stream:
        nodes = list(csv.DictReader(stream))
    with (data / "edges.csv").open() as stream:
        edges = list(csv.DictReader(stream))
    lookup = {node["root_id"]: node for node in nodes}
    if len(lookup) != len(nodes):
        raise RuntimeError("Duplicate input neuron ID")
    if (
        not isinstance(branch_roots, list)
        or not all(isinstance(root, str) for root in branch_roots)
        or len(set(branch_roots)) != len(branch_roots)
        or not set(branch_roots).issubset(lookup)
    ):
        raise RuntimeError("Branch recording IDs must be a unique list of selected Minnie neurons")
    geometry = json.loads((morphology / "manifest.json").read_text())
    return ColumnData(
        nodes, edges, lookup, graph_summary(nodes, edges),
        {row["root_id"]: row for row in geometry["records"]}, geometry,
    )


def checked_swc(column, morphology, node):
    """Match nucleus identity and file checksum before importing a cell."""
    root = node["root_id"]
    path = morphology / "swc" / (root + ".swc")
    entry = column.skeletons[root]
    if entry["nucleus_id"] != node["nucleus_id"] or sha256(path) != entry["import_sha256"]:
        raise RuntimeError("Reconstruction identity or checksum mismatch")
    return path
