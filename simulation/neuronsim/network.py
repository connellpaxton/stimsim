"""Construct cells, synapses, and whole-cell native thread partitions."""

from collections import Counter
from dataclasses import dataclass
import csv
import gzip
import json
import time

from .config import PARAMETERS, PROFILE, SOURCE
from .data import checked_swc, sha256


def load_published_model(h, neuron):
    """Verify and load the unmodified ModelDB 241165 mechanism/profile files."""
    manifest = json.loads((SOURCE / "SOURCE_MANIFEST.json").read_text())
    for item in manifest["files"]:
        if sha256(SOURCE / item["path"]) != item["sha256"]:
            raise RuntimeError("Published reference source changed")
    neuron.load_mechanisms(str(SOURCE / "mechanisms"))
    h.load_file("stdrun.hoc")
    h.load_file("import3d.hoc")
    h.load_file(str(PROFILE))
    profile = h.cACint209_biophys()
    h.load_file(str(SOURCE / "cells/L4_LBC_cACint209_1/constants.hoc"))
    return manifest, profile


@dataclass
class Cells:
    """Keep Python references to native sections and spike detectors alive."""

    by_root: dict
    sections: list
    section_names: dict
    placement_choices: dict
    removed_zero_channels: Counter


def build_cells(h, pc, column, morphology, profile, started, global_gids=None):
    """Create one recorded tree per root ID and register its somatic spike GID."""
    import numpy as np
    from .cell import PublishedSWCCell

    cells = {}
    sections = []
    section_names = {}
    choices = {}
    removed = Counter()
    for index, node in enumerate(sorted(column.nodes, key=lambda row: row["root_id"])):
        root = node["root_id"]
        gid = global_gids[root] if global_gids is not None else index
        path = checked_swc(column, morphology, node)
        cell = PublishedSWCCell(
            path, profile, original_discretization=True, coordinate_cache=False
        )
        cell.gid = gid
        cells[root] = cell
        sections.extend(cell.sections)
        removed.update(cell.removed_zero_channels)

        names = {}
        for prefix in ("soma", "axon", "dend", "apic"):
            for index, section in enumerate(list(getattr(cell, prefix, []))):
                names[f"{prefix}_{index}"] = section
        section_names[root] = names
        dendrites = [
            (name, section) for name, section in names.items()
            if name.startswith(("dend_", "apic_"))
        ]
        available = dendrites or [("soma_0", cell.soma[0])]
        areas = np.array([
            sum(float(h.area(segment.x, sec=section)) for segment in section)
            for _, section in available
        ])
        choices[root] = ([name for name, _ in available], areas / areas.sum())

        # GIDs identify spike sources. A source spike schedules NetCon events;
        # Python does not scan every synapse at every integration step.
        pc.set_gid2node(gid, int(pc.id()))
        cell.detector = h.NetCon(
            cell.soma[0](0.5)._ref_v, None, sec=cell.soma[0]
        )
        cell.detector.threshold = 0.0
        pc.cell(gid, cell.detector)
        if (index + 1) % 25 == 0:
            elapsed = time.perf_counter() - started
            prefix = f"Rank {int(pc.id())}: " if int(pc.nhost()) > 1 else ""
            print(f"{prefix}Published regional cells {index + 1}/{len(column.nodes)}; elapsed {elapsed:.1f} s", flush=True)
    print(f"All {len(cells)} cells use the unchanged published cACint209 regional assignments", flush=True)
    return Cells(cells, sections, section_names, choices, removed)


def synaptic_reversal(node, inhibition=False, inhibitory_reversal=-80.0):
    """Presynaptic Minnie class selects polarity; weights remain positive."""
    if node["ei_class"] not in {"E", "I"}:
        raise ValueError("Synapse source requires an unambiguous E/I class")
    if inhibition and node["ei_class"] == "I":
        return float(inhibitory_reversal)
    return 0.0


@dataclass
class Contacts:
    synapses: list
    connections: list
    by_source_class: Counter


def connect_synapses(h, pc, column, cells, output, inhibition, inhibitory_reversal, global_gids=None):
    """Expand each observed n_syn count into that many native Exp2Syn objects."""
    import numpy as np

    settings = PARAMETERS
    rng = np.random.default_rng(settings.placement_seed)
    synapses = []
    connections = []
    counts = Counter()
    next_contact = 0
    with gzip.open(output / "modeled_sites.csv.gz", "wt", newline="") as stream:
        writer = csv.writer(stream)
        header = [
            "pre_root_id", "post_root_id", "contact_ordinal", "model_section",
            "model_x", "conductance_us", "reversal_mv", "delay_ms",
        ]
        writer.writerow(header + (["global_contact_index"] if global_gids is not None else []))
        # Preserve input row order and random-number call order for exact replay.
        for index, edge in enumerate(column.edges, 1):
            offset = int(edge.get("global_contact_offset", next_contact))
            if offset < next_contact:
                raise ValueError("Contact offsets must preserve global input order")
            if offset > next_contact:
                rng.bit_generator.advance(2 * (offset - next_contact))
            pre, post = edge["pre_pt_root_id"], edge["post_pt_root_id"]
            source = column.lookup[pre]
            reversal = synaptic_reversal(source, inhibition, inhibitory_reversal)
            counts[source["ei_class"]] += int(edge["n_syn"])
            names, probabilities = cells.placement_choices[post]
            for ordinal in range(int(edge["n_syn"])):
                name = names[int(rng.choice(len(names), p=probabilities))]
                position = float(rng.uniform(0.1, 0.9))
                synapse = h.Exp2Syn(cells.section_names[post][name](position))
                synapse.tau1 = settings.synapse_rise_ms
                synapse.tau2 = settings.synapse_decay_ms
                synapse.e = reversal
                source_gid = global_gids[pre] if global_gids is not None else cells.by_root[pre].gid
                connection = pc.gid_connect(source_gid, synapse)
                connection.delay = settings.delay_ms
                connection.weight[0] = settings.contact_conductance_us
                connections.append(connection)
                synapses.append(synapse)
                row = [
                    pre, post, ordinal, name, position,
                    settings.contact_conductance_us, reversal, settings.delay_ms,
                ]
                writer.writerow(row + ([offset + ordinal] if global_gids is not None else []))
            next_contact = offset + int(edge["n_syn"])
            if index % 10000 == 0:
                print(f"Connected {index:,}/{len(column.edges):,} observed pairs", flush=True)
    if int(h.List("Exp2Syn").count()) != column.graph["anatomical_contacts"]:
        raise RuntimeError("Synapse object count differs from observed contacts")
    return Contacts(synapses, connections, counts)


@dataclass
class ThreadPartition:
    roots: list
    loads: list

    def verify(self, h, sections):
        actual = [0] * len(self.loads)
        for section in sections:
            actual[int(h.original_family_thread(sec=section))] += section.nseg
        if actual != self.loads:
            raise RuntimeError("Native thread partition mismatch")
        return actual


def partition_cells(h, pc, cells, thread_count):
    """Greedily balance whole cells; never split a cell's cable tree."""
    partitions = [h.SectionList() for _ in range(thread_count)]
    loads = [0] * thread_count
    for cell in sorted(cells.by_root.values(), key=lambda c: (-len(c.segments), c.cell_id)):
        thread = min(range(thread_count), key=lambda index: loads[index])
        for section in cell.sections:
            if section.parentseg() is None:
                partitions[thread].append(sec=section)
        loads[thread] += len(cell.segments)
    pc.nthread(thread_count)
    for index, roots in enumerate(partitions):
        pc.partition(index, roots)

    # NEURON 9's Python binding coerces this result to bool. HOC preserves IDs
    # greater than one so verification measures the actual native partition.
    h("objref original_family_pc")
    h.original_family_pc = pc
    h("func original_family_thread() { return original_family_pc.sec_in_thread() }")
    return ThreadPartition(partitions, loads)


def select_stimulated_cells(cells, column, reference=None):
    """Reuse a reference's IDs, or select high-output Minnie layer-5 IT cells."""
    if reference:
        stimulated = json.loads(reference.resolve().read_text())["stimulus"]["cells"]
        if not set(stimulated).issubset(cells.by_root):
            raise RuntimeError("Reference stimulus neuron is absent")
    else:
        eligible = [
            root for root in cells.by_root
            if column.lookup[root]["ei_class"] == "E"
            and column.lookup[root]["manual_cell_type"] == "5P-IT"
        ]
        stimulated = sorted(
            eligible,
            key=lambda root: (-column.graph["outgoing_contacts"].get(root, 0), root),
        )[:32]
    if not stimulated:
        raise RuntimeError("No stimulus cell is available")
    return stimulated
