"""Orchestrate the current Minnie model and advance NEURON's native solver."""

from collections import Counter
import json
from pathlib import Path
import platform
import time

from .config import PARAMETERS, PROFILE, make_protocols
from .data import load_column, sha256
from .network import (
    build_cells,
    connect_synapses,
    load_published_model,
    partition_cells,
    select_stimulated_cells,
)


def resident_memory_bytes():
    """Linux resident-memory snapshot; this is not a peak-memory measurement."""
    return next(
        int(line.split()[1]) * 1024
        for line in Path("/proc/self/status").read_text().splitlines()
        if line.startswith("VmRSS:")
    )


def configure_solver(h, pc):
    """Use the benchmarked fixed timestep and native cache-efficient layout."""
    h.CVode().active(0)
    h.CVode().cache_efficient(1)
    h.dt = PARAMETERS.dt_ms
    h.steps_per_ms = 1.0 / PARAMETERS.dt_ms
    h.celsius = PARAMETERS.temperature_c
    pc.set_maxstep(PARAMETERS.delay_ms)


def make_summary(args, neuron, pc, column, cells, contacts, partition, source, started):
    """Export model assumptions and input provenance alongside every result."""
    p = PARAMETERS
    reconstruction_version = column.geometry_manifest.get("reconstruction_version", "unspecified")
    wiring_version = column.geometry_manifest.get("wiring_version", "unspecified")
    polarity_assumption = (
        "Minnie presynaptic E/I labels set excitatory versus inhibitory reversal potentials; synaptic strengths, kinetics and delay remain assumed. Having both polarities does not establish E/I balance."
        if args.inhibition else
        "All contacts have e=0 mV, including contacts from anatomically inhibitory cells, exactly as the original SynapseConstructor default; this is not a balanced E/I network."
    )
    return {
        "scope": "Published ModelDB 241165 cACint209 regional adaptation of the supplied script on Minnie geometry; not an exact reproduction of missing custom MOD variants",
        "hostname": platform.node(), "neuron_version": neuron.__version__,
        "modeldb_id": 241165, "source_commit": source["commit"],
        "profile_sha256": sha256(PROFILE),
        "geometry_manifest_sha256": sha256(args.morphologies / "manifest.json"),
        "reconstruction_materialization": reconstruction_version,
        "connectivity_materialization": wiring_version,
        "trees_with_recorded_axon": column.geometry_manifest.get("with_recorded_axon"),
        "root_compartment_normalizations": column.geometry_manifest.get("root_compartment_normalizations", 0),
        "input_hashes": {
            name: sha256(args.data / name)
            for name in ("nodes.csv", "edges.csv", "summary.json")
        },
        "neurons": len(cells.by_root), "sections": len(cells.sections),
        "segments": sum(partition.loads),
        "anatomical_contact_count": column.graph["anatomical_contacts"],
        "synaptic_receptor_objects": len(contacts.synapses),
        "directed_pairs": len(column.edges), "threads": args.threads,
        "mpi_ranks": int(pc.nhost()),
        "temperature_c": p.temperature_c, "dt_ms": p.dt_ms,
        "v_init_mv": p.initial_voltage_mv, "conductance_us": p.contact_conductance_us,
        "inhibitory_outputs_enabled": args.inhibition,
        "contacts_by_presynaptic_anatomical_class": dict(contacts.by_source_class),
        "synaptic_reversal_mv": {
            "E": 0.0, "I": args.inhibitory_reversal if args.inhibition else 0.0,
        },
        "synaptic_tau1_ms": p.synapse_rise_ms, "synaptic_tau2_ms": p.synapse_decay_ms,
        "synaptic_delay_ms": p.delay_ms, "placement_seed": p.placement_seed,
        "discretization": "Original script d_lambda=0.01 at 100 Hz, followed by factor 3**2 reduction and odd segment count",
        "zero_conductance_segment_instances_removed": dict(cells.removed_zero_channels),
        "construction_seconds": time.perf_counter() - started,
        "assumptions": [
            "One published L4 basket-cell profile is applied to all trees to follow the original single-profile Cell design; this is not cell-type-specific V1 physiology.",
            polarity_assumption,
            "The aggregate graph lacks synaptic coordinates: postsynaptic sections are area sampled on recorded dendrites, and presynaptic events come from the soma rather than measured contact segments.",
            "Temperature is the published 34 C; the original wrapper did not explicitly set temperature.",
            "No extracellular stimulus is applied. Inactive xtra/extracellular mechanisms are omitted.",
            f"Recorded v{reconstruction_version} trees are matched to v{wiring_version} wiring with verified root/nucleus identity; no surrogate shapes or invented axons are used.",
        ],
    }


def run_protocol(h, pc, protocol, output, shared, column, cells, contacts,
                 partition, clamps, recording, branches, first):
    """Reset states, inject current, solve, and save one condition."""
    import numpy as np

    p = PARAMETERS
    destination = output / protocol.label
    destination.mkdir(parents=True, exist_ok=True)
    for root, clamp in clamps.items():
        clamp.amp = protocol.amplitude_na if root in protocol.cells else 0.0
        clamp.delay = protocol.onset_ms
        clamp.dur = protocol.pulse_ms
    for connection in contacts.connections:
        connection.weight[0] = p.contact_conductance_us * protocol.weight_factor
    recording.clear_spikes()

    initialization = time.perf_counter()
    h.finitialize(p.initial_voltage_mv)
    initialization = time.perf_counter() - initialization
    actual_loads = partition.verify(h, cells.sections)

    # NEURON advances all 0.025-ms steps inside C++ between these progress calls.
    if int(pc.nhost()) > 1:
        pc.barrier()
    solve = time.perf_counter()
    stop = 0.0
    while stop < protocol.duration_ms:
        stop = min(protocol.duration_ms, stop + p.progress_interval_ms)
        pc.psolve(stop)
        if int(pc.id()) == 0:
            print(f"{protocol.label}: {stop:g}/{protocol.duration_ms:g} ms; solve {time.perf_counter() - solve:.2f} s", flush=True)
    elapsed = time.perf_counter() - solve

    times, voltages, spikes, roots = recording.read(p.dt_ms, p.record_stride)
    final = np.fromiter(
        (segment.v for section in cells.sections for segment in section), dtype=float
    )
    if (
        not np.isfinite(voltages).all() or not np.isfinite(final).all()
        or abs(float(h.t) - protocol.duration_ms) > max(1e-7, protocol.duration_ms * 1e-10)
    ):
        raise RuntimeError("Invalid simulation output")
    result = {
        **shared,
        "construction_seconds": shared["construction_seconds"] if first else 0.0,
        "shared_construction": str(output), "stimulus": protocol.stimulus_metadata(),
        "duration_ms": protocol.duration_ms, "weight_factor": protocol.weight_factor,
        "initialization_seconds": initialization, "solve_seconds": elapsed,
        "actual_final_time_ms": float(h.t), "segments_per_thread": actual_loads,
        "recorded_neuronal_spikes": len(roots),
        "cells_with_spikes": len(set(roots.tolist())),
        "nonstimulated_cell_spikes": sum(root not in protocol.cells for root in roots),
        "spikes_by_anatomical_ei_class": dict(Counter(column.lookup[root]["ei_class"] for root in roots)),
        "resident_ram_bytes_after_simulation": resident_memory_bytes(),
        "all_recorded_voltages_and_final_compartment_voltages_finite": True,
    }
    recording.save(destination, times, voltages, spikes, roots)
    if branches:
        result["branch_recording"] = branches.save(destination, times, p.record_stride)
    (destination / "summary.json").write_text(json.dumps(result, indent=2) + "\n")
    if int(pc.id()) == 0:
        print(json.dumps({key: result[key] for key in (
            "solve_seconds", "recorded_neuronal_spikes", "nonstimulated_cell_spikes",
            "resident_ram_bytes_after_simulation",
        )}, indent=2), flush=True)


def run(args):
    """The complete execution order, kept visible for the code walkthrough."""
    import neuron
    from neuron import h
    from .recording import BranchRecording, PopulationRecording, write_cell_assignments

    started = time.perf_counter()
    source, profile = load_published_model(h, neuron)
    branch_roots = json.loads(args.branch_roots.read_text()) if args.branch_roots else []
    column = load_column(args.data, args.morphologies, branch_roots)
    pc = h.ParallelContext()
    pc.nthread(1)
    cells = build_cells(h, pc, column, args.morphologies, profile, started)
    contacts = connect_synapses(
        h, pc, column, cells, args.output, args.inhibition, args.inhibitory_reversal
    )
    partition = partition_cells(h, pc, cells, args.threads)
    configure_solver(h, pc)
    stimulated = select_stimulated_cells(cells, column, args.stimulus_reference)

    clamps = {}
    for root in stimulated:
        clamp = h.IClamp(cells.by_root[root].soma[0](0.5))
        clamp.amp = clamp.dur = 0.0
        clamps[root] = clamp
    recording = PopulationRecording(h, pc, cells.by_root)
    branches = None
    if branch_roots:
        branches = BranchRecording(
            cells.by_root, branch_roots, cells.section_names,
            args.output / "branch_geometry.json",
        )
    write_cell_assignments(
        args.output, cells.by_root, column.lookup, recording.ordered, stimulated
    )
    shared = make_summary(args, neuron, pc, column, cells, contacts, partition, source, started)
    print(json.dumps({key: shared[key] for key in (
        "neurons", "segments", "anatomical_contact_count", "construction_seconds"
    )}, indent=2), flush=True)
    for index, protocol in enumerate(make_protocols(args, stimulated)):
        run_protocol(
            h, pc, protocol, args.output, shared, column, cells, contacts,
            partition, clamps, recording, branches, first=index == 0,
        )
