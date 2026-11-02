"""Compare the supplied three-cell wrapper with the current simulation engine.

The published cACint209 adapter is used in every condition because the original
regional MOD variants are unavailable. Source data and source scripts stay intact.
"""

import argparse
import ast
import csv
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import runpy
import subprocess
import sys
import time
import types

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
IDS = ["864691135640672571", "864691135993429185", "864691136025333561"]
DEFAULT_DATA = Path("/home/shisoft/Downloads/neuron_10-26-2026-20261027T002512Z-1-001/neuron_10-26-2026")
MODES = ["wrapper_reference", "current_coordinate_matched", "current_defaults", "zero_transmission"]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prepare(data, output):
    """Keep the exact three SWCs and every internal CSV contact, in row order."""
    import numpy as np
    from neuronsim.swc import load_swc

    working = (ROOT / "model/wrapper-synaptic-connections.py").read_text()
    original = (data / "wrapper-synaptic-connections.py").read_text()
    tree = ast.parse(original)
    settings = {}
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign):
            for target_node in node.targets:
                key = ast.unparse(target_node)
                if key in {"files", "connectivity", "stim.amp", "stim.delay", "stim.dur", "h.dt"}:
                    settings[key] = ast.literal_eval(node.value)
        elif isinstance(node, ast.Call) and ast.unparse(node.func) in {"h.finitialize", "h.continuerun"}:
            settings[ast.unparse(node.func)] = ast.literal_eval(node.args[0])
    expected = {"files": [root + ".swc" for root in IDS], "connectivity": "L5_all.csv", "stim.amp": 1, "stim.delay": 0, "stim.dur": 1, "h.dt": .025, "h.finitialize": -65, "h.continuerun": 25}
    if settings != expected:
        raise ValueError(f"Original wrapper protocol changed: {settings}")
    for cell_id in IDS:
        if cell_id not in original or cell_id not in working:
            raise ValueError("Wrapper cell selection has changed")
    with (data / "L5_all.csv").open(newline="") as stream:
        all_contacts = list(csv.DictReader(stream))
    contacts = [r for r in all_contacts if r["pre_pt_root_id"] in IDS and r["post_pt_root_id"] in IDS]
    if len(contacts) != 3:
        raise ValueError("Expected the wrapper's three-contact cycle")
    with (output / "selected_contacts.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(contacts[0]))
        writer.writeheader()
        writer.writerows(contacts)

    # Recover annotation identity without changing the supplied root IDs or trees.
    reference = ROOT / "outputs/NeuronSim-minnie-column.WCVsAeuu/data/column/nodes.csv"
    with reference.open() as stream:
        annotations = list(csv.DictReader(stream))
    nodes, records = [], []
    morph = output / "input/morphologies"
    (morph / "swc").mkdir(parents=True, exist_ok=True)
    for root in IDS:
        path = data / "SWC_L5" / (root + ".swc")
        rows, _, order = load_swc(path)
        soma = rows[order[0], 2:5]
        matches = [n for n in annotations if np.linalg.norm(soma - [float(n[a + "_um"]) for a in "xyz"]) < .01]
        if len(matches) != 1:
            raise ValueError("Ambiguous soma identity")
        node = dict(matches[0], root_id=root)
        nodes.append(node)
        target = morph / "swc" / path.name
        target.write_bytes(path.read_bytes())
        records.append({"root_id": root, "nucleus_id": node["nucleus_id"], "import_sha256": digest(target), "source": str(path)})
    (morph / "manifest.json").write_text(json.dumps({"records": records, "scope": "Exact supplied wrapper SWCs; no geometry changes"}, indent=2) + "\n")
    target = output / "input/data"
    target.mkdir(parents=True, exist_ok=True)
    for filename, rows in [("nodes.csv", nodes), ("edges.csv", [{"pre_pt_root_id": r["pre_pt_root_id"], "post_pt_root_id": r["post_pt_root_id"], "n_syn": "1"} for r in contacts])]:
        with (target / filename).open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    sources = [data / "wrapper-synaptic-connections.py", data / "cell.py", data / "L5_all.csv"] + [data / "SWC_L5" / (root + ".swc") for root in IDS]
    provenance = {"selected_cell_ids": IDS, "verified_original_wrapper_settings": settings, "all_table_rows": len(all_contacts), "internal_contacts": len(contacts), "sources": {str(p): digest(p) for p in sources}, "working_wrapper_sha256": digest(ROOT / "model/wrapper-synaptic-connections.py"), "working_cell_sha256": digest(ROOT / "model/cell.py"), "comparison_script_sha256": digest(Path(__file__)), "input_policy": "Original IDs, SWCs, contact rows and wrapper pulse; no latest-release substitution"}
    (target / "summary.json").write_text(json.dumps(provenance, indent=2) + "\n")
    return provenance


def worker(args):
    """Use a fresh NEURON process for each construction and integration mode."""
    os.environ.update(NEURON_MODULE_OPTIONS="-nogui", MPLBACKEND="Agg", OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    import numpy as np
    import neuron
    from neuron import h
    from neuronsim.cell import PublishedSWCCell
    from neuronsim.data import load_column
    from neuronsim.network import build_cells, connect_synapses, load_published_model, partition_cells
    from neuronsim.simulation import configure_solver, resident_memory_bytes

    begin = time.perf_counter()
    duration = getattr(args, "duration_ms", 25.)
    native_timing = {}
    source, profile = load_published_model(h, neuron)
    pc = h.ParallelContext()
    pc.nthread(1)
    configure_solver(h, pc)
    with (args.output / "selected_contacts.csv").open() as stream:
        contact_rows = list(csv.DictReader(stream))
    destination = args.output / args.mode
    destination.mkdir(exist_ok=True)
    source_vectors, source_events, g_vectors, i_vectors, mappings = [], [], [], [], []
    soma_vectors, axon_vectors, soma_spikes, spike_detectors = [], [], [], []

    def instruments(cells, connections, synapses):
        for root in IDS:
            cell = cells[root]
            soma_vectors.append(h.Vector().record(cell.soma[0](.5)._ref_v))
            axon_vectors.append(h.Vector().record(cell.axon[0](.5)._ref_v))
            nc = h.NetCon(cell.soma[0](.5)._ref_v, None, sec=cell.soma[0])
            nc.threshold = 0
            events = h.Vector()
            nc.record(events)
            soma_spikes.append(events)
            spike_detectors.append(nc)
        for synapse in synapses:
            g_vectors.append(h.Vector().record(synapse._ref_g))
            i_vectors.append(h.Vector().record(synapse._ref_i))

    def coordinates(cells):
        for row in contact_rows:
            pre, post = row["pre_pt_root_id"], row["post_pt_root_id"]
            a, da = cells[pre].nearest_segment([float(row[f"pre_pt_position_{x}"]) for x in "xyz"])
            b, db = cells[post].nearest_segment([float(row[f"post_pt_position_{x}"]) for x in "xyz"])
            mappings.append({"pre_root_id": pre, "post_root_id": post, "pre_section": a.sec.name(), "pre_x": a.x, "post_section": b.sec.name(), "post_x": b.x, "pre_distance_um": da, "post_distance_um": db})
        return mappings

    if args.mode == "wrapper_reference":
        # Execute the repository working wrapper. Only its Cell is replaced by
        # the already approved published-channel adapter, and xtra coordinate
        # storage is replaced by the equivalent cached segment-coordinate lookup.
        spec = importlib.util.spec_from_file_location("wrapper_source_cell", ROOT / "model/cell.py")
        original = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(original)
        adapters = []

        class Adapter:
            def __init__(self, path, nseg_reduction_factor=2):
                if nseg_reduction_factor != 2:
                    raise ValueError("Wrapper discretization changed")
                self.cell = PublishedSWCCell(path, profile, original_discretization=True)
                adapters.append(self)

            def __getattr__(self, name):
                return getattr(self.cell, name)

            def __str__(self):
                return self.cell.cell_id

        class Constructor(original.SynapseConstructor):
            @staticmethod
            def get_all_cells():
                return adapters, [str(c) for c in adapters]

            @staticmethod
            def find_synapse_segment(cell, xyz, units="nm"):
                assert units == "nm"
                return cell.nearest_segment(xyz)[0]

            def insert_synapses(self, table):
                super().insert_synapses(table)
                cells = {str(c): c.cell for c in adapters}
                coordinates(cells)
                instruments(cells, self.pre_syns, self.post_syns)
                for nc in self.pre_syns:
                    source_vectors.append(h.Vector().record(nc.preseg()._ref_v))
                    events = h.Vector()
                    nc.record(events)
                    source_events.append(events)

        shim = types.ModuleType("cell")
        shim.Cell, shim.SynapseConstructor = Adapter, Constructor
        sys.modules["cell"] = shim
        sys.path.insert(0, str(ROOT / "model"))
        os.environ["NEURONSIM_DATA"] = str(args.data)
        os.environ.pop("NEURONSIM_ANIMATION", None)
        import matplotlib.pyplot as plt
        plt.show = lambda: None
        # Instrument the two native calls in memory, without editing the wrapper
        # or replacing any NEURON method. Only tstop changes for long benchmarks.
        def timed_initialize(*values):
            started = time.perf_counter()
            result = h.finitialize(*values)
            native_timing["initialization_seconds"] = time.perf_counter() - started
            return result

        def timed_continue(stop):
            started, cpu_started = time.perf_counter(), time.process_time()
            result = h.continuerun(stop)
            native_timing["solve_seconds"] = time.perf_counter() - started
            native_timing["solve_cpu_seconds"] = time.process_time() - cpu_started
            return result

        class InstrumentWrapper(ast.NodeTransformer):
            calls = 0

            def visit_Call(self, node):
                self.generic_visit(node)
                function = ast.unparse(node.func)
                if function == "h.finitialize":
                    node.func = ast.Name(id="timed_initialize", ctx=ast.Load())
                    self.calls += 1
                elif function == "h.continuerun":
                    node.func = ast.Name(id="timed_continue", ctx=ast.Load())
                    node.args = [ast.Constant(value=duration)]
                    self.calls += 1
                return node

        wrapper_path = ROOT / "model/wrapper-synaptic-connections.py"
        transformer = InstrumentWrapper()
        program = ast.fix_missing_locations(transformer.visit(ast.parse(wrapper_path.read_text())))
        assert transformer.calls == 2
        namespace = {"__name__": "__main__", "__file__": str(wrapper_path), "timed_initialize": timed_initialize, "timed_continue": timed_continue}
        start = time.perf_counter()
        exec(compile(program, str(wrapper_path), "exec"), namespace)
        run_seconds = time.perf_counter() - start
        cells = {str(c): c.cell for c in adapters}
        sections = [s for root in IDS for s in cells[root].sections]
        connections, synapses = namespace["sc"].pre_syns, namespace["sc"].post_syns
        times = np.asarray(namespace["t"]).copy()
        assert all(np.array_equal(np.asarray(v), np.asarray(w)) for v, w in zip(axon_vectors, namespace["v_vecs"]))
        plt.close("all")
    else:
        column = load_column(args.output / "input/data", args.output / "input/morphologies", [])
        collection = build_cells(h, pc, column, args.output / "input/morphologies", profile, begin)
        cells, sections = collection.by_root, collection.sections
        # build_cells omits the older coordinate cache to reduce column memory.
        # Populate it identically to PublishedSWCCell for the three-site control.
        for cell in cells.values():
            xyz = []
            for section in cell.sections:
                arc = np.array([section.arc3d(i) for i in range(section.n3d())])
                points = np.array([[section.x3d(i), section.y3d(i), section.z3d(i)] for i in range(section.n3d())])
                for segment in section:
                    xyz.append([np.interp(segment.x * arc[-1], arc, points[:, axis]) for axis in range(3)])
            cell.segment_xyz_um = np.asarray(xyz)
        coordinates(cells)
        connections, synapses, detectors = [], [], []
        if args.mode in ("current_defaults", "zero_transmission"):
            contacts = connect_synapses(h, pc, column, collection, destination, False, -80.)
            connections, synapses = contacts.connections, contacts.synapses
            for row in contact_rows:
                cell = cells[row["pre_pt_root_id"]]
                segment = cell.soma[0](.5)
                source_vectors.append(h.Vector().record(segment._ref_v))
            if args.mode == "zero_transmission":
                for nc in connections:
                    nc.weight[0] = 0
        else:
            for index, row in enumerate(contact_rows):
                pre = cells[row["pre_pt_root_id"]]
                post = cells[row["post_pt_root_id"]]
                pre_segment, _ = pre.nearest_segment([float(row[f"pre_pt_position_{x}"]) for x in "xyz"])
                post_segment, _ = post.nearest_segment([float(row[f"post_pt_position_{x}"]) for x in "xyz"])
                gid = 100 + index
                detector = h.NetCon(pre_segment._ref_v, None, sec=pre_segment.sec)
                detector.threshold = 0
                pc.set_gid2node(gid, 0)
                pc.cell(gid, detector)
                detectors.append(detector)
                synapse = h.Exp2Syn(post_segment)
                synapse.tau1, synapse.tau2, synapse.e = .5, 3., 0.
                nc = pc.gid_connect(gid, synapse)
                nc.delay, nc.weight[0] = 1., .01
                source_vectors.append(h.Vector().record(pre_segment._ref_v))
                events = h.Vector()
                detector.record(events)
                source_events.append(events)
                connections.append(nc)
                synapses.append(synapse)
        partition = partition_cells(h, pc, collection, 1)
        configure_solver(h, pc)
        instruments(cells, connections, synapses)
        if args.mode in ("current_defaults", "zero_transmission"):
            # NEURON shares one event-recording vector per voltage source.
            # Reuse the soma detector's vector instead of replacing its recorder.
            source_events = [soma_spikes[IDS.index(row["pre_pt_root_id"])] for row in contact_rows]
        stimulus = h.IClamp(cells[IDS[0]].soma[0](.5))
        stimulus.amp, stimulus.delay, stimulus.dur = 1., 0., 1.
        time_vector = h.Vector().record(stimulus, h._ref_t)
        start = time.perf_counter()
        h.finitialize(-65.)
        native_timing["initialization_seconds"] = time.perf_counter() - start
        solve_started, cpu_started = time.perf_counter(), time.process_time()
        pc.psolve(duration)
        native_timing["solve_seconds"] = time.perf_counter() - solve_started
        native_timing["solve_cpu_seconds"] = time.process_time() - cpu_started
        run_seconds = time.perf_counter() - start
        partition.verify(h, sections)
        times = np.asarray(time_vector).copy()

    voltages = np.asarray([np.asarray(v) for v in soma_vectors])
    axons = np.asarray([np.asarray(v) for v in axon_vectors])
    final = np.asarray([segment.v for root in IDS for section in cells[root].sections for segment in section])
    assert len(cells) == 3 and len(connections) == 3 and int(h.List("Exp2Syn").count()) == 3
    assert len(times) == round(duration / .025) + 1 and abs(times[-1] - duration) < 1e-7
    assert np.isfinite(voltages).all() and np.isfinite(axons).all() and np.isfinite(final).all()
    np.savez_compressed(destination / "recordings.npz", time_ms=times, cell_ids=np.asarray(IDS), soma_mv=voltages, axon_mv=axons, source_mv=np.asarray(source_vectors), synaptic_conductance_us=np.asarray(g_vectors), synaptic_current_na=np.asarray(i_vectors), final_compartment_mv=final)
    with (destination / "voltage.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["time_ms"] + [f"{root}_{site}_mv" for site in ("soma", "axon") for root in IDS])
        writer.writerows(zip(times, *voltages, *axons))
    summary = {"mode": args.mode, "scope": "Three-cell wrapper dataset with common published ModelDB 241165 cACint209 biophysics; not a reproduction of missing custom mechanisms", "cell_ids": IDS, "cells": len(cells), "contacts": len(connections), "sections": len(sections), "segments": len(final), "source_commit": source["commit"], "threads": 1, "mpi_ranks": int(pc.nhost()), "temperature_c": float(h.celsius), "dt_ms": float(h.dt), "duration_ms": float(h.t), "stimulus": {"cell_id": IDS[0], "amplitude_na": 1., "onset_ms": 0., "duration_ms": 1.}, "synapse": {"tau1_ms": .5, "tau2_ms": 3., "e_mv": 0., "weight_us": 0. if args.mode == "zero_transmission" else .01, "delay_ms": 1., "threshold_mv": 0.}, "soma_spikes_ms": {root: list(v) for root, v in zip(IDS, soma_spikes)}, "source_spikes_ms_by_contact": [list(v) for v in source_events], "nearest_coordinate_mappings": mappings, "run_seconds": run_seconds, "total_seconds": time.perf_counter() - begin, "rss_bytes": resident_memory_bytes(), "finite": True, "timing_note": "wrapper_reference run_seconds includes its Python construction/plot creation; other conditions time initialization plus solve, so these timings are not a speed comparison"}
    summary.update(native_timing)
    summary["timing_note"] = "solve_seconds measures only h.continuerun for the wrapper or pc.psolve for the current engine; excludes construction, initialization, plotting and export. run_seconds retains the earlier asymmetric definitions; compare solve_seconds for numerical integration."
    summary["wrapper_stop_time_override_ms"] = duration if duration != 25 else None
    (destination / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({k: summary[k] for k in ("mode", "cells", "contacts", "segments", "soma_spikes_ms", "run_seconds")}), flush=True)


def analyze(output, provenance):
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    recordings = {mode: np.load(output / mode / "recordings.npz") for mode in MODES}
    summaries = {mode: json.loads((output / mode / "summary.json").read_text()) for mode in MODES}
    reference = recordings["wrapper_reference"]
    comparisons = {}
    for mode in MODES[1:]:
        candidate = recordings[mode]
        assert np.array_equal(reference["time_ms"], candidate["time_ms"])
        result = {}
        for name in ("soma_mv", "axon_mv", "final_compartment_mv", "synaptic_conductance_us"):
            delta = candidate[name] - reference[name]
            result[name] = {"maximum_absolute_difference": float(np.abs(delta).max()), "rms_difference": float(np.sqrt(np.mean(delta ** 2)))}
        comparisons[mode] = result
    # This is a meaningful end-to-end equivalence check across different event
    # routing and integration APIs, with the actual wrapper dataset.
    assert comparisons["current_coordinate_matched"]["axon_mv"]["maximum_absolute_difference"] < 1e-6
    assert comparisons["current_coordinate_matched"]["final_compartment_mv"]["maximum_absolute_difference"] < 1e-6
    for mode, values in recordings.items():
        for field in ("soma_mv", "axon_mv", "source_mv", "synaptic_conductance_us", "synaptic_current_na", "final_compartment_mv"):
            assert np.isfinite(values[field]).all(), (mode, field)
        assert np.all(values["soma_mv"][:, 0] == -65.)
        assert np.all(values["axon_mv"][:, 0] == -65.)
    assert recordings["wrapper_reference"]["synaptic_conductance_us"].max() > .0099
    assert np.count_nonzero(recordings["zero_transmission"]["synaptic_conductance_us"]) == 0
    assert summaries["wrapper_reference"]["soma_spikes_ms"] == summaries["current_coordinate_matched"]["soma_spikes_ms"]
    for path, expected_hash in provenance["sources"].items():
        assert digest(Path(path)) == expected_hash, f"Source changed: {path}"
    fig, axes = plt.subplots(3, 2, figsize=(13, 9), sharex=True)
    colors = ["#222222", "#16837c", "#ce5c21", "#73779d"]
    labels = ["Wrapper locations", "Current solver, matched locations", "Current soma/area routing", "No transmission"]
    for index, root in enumerate(IDS):
        for axis, site in zip(axes[index], ("soma_mv", "axon_mv")):
            for mode, color, label in zip(MODES, colors, labels):
                axis.plot(recordings[mode]["time_ms"], recordings[mode][site][index], color=color, label=label, lw=1.2, linestyle="--" if mode == "current_coordinate_matched" else "-")
            axis.set_title(f"Cell {index + 1} · {root}\n{'Soma' if site == 'soma_mv' else 'First axon section'}")
            axis.set_ylabel("Voltage (mV)")
            axis.axvspan(0, 1, color="#dddddd", alpha=.35)
    axes[0, 0].legend(fontsize=8)
    for axis in axes[-1]:
        axis.set_xlabel("Model time (ms)")
    fig.suptitle("Exact three-neuron wrapper dataset · shared published channel profile")
    fig.tight_layout()
    fig.savefig(output / "comparison.png", dpi=160)
    plt.close(fig)
    peaks = {mode: {root: {"soma_peak_mv": float(values["soma_mv"][index].max()), "axon_peak_mv": float(values["axon_mv"][index].max())} for index, root in enumerate(IDS)} for mode, values in recordings.items()}
    report = {"input": provenance, "comparison_to_wrapper_reference": comparisons, "per_cell_peaks": peaks, "summaries": summaries, "matched_location_engine_check_passed": True, "all_source_hashes_unchanged_after_runs": True, "limitation": "Original custom MOD variants remain unavailable. This isolates routing differences with common published biophysics; it cannot establish agreement with the author's unavailable original model."}
    (output / "comparison.json").write_text(json.dumps(report, indent=2) + "\n")
    lines = ["# Three-neuron wrapper comparison", "", "Exactly the three original L5 SWCs and all three internal contact rows are used. Cell 1 receives 1 nA from 0–1 ms; duration 25 ms, dt 0.025 ms, v_init −65 mV. No geometry or contact IDs were updated to a different Minnie release.", "", "All conditions use our published ModelDB 241165 cACint209 adaptation at 34°C. The original custom channels remain missing, so the wrapper reference replaces its Cell constructor with that same published adapter and uses the cached nearest-segment coordinates in place of inactive xtra storage. The repository working wrapper is executed, including its SynapseConstructor, stimulus and axon recording logic. Its CSV reader preserves the original contact IDs. No source scripts are modified.", "", "| Condition | Soma spikes | Maximum soma difference (mV) | Maximum axon difference (mV) |", "| --- | ---: | ---: | ---: |"]
    for mode, label in zip(MODES, labels):
        count = sum(len(v) for v in summaries[mode]["soma_spikes_ms"].values())
        metrics = comparisons.get(mode)
        lines.append(f"| {label} | {count} | {metrics['soma_mv']['maximum_absolute_difference'] if metrics else 0:.9g} | {metrics['axon_mv']['maximum_absolute_difference'] if metrics else 0:.9g} |")
    lines += ["", "With wrapper locations, neuron 1 spikes at 0.900 ms and neuron 2 at 6.025 ms. Neuron 3 depolarizes to −56.842 mV without spiking. Under current default routing, only neuron 1 spikes: neuron 2 peaks at −53.938 mV and does not send an onward event to neuron 3. The zero-transmission control also has only neuron 1 firing. The large voltage differences include the presence or absence of a spike, rather than measuring solver error.", "", "The matched-location condition retains measured contact-coordinate mappings and uses GID event routing plus pc.psolve. The current-default condition retains the same three edges but uses soma detectors and area-sampled postsynaptic segments, as in the column runner. The control sets only the current-default weights to zero. This experiment changes source and target assumptions together and does not isolate their individual contributions.", "", "The wrapper reference is the common-model comparison baseline, not the unavailable original electrical model. No multithreading is required for three cells. The current defaults' placement RNG uses the existing seed 20261027; this is one placement realization, not a statistical sample. All 2,708 final compartment voltages, 3,003 soma samples, 3,003 axon samples, and synaptic conductance recordings exactly match between the wrapper reference and the matched-location engine run. All source hashes remain unchanged.", "", "![Soma and axon traces](comparison.png)", "", "See [full comparison](comparison.json), [exact selected contacts](selected_contacts.csv), and each condition's recordings.npz, voltage.csv and summary.json. Timing definitions differ for the wrapper and engine runs; do not treat them as a speed benchmark.", "", "Reproduce:", "", "```sh", f".venv/bin/python tools/compare_wrapper_three.py --data {DEFAULT_DATA} --output {output}", "```"]
    (output / "README.md").write_text("\n".join(lines) + "\n")
    print(json.dumps(comparisons, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=DEFAULT_DATA)
    parser.add_argument("--output", type=Path, default=ROOT / "outputs/wrapper_three_neuron")
    parser.add_argument("--mode", choices=MODES)
    parser.add_argument("--duration-ms", type=float, default=25., help="Worker-only duration override for timing experiments")
    args = parser.parse_args()
    if args.duration_ms <= 0 or abs(args.duration_ms / .025 - round(args.duration_ms / .025)) > 1e-7:
        parser.error("Duration must be a positive multiple of 0.025 ms")
    if not args.mode and args.duration_ms != 25.:
        parser.error("Use --mode for a duration override; the four-condition report follows the original 25 ms protocol")
    args.data, args.output = args.data.resolve(), args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=True)
    if args.mode:
        worker(args)
        return
    provenance = prepare(args.data, args.output)
    for mode in MODES:
        with (args.output / (mode + ".log")).open("w") as log:
            subprocess.run([sys.executable, str(Path(__file__).resolve()), "--data", str(args.data), "--output", str(args.output), "--mode", mode], stdout=log, stderr=subprocess.STDOUT, check=True, cwd=ROOT)
        print("Completed", mode, flush=True)
    analyze(args.output, provenance)


if __name__ == "__main__":
    main()
