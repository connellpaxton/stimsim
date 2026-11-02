"""Build and simulate independent cell subsets concurrently with NEURON MPI."""

from copy import copy
import csv
import json
import time


def run(args):
    # MPI must be initialized before the first NEURON ParallelContext.
    from mpi4py import MPI
    import neuron
    from neuron import h
    from .config import make_protocols
    from .data import ColumnData
    from .network import build_cells, connect_synapses, load_published_model, partition_cells
    from .recording import PopulationRecording, write_cell_assignments
    from .simulation import configure_solver, make_summary, run_protocol

    comm = MPI.COMM_WORLD
    rank, ranks = comm.Get_rank(), comm.Get_size()
    try:
        started = time.perf_counter()
        manifest = json.loads((args.mpi_shards / "manifest.json").read_text())
        if manifest["ranks"] != ranks:
            raise ValueError("Prepared shard count differs from MPI rank count")
        if args.branch_roots:
            raise ValueError("MPI branch recording is not implemented; omit --branch-roots")
        h.nrnmpi_init()
        pc = h.ParallelContext()
        if int(pc.id()) != rank or int(pc.nhost()) != ranks:
            raise RuntimeError("NEURON and mpi4py communicators differ")
        pc.nthread(1)
        source, profile = load_published_model(h, neuron)
        local = copy(args)
        directory = args.mpi_shards / f"rank{rank:04d}"
        local.data, local.morphologies = directory / "data", directory / "morphologies"
        local.output = args.output / f"rank{rank:04d}"
        local.output.mkdir(parents=True, exist_ok=True)
        with (args.mpi_shards / "global_nodes.csv").open() as stream:
            global_nodes = sorted(csv.DictReader(stream), key=lambda row: row["root_id"])
        gids = {node["root_id"]: gid for gid, node in enumerate(global_nodes)}
        with (local.data / "nodes.csv").open() as stream:
            nodes = list(csv.DictReader(stream))
        with (local.data / "edges.csv").open() as stream:
            edges = list(csv.DictReader(stream))
        geometry = json.loads((local.morphologies / "manifest.json").read_text())
        column = ColumnData(
            nodes, edges, {node["root_id"]: node for node in global_nodes},
            {"anatomical_contacts": sum(int(row["n_syn"]) for row in edges)},
            {row["root_id"]: row for row in geometry["records"]}, geometry,
        )
        cells = build_cells(h, pc, column, local.morphologies, profile, started, global_gids=gids)
        cell_seconds = time.perf_counter() - started
        contacts = connect_synapses(
            h, pc, column, cells, local.output, args.inhibition, args.inhibitory_reversal, global_gids=gids
        )
        if comm.allreduce(len(contacts.synapses), op=MPI.SUM) != manifest["contacts"]:
            raise RuntimeError("Distributed contact count differs from source graph")
        partition = partition_cells(h, pc, cells, args.threads)
        comm.Barrier()
        configure_solver(h, pc)
        stimulated = manifest["stimulus"]["cells"]
        clamps = {}
        for root in stimulated:
            if root not in cells.by_root:
                continue
            clamp = h.IClamp(cells.by_root[root].soma[0](0.5))
            clamp.amp = clamp.dur = 0.0
            clamps[root] = clamp
        recording = PopulationRecording(h, pc, cells.by_root)
        write_cell_assignments(local.output, cells.by_root, column.lookup, recording.ordered, stimulated)
        shared = make_summary(local, neuron, pc, column, cells, contacts, partition, source, started)
        shared["mpi_rank"] = rank
        shared["cell_construction_seconds"] = cell_seconds
        shared["global_input_hashes"] = manifest["input_hashes"]
        shared["global_geometry_manifest_sha256"] = manifest["geometry_manifest_sha256"]
        timing = comm.gather({"rank": rank, "cell_seconds": cell_seconds, "construction_seconds": shared["construction_seconds"],
                              "neurons": len(cells.by_root), "segments": sum(partition.loads), "contacts": len(contacts.synapses)}, root=0)
        if rank == 0:
            (args.output / "construction.json").write_text(json.dumps({"ranks": ranks, "threads_per_rank": args.threads, "rank_statistics": timing,
                "maximum_cell_construction_seconds": max(row["cell_seconds"] for row in timing),
                "maximum_total_construction_seconds": max(row["construction_seconds"] for row in timing)}, indent=2) + "\n")
            print(f"MPI construction complete: {manifest['neurons']:,} cells, {manifest['contacts']:,} contacts across {ranks} ranks", flush=True)
        for index, protocol in enumerate(make_protocols(args, stimulated)):
            run_protocol(h, pc, protocol, local.output, shared, column, cells, contacts,
                         partition, clamps, recording, None, first=index == 0)
            comm.Barrier()
        if rank == 0:
            (args.output / "MPI_COMPLETE.json").write_text(json.dumps({"status": "complete", "ranks": ranks, "threads_per_rank": args.threads,
                "neurons": manifest["neurons"], "contacts": manifest["contacts"], "shards": str(args.mpi_shards)}) + "\n")
        comm.Barrier()
        pc.done()
    except BaseException:
        # A single-rank error must not leave peers blocked in collective calls.
        import traceback
        traceback.print_exc()
        comm.Abort(1)
        raise
