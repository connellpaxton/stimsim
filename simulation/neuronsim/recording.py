"""Owned native voltage/spike recorders and optional branch geometry export."""

import csv
import json

import numpy as np


class PopulationRecording:
    """Record every step on the owning cell thread, then downsample by stride."""

    def __init__(self, h, pc, cells):
        self.ordered = sorted(cells)
        self.owners = []
        self.vectors = []
        for root in self.ordered:
            section = cells[root].soma[0]
            anchor = h.IClamp(section(0.5))
            anchor.amp = anchor.dur = 0.0
            self.owners.append(anchor)
            self.vectors.append(h.Vector().record(anchor, section(0.5)._ref_v))
        self.spike_times = h.Vector()
        self.spike_gids = h.Vector()
        pc.spike_record(-1, self.spike_times, self.spike_gids)
        self.gid_to_root = {cell.gid: root for root, cell in cells.items()}

    def clear_spikes(self):
        self.spike_times.resize(0)
        self.spike_gids.resize(0)

    def read(self, dt_ms, stride):
        voltages = np.array([np.asarray(vector)[::stride] for vector in self.vectors])
        times = np.arange(voltages.shape[1]) * (dt_ms * stride)
        event_times = np.asarray(self.spike_times).copy()
        gids = np.asarray(self.spike_gids, dtype=np.int64).copy()
        # Canonical order makes cross-thread comparisons meaningful.
        order = np.lexsort((gids, event_times))
        event_times, gids = event_times[order], gids[order]
        roots = np.asarray([self.gid_to_root[int(gid)] for gid in gids], dtype="U18")
        return times, voltages, event_times, roots

    def save(self, destination, times, voltages, event_times, roots):
        np.savez_compressed(
            destination / "recordings.npz",
            time_ms=times,
            soma_mv=voltages,
            root_ids=np.asarray(self.ordered),
            spike_times_ms=event_times,
            spike_root_ids=roots,
        )


class BranchRecording:
    """Record individual cable segments and their actual imported 3D bends."""

    def __init__(self, cells, roots, section_names, output):
        from neuron import h

        self.owners = []
        self.vectors = []
        rows = []
        for root in roots:
            cell = cells[root]
            owner = h.IClamp(cell.soma[0](0.5))
            owner.amp = owner.dur = 0.0
            self.owners.append(owner)
            for name, section in section_names[root].items():
                xyz = np.array([
                    [section.x3d(i), section.y3d(i), section.z3d(i)]
                    for i in range(section.n3d())
                ])
                arcs = np.array([section.arc3d(i) for i in range(section.n3d())])
                if len(arcs) < 2 or arcs[-1] <= 0:
                    raise ValueError("Invalid recorded branch geometry")
                for index, segment in enumerate(section):
                    low = index / section.nseg * arcs[-1]
                    high = (index + 1) / section.nseg * arcs[-1]
                    samples = np.r_[low, arcs[(arcs > low) & (arcs < high)], high]
                    points = np.column_stack([
                        np.interp(samples, arcs, xyz[:, axis]) for axis in range(3)
                    ])
                    rows.append({
                        "root_id": root, "section": name, "x": float(segment.x),
                        "points_um": points.tolist(),
                    })
                    # Every section of this cell is on the owner's native thread.
                    self.vectors.append(h.Vector().record(owner, segment._ref_v))
        self.geometry = {
            "roots": roots, "units": "um",
            "recording": "Actual per-segment membrane voltage; imported 3D cable polylines",
            "rows": rows,
        }
        output.write_text(json.dumps(self.geometry) + "\n")

    def save(self, destination, times, stride=4):
        values = np.array([
            np.asarray(vector)[::stride] for vector in self.vectors
        ], dtype=np.float32)
        if values.shape != (len(self.vectors), len(times)) or not np.isfinite(values).all():
            raise RuntimeError("Invalid branch voltage recording")
        np.savez_compressed(
            destination / "branch_recordings.npz", time_ms=times, segment_mv=values
        )
        return {
            "segments": len(self.vectors), "samples": int(values.size),
            "storage_dtype": "float32",
            "interpretation": "Visualization copies of native voltage recordings, sampled every four fixed steps",
        }


def write_cell_assignments(output, cells, lookup, ordered, stimulated):
    with (output / "cell_models.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow([
            "root_id", "gid", "manual_type", "ei_class", "electrical_profile", "stimulated"
        ])
        for root in ordered:
            writer.writerow([
                root, cells[root].gid, lookup[root]["manual_cell_type"],
                lookup[root]["ei_class"], "cACint209", root in stimulated,
            ])
