"""Render recorded soma activity at measured positions, without branch recordings."""

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import csv
import hashlib
import html
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
import numpy as np


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def prepare(args, cache):
    recording = args.run / "recordings.npz"
    request = {
        "recording": str(recording), "recording_bytes": recording.stat().st_size,
        "recording_mtime_ns": recording.stat().st_mtime_ns,
        "nodes_sha256": digest(args.nodes), "summary_sha256": digest(args.run / "summary.json"),
        "start_ms": args.start_ms, "duration_ms": args.duration_ms,
        "fps": args.fps, "playback_seconds": args.playback_seconds,
        "width": args.width, "height": args.height,
    }
    metadata = cache / "metadata.json"
    if metadata.exists():
        result = json.loads(metadata.read_text())
        if result["request"] != request:
            raise ValueError("Render cache belongs to different inputs or settings; use a fresh output directory")
        return result
    cache.mkdir(parents=True, exist_ok=True)
    summary = json.loads((args.run / "summary.json").read_text())
    with args.nodes.open() as stream:
        nodes = {row["root_id"]: row for row in csv.DictReader(stream)}
    with np.load(recording) as archive:
        roots = archive["root_ids"]
        times = archive["time_ms"]
        count = int(round(args.fps * args.playback_seconds))
        stop = args.start_ms + args.duration_ms
        if args.start_ms < times[0] or stop > times[-1] or count < 2:
            raise ValueError("Requested interval is outside the recording, or too few video frames")
        targets = np.linspace(args.start_ms, stop, count)
        right = np.clip(np.searchsorted(times, targets), 1, len(times) - 1)
        left = right - 1
        indices = np.where(np.abs(times[left] - targets) <= np.abs(times[right] - targets), left, right)
        error = float(np.max(np.abs(times[indices] - targets)))
        if error > np.diff(times).max() / 2 + 1e-7:
            raise RuntimeError("Video sample times are misaligned")
        index = {str(root): i for i, root in enumerate(roots)}
        if len(index) != summary["neurons"] or not set(index).issubset(nodes):
            raise ValueError("Recording and position identities differ")
        xyz = np.array([[float(nodes[str(root)][axis + "_um"]) for axis in "xyz"] for root in roots])
        ei = np.array([nodes[str(root)]["ei_class"] for root in roots])
        stimulated = np.array([index[root] for root in summary["stimulus"]["cells"]], dtype=np.int32)
        spikes = archive["spike_times_ms"]
        spike_roots = archive["spike_root_ids"]
        selected = (spikes >= args.start_ms - 2) & (spikes <= stop)
        spikes = spikes[selected]
        spike_cells = np.array([index[str(root)] for root in spike_roots[selected]], dtype=np.int32)
        print(f"Reading {summary['neurons']:,} recorded soma traces; selecting {count} frames", flush=True)
        soma = archive["soma_mv"]
        if soma.shape != (len(roots), len(times)):
            raise ValueError("Unexpected recording dimensions")
        values = soma[:, indices].T
        display = values.astype(np.float32)
        conversion_error = float(np.max(np.abs(values - display)))
        if not np.isfinite(display).all() or not np.isfinite(xyz).all():
            raise ValueError("Nonfinite video data")
        for name, value in {
            "voltages": display, "times": times[indices], "sample_indices": indices,
            "xyz": xyz, "roots": roots, "ei": ei, "stimulated": stimulated,
            "spike_times": spikes, "spike_cells": spike_cells,
        }.items():
            np.save(cache / (name + ".npy"), value)
    executable = shutil.which("ffmpeg")
    if executable is None:
        import imageio_ffmpeg
        executable = imageio_ffmpeg.get_ffmpeg_exe()
    result = {
        "request": request, "summary": summary, "frames": count, "ffmpeg": executable,
        "recording_sha256": digest(recording), "source_interval_ms": [args.start_ms, stop],
        "source_recording_interval_ms": float(np.diff(times).min()),
        "maximum_video_sample_timing_error_ms": error,
        "visualization_dtype": "float32", "maximum_visualization_conversion_error_mv": conversion_error,
        "data": "All recorded neurons; soma voltages at the nearest saved sample. Two projections of measured soma coordinates. White rings show recorded spikes in the preceding 2 ms. No axonal or branch voltage propagation is synthesized.",
    }
    metadata.write_text(json.dumps(result, indent=2) + "\n")
    return result


class Scene:
    def __init__(self, cache):
        self.cache = Path(cache)
        self.config = json.loads((self.cache / "metadata.json").read_text())
        self.arrays = {name: np.load(self.cache / (name + ".npy"), mmap_mode="r") for name in
                       ["voltages", "times", "xyz", "roots", "stimulated", "spike_times", "spike_cells"]}
        request = self.config["request"]
        self.summary = self.config["summary"]
        matplotlib.rcParams["animation.ffmpeg_path"] = self.config["ffmpeg"]
        plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11,
            "text.color": "#e7eef8", "axes.labelcolor": "#b7c8dc", "xtick.color": "#9aadc5",
            "ytick.color": "#9aadc5", "axes.edgecolor": "#394b63", "axes.facecolor": "#091522",
            "savefig.facecolor": "#091522"})
        self.fig = plt.figure(figsize=(request["width"] / 100, request["height"] / 100), facecolor="#091522")
        self.fig.text(.055, .958, "MINNIE65  |  RECORDED POPULATION ACTIVITY", fontsize=24, weight="bold")
        self.fig.text(.055, .919, f"{self.summary['neurons']:,} neurons  ·  {self.summary['anatomical_contact_count']:,} modeled contacts  ·  soma voltage and spikes", fontsize=14, color="#b7c8dc")
        self.clock = self.fig.text(.055, .868, "", fontsize=17, color="#ffd166")
        self.activity = self.fig.text(.945, .868, "", fontsize=12, ha="right", color="#b7c8dc")
        norm = Normalize(-85, 40, clip=True)
        cmap = plt.get_cmap("turbo")
        self.populations = []
        stimulated = self.arrays["stimulated"]
        for bounds, axes, title in [([.065, .16, .43, .64], (0, 1), "Front projection · X–Y"),
                                   ([.55, .46, .37, .34], (0, 2), "Top projection · X–Z")]:
            ax = self.fig.add_axes(bounds)
            xy = self.arrays["xyz"][:, axes]
            points = ax.scatter(xy[:, 0], xy[:, 1], c=self.arrays["voltages"][0], s=1.6,
                                cmap=cmap, norm=norm, linewidths=0, alpha=.8, rasterized=True)
            ax.scatter(xy[stimulated, 0], xy[stimulated, 1], s=42, facecolors="none", edgecolors="#ffd166", linewidths=.9)
            halo = ax.scatter([], [], s=10, facecolors="none", edgecolors="#f7fbff", linewidths=.4, alpha=.9)
            ax.set_aspect("equal")
            ax.invert_yaxis()
            ax.set(title=title, xlabel="X (µm)", ylabel=("Y" if axes[1] == 1 else "Z") + " (µm)")
            ax.margins(.035)
            self.populations.append((points, halo, xy))
        self.fig.legend(handles=[Line2D([], [], marker="o", linestyle="", markerfacecolor="none", markeredgecolor="#ffd166", label="Injected cells"),
                                 Line2D([], [], marker="o", linestyle="", markerfacecolor="none", markeredgecolor="white", label="Spike in preceding 2 ms")],
                        loc="lower left", bbox_to_anchor=(.065, .085), ncol=2, frameon=False, fontsize=12)
        color_axis = self.fig.add_axes([.95, .485, .009, .28])
        self.fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=color_axis, ticks=[-85, -65, -40, 0, 40], extend="both")
        color_axis.set_title("mV", fontsize=11)
        self.rate_axis = self.fig.add_axes([.55, .17, .37, .18])
        start, stop = self.config["source_interval_ms"]
        edges = np.r_[np.arange(start, stop, 5.), stop]
        counts, _ = np.histogram(self.arrays["spike_times"], edges)
        self.rate_times = (edges[:-1] + edges[1:]) / 2
        self.rate_ends = edges[1:]
        self.rates = counts / (np.diff(edges) / 1000 * self.summary["neurons"])
        self.rate_line, = self.rate_axis.plot([], [], color="#6ee7da", linewidth=1.5)
        self.cursor = self.rate_axis.axvline(start, color="white", linewidth=.8, alpha=.8)
        self.rate_axis.set(xlim=(start, stop), ylim=(0, max(1, self.rates.max() * 1.15)),
                           xlabel="Model time (ms)", ylabel="Spikes / cell / s", title="Population firing rate · recorded 5 ms bins")
        stimulus = self.summary["stimulus"]
        self.rate_axis.axvspan(stimulus["onset_ms"], stimulus["onset_ms"] + stimulus["duration_ms"], color="#ffd166", alpha=.12)
        self.rate_axis.grid(alpha=.12)
        slowdown = request["playback_seconds"] / (request["duration_ms"] / 1000)
        self.fig.text(.55, .10, f"Playback slowed {slowdown:g}×  ·  clock shows recorded model time", fontsize=11, color="#b7c8dc")
        self.fig.text(.055, .045, "Recorded soma activity at measured positions; this run has no branch-voltage recording. White rings indicate spikes, not axonal travel.", fontsize=11, color="#9aadc5")
        self.fig.text(.055, .021, "Minnie65 v1300 partial cell shapes and wiring · published ModelDB 241165 electrical profile · assumed synaptic strengths, locations and delays", fontsize=10, color="#8296b2")

    def update(self, frame):
        t = float(self.arrays["times"][frame])
        spikes = self.arrays["spike_times"]
        lo, hi = np.searchsorted(spikes, [t - 2, t], side="right")
        active = np.unique(self.arrays["spike_cells"][lo:hi])
        for points, halo, xy in self.populations:
            points.set_array(self.arrays["voltages"][frame])
            halo.set_offsets(xy[active] if len(active) else np.empty((0, 2)))
        stimulus = self.summary["stimulus"]
        on = stimulus["onset_ms"] <= t < stimulus["onset_ms"] + stimulus["duration_ms"]
        stop = self.config["source_interval_ms"][1]
        self.clock.set_text(f"Model time {t:6.1f} / {stop:g} ms    |    Stimulus {'ON' if on else 'OFF'} · {stimulus['amplitude_na']:g} nA into {len(stimulus['cells'])} cells")
        self.activity.set_text(f"{hi - lo:,} spikes in preceding 2 ms")
        past = self.rate_ends <= t
        self.rate_line.set_data(self.rate_times[past], self.rates[past])
        self.cursor.set_xdata([t, t])


def render_chunk(job):
    cache, first, last, path = job
    scene = Scene(cache)
    writer = FFMpegWriter(fps=scene.config["request"]["fps"], codec="libx264",
                         extra_args=["-crf", "20", "-preset", "fast", "-pix_fmt", "yuv420p", "-threads", "1"])
    with writer.saving(scene.fig, str(path), dpi=100):
        for frame in range(first, last):
            scene.update(frame)
            writer.grab_frame()
    plt.close(scene.fig)
    return {"first": first, "last_exclusive": last, "path": str(path)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True, help="Directory containing recordings.npz and summary.json")
    parser.add_argument("--nodes", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True, help="New MP4 path")
    parser.add_argument("--start-ms", type=float, default=0)
    parser.add_argument("--duration-ms", type=float, default=500)
    parser.add_argument("--playback-seconds", type=float, default=30)
    parser.add_argument("--fps", type=int, default=30)
    parser.add_argument("--workers", type=int, default=12)
    parser.add_argument("--width", type=int, default=1920)
    parser.add_argument("--height", type=int, default=1080)
    parser.add_argument("--preview-only", action="store_true")
    args = parser.parse_args()
    if min(args.duration_ms, args.playback_seconds, args.fps, args.workers, args.width, args.height) <= 0 or args.width % 2 or args.height % 2:
        parser.error("Positive settings and even video dimensions are required")
    args.run, args.nodes, args.output = args.run.resolve(), args.nodes.resolve(), args.output.resolve()
    if args.output.exists() and not args.preview_only:
        raise FileExistsError(args.output)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    cache = args.output.parent / "render_cache"
    started = time.perf_counter()
    metadata = prepare(args, cache)
    scene = Scene(cache)
    preview_frames = sorted(set(int(np.abs(scene.arrays["times"] - target).argmin())
                                for target in [args.start_ms, 15, 35, 100, 250, args.start_ms + args.duration_ms]))
    previews = []
    for frame in preview_frames:
        scene.update(frame)
        path = args.output.parent / f"preview_{scene.arrays['times'][frame]:07.1f}ms.png"
        if not path.exists():
            scene.fig.savefig(path, dpi=100)
        previews.append(path.name)
    plt.close(scene.fig)
    print(json.dumps({"previews": previews, "frames": metadata["frames"], "preparation_seconds": time.perf_counter() - started}), flush=True)
    if args.preview_only:
        return
    workers = min(args.workers, metadata["frames"])
    bounds = np.linspace(0, metadata["frames"], workers + 1, dtype=int)
    jobs = [(cache, int(bounds[i]), int(bounds[i + 1]), cache / f"chunk_{i:03d}.mp4") for i in range(workers)]
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(render_chunk, job) for job in jobs]
        for index, future in enumerate(as_completed(futures), 1):
            result = future.result()
            print(f"Rendered chunk {index}/{workers}: frames {result['first']}–{result['last_exclusive'] - 1}; elapsed {time.perf_counter() - started:.1f} s", flush=True)
    listing = cache / "concat.txt"
    listing.write_text("".join(f"file '{job[3].name}'\n" for job in jobs))
    subprocess.run([metadata["ffmpeg"], "-hide_banner", "-loglevel", "error", "-f", "concat", "-safe", "0", "-i", str(listing),
                    "-c", "copy", "-movflags", "+faststart", str(args.output)], check=True)
    player = args.output.with_suffix(".html")
    player.write_text(f'''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Minnie population activity</title><style>body{{background:#091522;color:#e7eef8;font:16px system-ui;margin:30px auto;max-width:1500px;padding:0 24px}}video{{width:100%;background:#091522}}p{{color:#b7c8dc;line-height:1.6}}</style>
<h1>Recorded Minnie population activity</h1><p>{metadata['summary']['neurons']:,} neurons · {metadata['summary']['anatomical_contact_count']:,} contacts · model time {args.start_ms:g}–{args.start_ms + args.duration_ms:g} ms · {args.playback_seconds:g} seconds of slowed playback</p>
<video controls playsinline preload="metadata" poster="{html.escape(previews[min(1,len(previews)-1)],quote=True)}"><source src="{html.escape(args.output.name,quote=True)}" type="video/mp4"></video>
<p>Each point is a neuron at its measured soma position. Colors show recorded soma voltage. White rings indicate spikes in the preceding 2 ms; gold rings mark injected cells. Both views show the same neurons from different projections.</p>
<p>This video shows the recorded soma activity. This run did not record voltage along branches. No axonal travel or activity was synthesized. Geometry and wiring: Minnie65 v1300; electrical parameters: the published ModelDB 241165 adaptation.</p></html>''')
    metadata.update(video=str(args.output), player=str(player), preview_files=previews,
                    encoded_playback_seconds=metadata["frames"] / args.fps, render_workers=workers,
                    render_seconds=time.perf_counter() - started, video_sha256=digest(args.output))
    args.output.with_suffix(".json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps({key: metadata[key] for key in ["video", "frames", "encoded_playback_seconds", "render_seconds", "video_sha256"]}, indent=2), flush=True)


if __name__ == "__main__":
    main()
