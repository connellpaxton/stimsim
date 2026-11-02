"""Download the public v1300 Minnie65 simulation inputs with source checksums.

Only SWCs, connectivity CSVs, cell metadata and schemas are selected. Objects are
pinned by GCS generation; final names appear only after size/MD5 verification.
"""

import argparse
import base64
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED, as_completed
import datetime
import fcntl
import hashlib
import http.client
import json
import os
from pathlib import Path
import shutil
import socket
import ssl
import threading
import time
import urllib.parse
import urllib.request

DOCUMENTATION = "https://tutorial.microns-explorer.org/static-repositories.html"
DEFAULT_OUTPUT = Path("/mnt/nas-storage/NeuronSim/datasets/minnie65/v1300")
LOCAL = threading.local()
TLS = ssl.create_default_context()


def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def atomic_json(path, value):
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n")
    temporary.replace(path)


def list_objects(bucket, prefix):
    """List only object metadata; no image/mesh buckets are enumerated."""
    items, token = [], None
    while True:
        query = {"prefix": prefix, "maxResults": 1000, "fields": "nextPageToken,items(name,size,md5Hash,crc32c,generation,updated,contentEncoding)"}
        if token:
            query["pageToken"] = token
        url = f"https://storage.googleapis.com/storage/v1/b/{bucket}/o?" + urllib.parse.urlencode(query)
        for attempt in range(5):
            try:
                with urllib.request.urlopen(url, timeout=60) as source:
                    result = json.load(source)
                break
            except Exception:
                if attempt == 4:
                    raise
                time.sleep(min(2 ** attempt, 15))
        items.extend(result.get("items", []))
        token = result.get("nextPageToken")
        if not token:
            return items


def make_catalog(output):
    path = output / "SOURCE_CATALOG.jsonl"
    if path.exists():
        return [json.loads(line) for line in path.open() if line.strip()]
    with ThreadPoolExecutor(max_workers=3) as executor:
        futures = {kind: executor.submit(list_objects, "microns-static-links", f"skel/swc/{kind}/") for kind in ("proofread", "dendrite")}
        table_future = executor.submit(list_objects, "mat_dbs", "public/minnie65_phase3_v1/v1300/")
        raw = {kind: future.result() for kind, future in futures.items()}
        tables = table_future.result()
    selected = {}
    for kind in ("dendrite", "proofread"):
        for item in raw[kind]:
            if item["name"].endswith(".swc"):
                root = Path(item["name"]).stem
                if not root.isdigit():
                    raise ValueError("Unexpected skeleton object name")
                selected[root] = dict(item, bucket="microns-static-links", source_kind=kind, local_path=f"swc/{root}.swc", root_id=root)
    records = list(selected.values())
    wanted = {"aibs_cell_info.csv.gz", "aibs_cell_info_header.csv", "connections_with_nuclei.csv.gz", "connections_with_nuclei_header.csv", "synapses_pni_2_v1_filtered_view.csv.gz", "synapses_pni_2_v1_filtered_view_header.csv"}
    found = set()
    for item in tables:
        name = Path(item["name"]).name
        if name in wanted:
            found.add(name)
            directory = "metadata" if name.startswith("aibs_") else "connectivity"
            records.append(dict(item, bucket="mat_dbs", source_kind="table", local_path=f"{directory}/{name}"))
    if found != wanted:
        raise RuntimeError(f"Missing table objects: {wanted - found}")
    for record in records:
        if not record.get("md5Hash") or not record.get("generation"):
            raise RuntimeError(f"Object lacks an immutable generation or MD5: {record['name']}")
        record["size"] = int(record["size"])
        record["url"] = "https://storage.googleapis.com/" + record["bucket"] + "/" + urllib.parse.quote(record["name"], safe="/") + "?generation=" + record["generation"]
    records.sort(key=lambda row: row["local_path"])
    temporary = path.with_suffix(".jsonl.tmp")
    with temporary.open("w") as stream:
        for record in records:
            stream.write(json.dumps(record, separators=(",", ":")) + "\n")
    temporary.replace(path)
    atomic_json(output / "INVENTORY.json", {"created_utc": now(), "materialization_version": 1300, "source_documentation": DOCUMENTATION, "files": len(records), "bytes": sum(r["size"] for r in records), "unique_swc_root_ids": len(selected), "swc_source_kind_counts": {kind: sum(r.get("source_kind") == kind for r in records) for kind in ("proofread", "dendrite")}, "selection": "One SWC per root ID, preferring the proofread file when both exports exist. Original bytes are retained.", "catalog_counts_before_deduplication": {kind: sum(i['name'].endswith('.swc') for i in items) for kind, items in raw.items()}, "excluded_categories": ["EM imagery", "segmentation voxels", "surface meshes", "meshwork HDF5", "functional imaging movies"]})
    return records


def hashes(path):
    md5, sha = hashlib.md5(), hashlib.sha256()
    with path.open("rb") as stream:
        while chunk := stream.read(4 * 1024 * 1024):
            md5.update(chunk)
            sha.update(chunk)
    return md5, sha


def connection():
    if not hasattr(LOCAL, "connection"):
        LOCAL.connection = http.client.HTTPSConnection("storage.googleapis.com", timeout=90, context=TLS)
    return LOCAL.connection


def reset_connection():
    if hasattr(LOCAL, "connection"):
        LOCAL.connection.close()
        del LOCAL.connection


class Progress:
    def __init__(self):
        self.lock = threading.Lock()
        self.transferred = 0
        self.active = {}

    def add(self, count, local_path, position):
        with self.lock:
            self.transferred += count
            self.active[local_path] = position if position is not None else self.active.get(local_path, 0) + count


def download_ranged(record, output, progress):
    """Write independent HTTP ranges directly into one resumable NAS file."""
    path = output / record["local_path"]
    partial = path.with_name(path.name + ".part")
    state_path = path.with_name(path.name + ".ranges.json")
    done_path = path.with_name(path.name + ".ranges-done.jsonl")
    if state_path.exists():
        state = json.loads(state_path.read_text())
        assert state["generation"] == record["generation"] and state["size"] == record["size"]
    else:
        prefix = partial.stat().st_size if partial.exists() else 0
        if prefix > record["size"]:
            raise RuntimeError("Existing partial exceeds the object size")
        block = 64 * 1024 ** 2
        ranges = [[start, min(start + block, record["size"]) - 1] for start in range(prefix, record["size"], block)]
        state = {"generation": record["generation"], "size": record["size"], "legacy_contiguous_prefix_bytes": prefix, "ranges": ranges}
        atomic_json(state_path, state)
    if not partial.exists():
        partial.touch()
    with partial.open("r+b") as stream:
        stream.truncate(record["size"])
    completed = set()
    if done_path.exists():
        for line in done_path.open():
            try:
                entry = json.loads(line)
                if entry["generation"] == record["generation"]:
                    completed.add((entry["start"], entry["end"]))
            except json.JSONDecodeError:
                # An interrupted final journal line is not evidence of a range.
                continue
    todo = [(start, end) for start, end in state["ranges"] if (start, end) not in completed]
    with progress.lock:
        progress.active[record["local_path"]] = state["legacy_contiguous_prefix_bytes"] + sum(end - start + 1 for start, end in completed)

    def transfer(bounds):
        start, end = bounds
        offset = start
        sha = hashlib.sha256()
        for attempt in range(6):
            try:
                parsed = urllib.parse.urlsplit(record["url"])
                conn = connection()
                conn.request("GET", parsed.path + "?" + parsed.query, headers={"Range": f"bytes={offset}-{end}", "Accept-Encoding": "gzip"})
                response = conn.getresponse()
                expected_range = f"bytes {offset}-{end}/{record['size']}"
                if response.status != 206 or response.getheader("Content-Range") != expected_range:
                    raise RuntimeError(f"Unexpected range response: {response.status} {response.getheader('Content-Range')}")
                fd = os.open(partial, os.O_WRONLY)
                try:
                    while chunk := response.read(1024 * 1024):
                        if offset + len(chunk) > end + 1:
                            raise RuntimeError("Range body exceeds expected boundary")
                        written = 0
                        while written < len(chunk):
                            count = os.pwrite(fd, chunk[written:], offset + written)
                            if count <= 0:
                                raise RuntimeError("Range write made no progress")
                            written += count
                        sha.update(chunk)
                        offset += len(chunk)
                        progress.add(len(chunk), record["local_path"], None)
                finally:
                    os.close(fd)
                if offset != end + 1:
                    raise RuntimeError("Truncated range response")
                return {"start": start, "end": end, "generation": record["generation"], "sha256": sha.hexdigest(), "completed_utc": now()}
            except Exception:
                reset_connection()
                if attempt == 5:
                    raise
                time.sleep(min(2 ** attempt, 30))

    with done_path.open("a", buffering=1) as journal, ThreadPoolExecutor(max_workers=32) as executor:
        futures = [executor.submit(transfer, bounds) for bounds in todo]
        for future in as_completed(futures):
            journal.write(json.dumps(future.result(), separators=(",", ":")) + "\n")
    md5, sha = hashes(partial)
    if partial.stat().st_size != record["size"] or md5.digest() != base64.b64decode(record["md5Hash"]):
        raise RuntimeError("Assembled ranged download failed whole-object size/MD5 verification")
    partial.replace(path)
    with progress.lock:
        progress.active.pop(record["local_path"], None)
    return dict(local_path=record["local_path"], size=record["size"], sha256=sha.hexdigest(), md5_base64=record["md5Hash"], source_generation=record["generation"], verified_utc=now(), resumed_existing=False, download_method="32 parallel immutable HTTP ranges; complete source MD5 verified")


def download(record, output, progress):
    path = output / record["local_path"]
    path.parent.mkdir(parents=True, exist_ok=True)
    expected = base64.b64decode(record["md5Hash"])
    if path.exists():
        if path.stat().st_size != record["size"]:
            raise RuntimeError(f"Existing final file has unexpected size: {path}")
        md5, sha = hashes(path)
        if md5.digest() != expected:
            raise RuntimeError(f"Existing final file checksum mismatch: {path}")
        return dict(local_path=record["local_path"], size=record["size"], sha256=sha.hexdigest(), md5_base64=record["md5Hash"], source_generation=record["generation"], verified_utc=now(), resumed_existing=True)
    if record["size"] >= 128 * 1024 ** 2:
        return download_ranged(record, output, progress)
    partial = path.with_name(path.name + ".part")
    error = None
    for attempt in range(6):
        try:
            position = partial.stat().st_size if partial.exists() else 0
            if position > record["size"]:
                raise RuntimeError("Partial file exceeds expected object size")
            md5, sha = hashes(partial) if position else (hashlib.md5(), hashlib.sha256())
            if position < record["size"]:
                headers = {"Accept-Encoding": "gzip"}
                if position:
                    headers["Range"] = f"bytes={position}-"
                parsed = urllib.parse.urlsplit(record["url"])
                conn = connection()
                conn.request("GET", parsed.path + "?" + parsed.query, headers=headers)
                response = conn.getresponse()
                expected_status = 206 if position else 200
                if response.status != expected_status:
                    status = response.status
                    response.read(4096)
                    reset_connection()
                    raise RuntimeError(f"HTTP {status}, expected {expected_status}")
                if position and not response.getheader("Content-Range", "").startswith(f"bytes {position}-"):
                    raise RuntimeError("Server returned an unexpected range")
                with partial.open("ab" if position else "wb") as target:
                    while chunk := response.read(1024 * 1024):
                        target.write(chunk)
                        md5.update(chunk)
                        sha.update(chunk)
                        position += len(chunk)
                        progress.add(len(chunk), record["local_path"], position)
                if position != record["size"]:
                    raise RuntimeError(f"Incomplete object: {position}/{record['size']}")
            if md5.digest() != expected:
                # This file belongs to this transfer; retain the failed bytes for
                # inspection while starting a clean retry of the same generation.
                partial.replace(partial.with_name(partial.name + f".bad-{time.time_ns()}"))
                raise RuntimeError("Source MD5 verification failed")
            partial.replace(path)
            with progress.lock:
                progress.active.pop(record["local_path"], None)
            return dict(local_path=record["local_path"], size=record["size"], sha256=sha.hexdigest(), md5_base64=record["md5Hash"], source_generation=record["generation"], verified_utc=now(), resumed_existing=False)
        except Exception as exc:
            error = f"{type(exc).__name__}: {exc}"
            reset_connection()
            if attempt < 5:
                time.sleep(min(2 ** attempt, 30))
    raise RuntimeError(f"{record['local_path']}: {error}")


def write_readme(output, summary):
    inventory = json.loads((output / "INVENTORY.json").read_text())
    text = f"""# Minnie65 v1300 simulation data

Status: **{summary['status']}**. Verified files: {summary['verified_files']:,}/{summary['expected_files']:,}.

This collection contains {inventory['unique_swc_root_ids']:,} unique SWC root IDs, connectivity and cell metadata from the publicly documented v1300 static exports. It was downloaded for NeuronSim; no anatomical data was resampled or repaired.

## Files

- `swc/<root_id>.swc`: one original cell-shape file per root ID. The proofread export is preferred where both catalogs contain that ID. Most remaining exports are dendrite-only.
- `connectivity/connections_with_nuclei.csv.gz`: aggregated presynaptic/postsynaptic cell-pair counts and nucleus IDs.
- `connectivity/synapses_pni_2_v1_filtered_view.csv.gz`: individual synapse IDs, pre/post/center positions and endpoint IDs.
- `connectivity/*_header.csv`: column names and SQL types. Data CSVs themselves do not contain column headers.
- `metadata/aibs_cell_info.csv.gz` and its header: stable nucleus identity, root ID, cell classification, visual area and axon/dendrite proofreading information.
- `SOURCE_CATALOG.jsonl`: original URLs, immutable GCS generations, MD5/CRC32C checksums, source kind and local paths.
- `VERIFIED_FILES.jsonl` and `SHA256SUMS`: verified local checksums.
- `DOWNLOAD_STATUS.json`: progress and completion summary.

Root and supervoxel IDs must stay exact integers or strings, never floating-point numbers. SWC geometry uses micrometers. Preserve the original table coordinate conventions; consult the accompanying CAVE metadata before converting coordinates. These data do not include fitted electrical channels, synaptic strengths or delays. NeuronSim's published mechanisms and electrical assumptions remain separate.

File integrity verification establishes that the downloads match their sources. Individual reconstructions still require morphology and proofreading checks before simulation. Raw source bytes remain intact, including any zero radii or incomplete branches.

Source documentation: {DOCUMENTATION}

MICrONS attribution and terms: https://www.microns-explorer.org/terms-and-conditions

Reproduce/resume from the NeuronSim repository:

```sh
.venv/bin/python tools/download_minnie_simulation_data.py --output {output}
```
"""
    (output / "README.md").write_text(text)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--workers", type=int, default=32)
    args = parser.parse_args()
    if not 1 <= args.workers <= 96:
        parser.error("Choose 1–96 download workers")
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    lock_file = (output / ".download.lock").open("a+")
    try:
        fcntl.flock(lock_file, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        raise SystemExit("Another downloader owns this destination")
    lock_file.seek(0)
    lock_file.truncate()
    lock_file.write(f"{socket.gethostname()} pid={os.getpid()} started={now()}\n")
    lock_file.flush()
    print("Inventorying simulation-only public objects", flush=True)
    records = make_catalog(output)
    needed = sum(r["size"] for r in records)
    if shutil.disk_usage(output).free < needed + 2 * 1024 ** 3:
        raise SystemExit("Insufficient free disk space for the selected inputs")
    started = time.monotonic()
    progress = Progress()
    verified, failed = {}, []
    status = {"status": "downloading", "started_utc": now(), "output": str(output), "expected_files": len(records), "expected_bytes": needed, "verified_files": 0, "verified_bytes": 0, "workers": args.workers}
    atomic_json(output / "DOWNLOAD_STATUS.json", status)
    write_readme(output, status)
    print(json.dumps(status), flush=True)
    # Start both large connection tables immediately; other workers fetch SWCs.
    ordered = sorted(records, key=lambda r: (r["source_kind"] != "table", -r["size"]))
    remaining = iter(ordered)
    last_report = 0.
    with (output / "VERIFIED_FILES.jsonl").open("a", buffering=1) as audit, ThreadPoolExecutor(max_workers=args.workers) as executor:
        pending = {}
        for _ in range(args.workers * 2):
            record = next(remaining, None)
            if record:
                pending[executor.submit(download, record, output, progress)] = record
        while pending:
            done, _ = wait(pending, timeout=2, return_when=FIRST_COMPLETED)
            for future in done:
                record = pending.pop(future)
                try:
                    result = future.result()
                    verified[result["local_path"]] = result
                    audit.write(json.dumps(result, separators=(",", ":")) + "\n")
                except Exception as exc:
                    failure = {"local_path": record["local_path"], "error": str(exc)}
                    failed.append(failure)
                    print("FAILED", json.dumps(failure), flush=True)
                record = next(remaining, None)
                if record:
                    pending[executor.submit(download, record, output, progress)] = record
            elapsed = time.monotonic() - started
            if elapsed - last_report >= 10 or not pending:
                with progress.lock:
                    transferred = progress.transferred
                    active = dict(progress.active)
                status.update(updated_utc=now(), elapsed_seconds=elapsed, verified_files=len(verified), verified_bytes=sum(r["size"] for r in verified.values()), transferred_bytes_this_run=transferred, transfer_MB_per_second=transferred / max(elapsed, .01) / 1e6, failed_files=len(failed), active_download_bytes=active)
                atomic_json(output / "DOWNLOAD_STATUS.json", status)
                print(json.dumps({k: status[k] for k in ("elapsed_seconds", "verified_files", "expected_files", "verified_bytes", "transferred_bytes_this_run", "transfer_MB_per_second", "failed_files")}), flush=True)
                last_report = elapsed
    with (output / "SHA256SUMS").open("w") as stream:
        for name, result in sorted(verified.items()):
            stream.write(result["sha256"] + "  " + name + "\n")
    status.update(status="complete" if len(verified) == len(records) and not failed else "incomplete", completed_utc=now(), failures=failed, verified_files=len(verified), verified_bytes=sum(r["size"] for r in verified.values()), integrity="Every final file matched its immutable source generation, source MD5 and byte count; local SHA256 retained")
    atomic_json(output / "DOWNLOAD_STATUS.json", status)
    write_readme(output, status)
    print("FINAL", json.dumps(status), flush=True)
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
