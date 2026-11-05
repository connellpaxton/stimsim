"""Index downloaded Minnie65 cell identities and original SWC provenance."""

import argparse
from collections import Counter
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import sqlite3
import struct


def decode_point(value):
    """Decode the exported 3D WKB point without a GIS dependency."""
    data = bytes.fromhex(value)
    byte_order = "<" if data[0] == 1 else ">"
    kind = struct.unpack_from(byte_order + "I", data, 1)[0]
    if kind != 0x80000001 or len(data) != 29:
        raise ValueError("Expected an EWKB 3D point without an SRID")
    point = struct.unpack_from(byte_order + "ddd", data, 5)
    if not all(math.isfinite(v) for v in point):
        raise ValueError("Non-finite cell position")
    return point


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=Path("/mnt/nas-storage/NeuronSim/datasets/minnie65/v1300"))
    args = parser.parse_args()
    base = args.data.resolve()
    metadata = base / "metadata"
    catalog = [json.loads(line) for line in (base / "SOURCE_CATALOG.jsonl").open()]
    skeletons = {row["root_id"]: row for row in catalog if row["local_path"].endswith(".swc")}
    columns = [row[0] for row in csv.reader((metadata / "aibs_cell_info_header.csv").open())]
    source = metadata / "aibs_cell_info.csv.gz"
    resolution = json.loads((metadata / "nucleus_detection_v0_table_metadata.json").read_text())["voxel_resolution"]
    with gzip.open(source, "rt") as stream:
        original = list(csv.DictReader(stream, fieldnames=columns))
    multiplicities = Counter(row["pt_root_id"] for row in original)
    fields = ["nucleus_id", "root_id", "supervoxel_id", "x_um", "y_um", "z_um", "broad_type", "broad_type_source", "cell_type", "cell_type_source", "mtype", "mtype_source", "dendrite_cleaned", "axon_cleaned", "dendrite_strategy", "axon_strategy", "visual_area", "root_nucleus_multiplicity", "swc_available_in_catalog", "swc_source_kind", "swc_path"]
    rows, errors = [], []
    for row in original:
        try:
            point = decode_point(row["pt_position"])
        except Exception as exc:
            errors.append({"nucleus_id": row["id"], "error": str(exc)})
            point = None
        root = row["pt_root_id"]
        skeleton = skeletons.get(root)
        result = {key: row.get(key, "") for key in fields}
        result.update(nucleus_id=row["id"], root_id=root, supervoxel_id=row["pt_supervoxel_id"], root_nucleus_multiplicity=multiplicities[root], swc_available_in_catalog=int(skeleton is not None), swc_source_kind=skeleton["source_kind"] if skeleton else "", swc_path=skeleton["local_path"] if skeleton else "")
        for axis, coordinate, voxel_size in zip("xyz", point or (None, None, None), resolution):
            result[axis + "_um"] = coordinate * voxel_size / 1000 if point else None
        rows.append(result)
    temporary = metadata / "cells.csv.gz.tmp"
    with gzip.open(temporary, "wt", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(metadata / "cells.csv.gz")
    temporary_db = metadata / "cells.sqlite.tmp"
    if temporary_db.exists():
        raise RuntimeError("An unfinished index exists; preserve it for inspection before retrying")
    connection = sqlite3.connect(temporary_db)
    numeric = {"x_um", "y_um", "z_um"}
    integer = {"root_nucleus_multiplicity", "swc_available_in_catalog"}
    definitions = [f'"{field}" {"REAL" if field in numeric else "INTEGER" if field in integer else "TEXT"}' for field in fields]
    connection.execute("CREATE TABLE cells (" + ",".join(definitions) + ")")
    connection.executemany("INSERT INTO cells VALUES (" + ",".join("?" for _ in fields) + ")", [[row[field] for field in fields] for row in rows])
    for field in ["nucleus_id", "root_id", "broad_type", "swc_available_in_catalog"]:
        connection.execute(f'CREATE INDEX cells_{field} ON cells ("{field}")')
    connection.commit()
    assert connection.execute("PRAGMA integrity_check").fetchone()[0] == "ok"
    connection.close()
    temporary_db.replace(metadata / "cells.sqlite")

    differences = []
    for row in rows:
        path = base / row["swc_path"] if row["swc_path"] else None
        if not path or not path.is_file() or row["x_um"] is None or row["root_nucleus_multiplicity"] != 1:
            continue
        with path.open() as stream:
            for line in stream:
                if line.strip() and not line.lstrip().startswith("#"):
                    values = line.split()
                    if len(values) == 7 and values[6] == "-1":
                        distance = math.sqrt(sum((float(values[index + 2]) - row[axis + "_um"]) ** 2 for index, axis in enumerate("xyz")))
                        differences.append(distance)
                    break
        if len(differences) >= 1000:
            break
    report = {"scope": "Derived cell metadata index; original compressed table and SWCs remain unchanged", "materialization_version": 1300, "cells": len(rows), "distinct_root_ids": len(multiplicities), "metadata_broad_types": dict(Counter(row["broad_type"] for row in rows)), "available_swc_root_ids": len(skeletons), "swc_roots_with_cell_metadata": len(set(skeletons) & set(multiplicities)), "cell_rows_with_swc": sum(row["swc_available_in_catalog"] for row in rows), "coordinate_resolution_nm": resolution, "coordinate_conversion": "EWKB point coordinate * voxel_resolution_nm / 1000 -> micrometers", "coordinate_parse_errors": errors, "sample_swc_soma_matches": {"checked": len(differences), "within_0_01_um": sum(x <= .01 for x in differences), "maximum_difference_um": max(differences) if differences else None}, "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(), "derived_files": {name: hashlib.sha256((metadata / name).read_bytes()).hexdigest() for name in ["cells.csv.gz", "cells.sqlite"]}, "notes": ["All root, nucleus and supervoxel IDs are stored as TEXT/string, never floating point.", "SWC availability describes catalog membership; consult DOWNLOAD_STATUS.json for transfer completion.", "Reconstruction completeness and nonpositive radii have not been repaired or hidden."]}
    (metadata / "index_validation.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
