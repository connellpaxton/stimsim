"""Check destination ownership and unchanged global contact offsets."""

import csv
import json
from pathlib import Path
import tempfile
import unittest

from neuronsim.sharding import prepare


class MPIShardingTests(unittest.TestCase):
    def test_each_contact_is_owned_once_and_retains_its_global_random_offset(self):
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            data, morphology = base / "data", base / "morphologies"
            data.mkdir()
            (morphology / "swc").mkdir(parents=True)
            nodes = [
                {"root_id": str(i), "nucleus_id": str(i + 10), "ei_class": "I" if i == 2 else "E", "manual_cell_type": "BC" if i == 2 else "5P-IT"}
                for i in range(3)
            ]
            edges = [
                {"pre_pt_root_id": "0", "post_pt_root_id": "1", "n_syn": "3"},
                {"pre_pt_root_id": "1", "post_pt_root_id": "0", "n_syn": "2"},
                {"pre_pt_root_id": "2", "post_pt_root_id": "1", "n_syn": "1"},
            ]
            for name, rows in (("nodes.csv", nodes), ("edges.csv", edges)):
                with (data / name).open("w", newline="") as stream:
                    writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
                    writer.writeheader()
                    writer.writerows(rows)
            (data / "summary.json").write_text(json.dumps({"internal_anatomical_synapses": 6}))
            records = [{"root_id": str(i), "nucleus_id": str(i + 10), "nodes": count} for i, count in enumerate([7, 2, 1])]
            (morphology / "manifest.json").write_text(json.dumps({"records": records}))
            result = prepare(data, morphology, base / "shards", 2)
            self.assertEqual(result["contacts"], 6)
            self.assertEqual(result["contacts_per_rank"], [2, 4])
            self.assertEqual(result["stimulus"]["cells"], ["0", "1"])
            local_edges, local_nodes = [], []
            for rank in range(2):
                directory = base / "shards" / f"rank{rank:04d}"
                with (directory / "data/edges.csv").open() as stream:
                    rows = list(csv.DictReader(stream))
                with (directory / "data/nodes.csv").open() as stream:
                    roots = {row["root_id"] for row in csv.DictReader(stream)}
                self.assertTrue(all(row["post_pt_root_id"] in roots for row in rows))
                self.assertEqual((directory / "morphologies/swc").resolve(), (morphology / "swc").resolve())
                local_edges.extend(rows)
                local_nodes.extend(roots)
            self.assertEqual(sorted(local_nodes), ["0", "1", "2"])
            ordered = sorted(local_edges, key=lambda row: int(row["global_contact_offset"]))
            self.assertEqual([int(row.pop("global_contact_offset")) for row in ordered], [0, 3, 5])
            self.assertEqual(ordered, edges)


if __name__ == "__main__":
    unittest.main()
