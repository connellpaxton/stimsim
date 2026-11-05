"""Check source integrity and deterministic parallel morphology preparation."""

import base64
from concurrent.futures import ProcessPoolExecutor
import hashlib
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
from prepare_nas_network import prepare_candidate


class NASPreparationTests(unittest.TestCase):
    def make_task(self, directory, name):
        source = directory / "source.swc"
        source.write_text("0 1 0 0 0 5 -1\n1 3 10 0 0 1 0\n2 3 10 0 0 1 1\n3 3 20 0 0 0 2\n")
        node = {"root_id": "864691135999999999", "nucleus_id": "123", "x_um": 0, "y_um": 0, "z_um": 0}
        record = {"md5Hash": base64.b64encode(hashlib.md5(source.read_bytes()).digest()).decode(), "url": "https://example.test/source.swc", "generation": "1", "source_kind": "dendrite"}
        return node, source, directory / name, record

    def test_parallel_preparation_preserves_geometry_and_source_bytes(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            serial_task = self.make_task(directory, "serial.swc")
            original = serial_task[1].read_bytes()
            serial = prepare_candidate(serial_task)
            parallel_task = (serial_task[0], serial_task[1], directory / "parallel.swc", serial_task[3])
            with ProcessPoolExecutor(max_workers=2) as pool:
                parallel = list(pool.map(prepare_candidate, [parallel_task]))[0]
            self.assertEqual(serial, parallel)
            self.assertIsNone(serial[2])
            self.assertEqual(serial[1]["nodes"], 3)
            self.assertEqual(serial[1]["cable_length_um"], 20)
            self.assertEqual(len(serial[1]["repairs"]), 1)
            self.assertEqual(len(serial[1]["collapsed_identical_vertices"]), 1)
            self.assertEqual(serial_task[1].read_bytes(), original)
            self.assertEqual(serial_task[2].read_bytes(), parallel_task[2].read_bytes())

    def test_changed_source_is_excluded_before_an_output_is_written(self):
        with tempfile.TemporaryDirectory() as directory:
            task = self.make_task(Path(directory), "output.swc")
            task[1].write_text(task[1].read_text() + "# changed source\n")
            _, record, error = prepare_candidate(task)
            self.assertIsNone(record)
            self.assertIn("published MD5", error["reason"])
            self.assertFalse(task[2].exists())

    def test_dendrite_root_is_normalized_only_with_verified_soma_identity(self):
        with tempfile.TemporaryDirectory() as directory:
            task = self.make_task(Path(directory), "normalized.swc")
            source = task[1].read_bytes().replace(b"0 1 0 0 0", b"0 3 0 0 0", 1)
            task[1].write_bytes(source)
            task[3]["md5Hash"] = base64.b64encode(hashlib.md5(source).digest()).decode()
            _, record, error = prepare_candidate(task)
            self.assertIsNone(error)
            self.assertEqual(record["root_compartment_normalization"]["original_type"], 3)
            self.assertTrue(task[2].read_text().startswith("0 1 0 0 0 5 -1"))
            self.assertEqual(task[1].read_bytes(), source)
            bad_node = {**task[0], "x_um": 10}
            bad_task = (bad_node, task[1], Path(directory) / "bad.swc", task[3])
            self.assertIn("Soma identity mismatch", prepare_candidate(bad_task)[2]["reason"])
            self.assertFalse(bad_task[2].exists())
            bad_task = (task[0], task[1], bad_task[2], {**task[3], "source_kind": "proofread"})
            self.assertIn("not rooted at a recorded soma", prepare_candidate(bad_task)[2]["reason"])


if __name__ == "__main__":
    unittest.main()
