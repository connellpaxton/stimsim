"""Protect resumable downloads and preservation of existing data."""

import base64
import hashlib
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
import download_minnie_simulation_data as downloader


class FakeConnection:
    def __init__(self, content):
        self.content = content
        self.calls = []

    def request(self, method, url, headers):
        self.calls.append((method, url, headers))
        bounds = headers.get("Range")
        if bounds:
            start, end = bounds.removeprefix("bytes=").split("-")
            self.start = int(start)
            self.end = int(end) if end else len(self.content) - 1
            self.status = 206
        else:
            self.start, self.end = 0, len(self.content) - 1
            self.status = 200
        self.cursor = self.start

    def getresponse(self):
        return self

    def getheader(self, key, default=None):
        if key == "Content-Range":
            return f"bytes {self.start}-{self.end}/{len(self.content)}"
        return default

    def read(self, size):
        end = min(self.cursor + size, self.end + 1)
        block = self.content[self.cursor:end]
        self.cursor = end
        return block

    def close(self):
        pass


class DownloadTests(unittest.TestCase):
    def record(self, content):
        return {"local_path": "swc/864691135539486066.swc", "size": len(content), "md5Hash": base64.b64encode(hashlib.md5(content).digest()).decode(), "generation": "123456789", "url": "https://storage.googleapis.com/example/object?generation=123456789"}

    def test_partial_transfer_resumes_without_duplicating_prefix(self):
        content = b"original source bytes\n" * 1000
        record = self.record(content)
        connection = FakeConnection(content)
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            partial = output / (record["local_path"] + ".part")
            partial.parent.mkdir()
            partial.write_bytes(content[:137])
            with patch.object(downloader, "connection", return_value=connection):
                result = downloader.download(record, output, downloader.Progress())
            self.assertEqual((output / record["local_path"]).read_bytes(), content)
            self.assertEqual(connection.calls[0][2]["Range"], "bytes=137-")
            self.assertIn("generation=123456789", connection.calls[0][1])
            self.assertEqual(result["sha256"], hashlib.sha256(content).hexdigest())
            self.assertFalse(partial.exists())

    def test_ranged_file_keeps_legacy_prefix_and_verifies_whole_object(self):
        content = b"neuron-data" * 200000
        record = self.record(content)
        connection = FakeConnection(content)
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            partial = output / (record["local_path"] + ".part")
            partial.parent.mkdir()
            partial.write_bytes(content[:917])
            with patch.object(downloader, "connection", return_value=connection):
                result = downloader.download_ranged(record, output, downloader.Progress())
            self.assertEqual((output / record["local_path"]).read_bytes(), content)
            self.assertEqual(connection.calls[0][2]["Range"], f"bytes=917-{len(content)-1}")
            self.assertEqual(result["sha256"], hashlib.sha256(content).hexdigest())
            self.assertTrue((output / (record["local_path"] + ".ranges-done.jsonl")).is_file())

    def test_existing_wrong_file_is_preserved_and_never_overwritten(self):
        record = self.record(b"expected bytes")
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            path = output / record["local_path"]
            path.parent.mkdir()
            original = b"different data"
            self.assertEqual(len(original), record["size"])
            path.write_bytes(original)
            with patch.object(downloader, "connection") as connection:
                with self.assertRaisesRegex(RuntimeError, "checksum mismatch"):
                    downloader.download(record, output, downloader.Progress())
            connection.assert_not_called()
            self.assertEqual(path.read_bytes(), original)


if __name__ == "__main__":
    unittest.main()
