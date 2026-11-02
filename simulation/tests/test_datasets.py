"""Dataset entry points must not silently run another cohort or contact policy."""

from dataclasses import replace
import csv
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from neuronsim.archive import prepare_inputs
from neuronsim.cli import parse_args
from neuronsim.datasets import COLUMN, L23, THREE_NEURONS, WHOLE_MINNIE, check_prepared
from neuronsim.merge_mpi import merge


def write_csv(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


class ArchiveDatasetTests(unittest.TestCase):
    def archive(self, base, dataset):
        folder = base / ("SWC_" + dataset.layer)
        folder.mkdir()
        for root in dataset.roots:
            (folder / (root + ".swc")).write_text("0 1 0 0 0 5 -1\n1 3 10 0 0 1 0\n")
        roots = list(dataset.roots)
        rows = [{"pre_pt_root_id": root, "post_pt_root_id": roots[(i + 1) % 3],
                 "pre_pt_position_x": "1234.5", "post_pt_position_x": "5678.9"}
                for i, root in enumerate(roots)]
        write_csv(base / (dataset.layer + "_all.csv"), rows)
        return rows

    def test_three_cell_identity_contact_order_and_coordinates_are_preserved(self):
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            expected = self.archive(base, THREE_NEURONS)
            result = prepare_inputs(base, THREE_NEURONS)
            self.assertEqual(result[-2], list(THREE_NEURONS.roots))
            self.assertEqual(result[-1], expected)

    def test_missing_wrapper_cell_is_not_replaced_by_another_swc(self):
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            self.archive(base, THREE_NEURONS)
            path = base / "SWC_L5" / (THREE_NEURONS.roots[0] + ".swc")
            path.rename(path.with_name("900000000000000001.swc"))
            with self.assertRaisesRegex(ValueError, "exact three wrapper SWCs"):
                prepare_inputs(base, THREE_NEURONS)

    def test_wrong_three_contact_cycle_is_rejected_even_if_counts_match(self):
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            rows = self.archive(base, THREE_NEURONS)
            rows[0]["post_pt_root_id"] = rows[0]["pre_pt_root_id"]
            write_csv(base / "L5_all.csv", rows)
            with self.assertRaisesRegex(ValueError, "three-contact cycle"):
                prepare_inputs(base, THREE_NEURONS)

    def test_l23_subset_is_selected_after_full_parent_cohort_validation(self):
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            fixture = replace(L23, roots=THREE_NEURONS.roots, neurons=3, contacts=3, available_swcs=3)
            self.archive(base, fixture)
            fixture = replace(fixture, roots=())
            result = prepare_inputs(base, fixture, limit=2)
            self.assertEqual(len(result[-2]), 2)
            self.assertEqual(len(result[1]), 3)
            with self.assertRaisesRegex(ValueError, "archive counts differ"):
                prepare_inputs(base, replace(fixture, contacts=4), limit=2)


class PreparedDatasetTests(unittest.TestCase):
    def fixture(self, base, spec=COLUMN):
        data, morphology = base / "data", base / "morphologies"
        data.mkdir(); morphology.mkdir()
        roots = THREE_NEURONS.roots[:2]
        nodes = [{"root_id": root, "nucleus_id": str(i + 10)} for i, root in enumerate(roots)]
        write_csv(data / "nodes.csv", nodes)
        write_csv(data / "edges.csv", [{"pre_pt_root_id": roots[0], "post_pt_root_id": roots[1], "n_syn": 3}])
        (data / "summary.json").write_text(json.dumps({"selected_neurons": 2, "internal_anatomical_synapses": 3, "cohort": "all"}))
        (morphology / "manifest.json").write_text(json.dumps({"ready": 2, "records": nodes,
            "reconstruction_version": spec.reconstruction_version, "wiring_version": spec.wiring_version}))
        return replace(spec, neurons=2, contacts=3), data, morphology

    def test_wrong_release_is_rejected_before_native_construction(self):
        with tempfile.TemporaryDirectory() as directory:
            spec, data, morphology = self.fixture(Path(directory))
            record = check_prepared(spec, data, morphology)
            self.assertEqual(record["reported_neurons"], 2)
            self.assertEqual(len(record["input_hashes"]), 3)
            with self.assertRaisesRegex(ValueError, "expected reconstruction/wiring"):
                check_prepared(replace(spec, reconstruction_version=1300, wiring_version=1300), data, morphology)

    def test_wrong_cohort_size_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            spec, data, morphology = self.fixture(Path(directory))
            with self.assertRaisesRegex(ValueError, "expected 1352 cells"):
                check_prepared(COLUMN, data, morphology)

    def test_node_nucleus_conflict_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            spec, data, morphology = self.fixture(Path(directory))
            with (data / "nodes.csv").open() as stream:
                rows = list(csv.DictReader(stream))
            rows[0]["nucleus_id"] = "wrong"
            write_csv(data / "nodes.csv", rows)
            with self.assertRaisesRegex(ValueError, "node/nucleus"):
                check_prepared(spec, data, morphology)

    def test_whole_minnie_does_not_accept_proofread_only_preparation(self):
        with tempfile.TemporaryDirectory() as directory:
            spec, data, morphology = self.fixture(Path(directory), WHOLE_MINNIE)
            path = data / "summary.json"
            summary = json.loads(path.read_text()); summary["cohort"] = "proofread"
            path.write_text(json.dumps(summary))
            with self.assertRaisesRegex(ValueError, "cohort all"):
                check_prepared(spec, data, morphology)

    def test_generic_polarity_default_is_preserved_but_can_be_overridden(self):
        argv = ["--data", "data", "--morphologies", "morphologies", "--output", "outputs/test"]
        self.assertFalse(parse_args(argv).inhibition)
        self.assertTrue(parse_args(argv, defaults={"inhibition": True}).inhibition)
        self.assertFalse(parse_args(argv + ["--no-inhibition"], defaults={"inhibition": True}).inhibition)

    def test_all_four_entry_files_describe_sources_without_loading_data(self):
        root = Path(__file__).resolve().parents[1]
        for filename, spec in [("run_three_neurons.py", THREE_NEURONS), ("run_l23.py", L23),
                               ("run_column.py", COLUMN), ("run_whole_minnie.py", WHOLE_MINNIE)]:
            with self.subTest(filename=filename):
                result = subprocess.run([sys.executable, str(root / filename), "--describe"], capture_output=True, text=True)
                self.assertEqual(result.returncode, 0, result.stderr)
                info = json.loads(result.stdout)
                self.assertEqual(info["key"], spec.key)
                self.assertEqual(info["neurons"], spec.neurons)
                self.assertTrue(all(url.startswith("https://") for url in info["sources"]))

    @unittest.skipUnless(importlib.util.find_spec("neuron") and
                         (Path(__file__).resolve().parents[1] / "vendor/aberra2018/mechanisms/x86_64/libnrnmech.so").exists(),
                         "Requires NEURON and compiled published mechanisms")
    def test_edited_contact_csv_cannot_hide_behind_correct_summary_counts(self):
        with tempfile.TemporaryDirectory() as directory:
            _, data, morphology = self.fixture(Path(directory))
            path = data / "edges.csv"
            path.write_text(path.read_text().replace(",3\n", ",2\n"))
            code = """
from dataclasses import replace
import sys
from neuronsim.datasets import COLUMN, launch
launch(replace(COLUMN, neurons=2, contacts=3), sys.argv[1:])
"""
            result = subprocess.run([sys.executable, "-c", code, "--data", str(data),
                "--morphologies", str(morphology), "--output", str(Path(directory) / "out")],
                cwd=Path(__file__).resolve().parents[1], capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("Observed cell/contact counts differ", result.stderr)

    def test_mpi_merge_rejects_mixed_dataset_labels_before_reading_recordings(self):
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory); shards = base / "shards"; run = base / "run"
            shards.mkdir(); run.mkdir()
            (shards / "manifest.json").write_text(json.dumps({"neurons": 2}))
            write_csv(shards / "global_nodes.csv", [{"root_id": root} for root in THREE_NEURONS.roots[:2]])
            (run / "MPI_COMPLETE.json").write_text(json.dumps({"ranks": 2, "shards": str(shards)}))
            # Matching numerical settings, but incompatible dataset identities.
            # Recordings are deliberately absent: rejection must precede reading them.
            common = dict.fromkeys(("duration_ms", "weight_factor", "source_commit", "profile_sha256",
                "global_input_hashes", "global_geometry_manifest_sha256", "stimulus", "dt_ms", "conductance_us",
                "synaptic_reversal_mv", "synaptic_tau1_ms", "synaptic_tau2_ms", "synaptic_delay_ms",
                "placement_seed", "threads", "mpi_ranks"), 0)
            for rank, key in enumerate(("column", "whole_minnie")):
                dest = run / f"rank{rank:04d}" / "long_protocol"; dest.mkdir(parents=True)
                (dest / "summary.json").write_text(json.dumps({**common, "dataset": {"key": key}}))
            with self.assertRaisesRegex(ValueError, "dataset identities differ"):
                merge(run, base / "merged")


if __name__ == "__main__":
    unittest.main()
