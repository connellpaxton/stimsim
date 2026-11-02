"""Dataset identities and shared launch/preflight logic; no neuron equations here."""

from dataclasses import asdict, dataclass
import csv
import json
import sys

from .data import sha256

MODEL_URL = "https://modeldb.science/241165"
STATIC_URL = "https://tutorial.microns-explorer.org/static-repositories.html"
ARCHIVE_REFERENCE = "https://github.com/connellpaxton/stimsim/tree/5a99045f470db2900d1bf8e9ad5e2ebdf0c04be3/neuron-python-share"
V943 = "https://storage.googleapis.com/mat_dbs/public/minnie65_phase3_v1/v943"
V1300 = "https://storage.googleapis.com/mat_dbs/public/minnie65_phase3_v1/v1300"


@dataclass(frozen=True)
class Dataset:
    key: str
    title: str
    kind: str
    neurons: int
    contacts: int
    default_duration_ms: float
    sources: tuple
    notes: tuple
    layer: str | None = None
    roots: tuple = ()
    available_swcs: int | None = None
    reconstruction_version: int | None = None
    wiring_version: int | None = None

    def metadata(self):
        return {**asdict(self), "electrical_model_url": MODEL_URL}


THREE_NEURONS = Dataset(
    "three_neurons", "Original three-neuron L5 wrapper experiment", "archive",
    3, 3, 25.0, (ARCHIVE_REFERENCE,),
    ("User-supplied neuron_10-26-2026 archive; original SWC download URL/version is unknown.",
     "The URL identifies reference scripts/CSVs, not a public download of these SWCs.",
     "Exact three IDs; supplied contact coordinates; all-excitatory wrapper assumptions; default 1 nA at 0–1 ms."),
    layer="L5", available_swcs=3,
    roots=("864691135640672571", "864691135993429185", "864691136025333561"),
)
L23 = Dataset(
    "l23", "Original archive: all valid L23 neurons", "archive",
    313, 5112, 1000.0, (ARCHIVE_REFERENCE,),
    ("User-supplied SWC_L23/ and L23_all.csv; original SWC download URL/version is unknown.",
     "349 original SWCs; 36 invalid zero-radius trees excluded, without repair.",
     "Supplied contact coordinates; all-excitatory wrapper assumptions; default 1 nA at 0–1 ms."),
    layer="L23", available_swcs=349,
)
COLUMN = Dataset(
    "column", "Minnie reconstructed cortical column", "prepared",
    1352, 124341, 1000.0,
    (STATIC_URL, V943 + "/allen_v1_column_types_slanted_ref_merged.csv.gz",
     V943 + "/connections_with_nuclei.csv.gz",
     "https://storage.googleapis.com/mat_dbs/public/minnie65_phase3_v1/v661/allen_v1_column_types_slanted_ref_merged.csv.gz",
     "https://bossdb-open-data.s3.amazonaws.com/iarpa_microns/minnie/minnie65/skeletons/v661/skeletons/"),
    ("v661 SWCs are linked to v943 annotations/wiring by stable nucleus ID and soma position.",
     "Soma spike sources; area-sampled dendritic targets; E/I enabled by default.",
     "32 selected 5P-IT cells receive 0.5 nA at 10–60 ms; anatomy is not guaranteed complete."),
    reconstruction_version=661, wiring_version=943,
)
WHOLE_MINNIE = Dataset(
    "whole_minnie", "Minnie65 v1300: full eligible partial-morphology cohort", "prepared",
    70995, 13309933, 1000.0,
    (STATIC_URL,
     V1300 + "/connections_with_nuclei.csv.gz?generation=1737139996006427",
     V1300 + "/aibs_cell_info.csv.gz?generation=1752194440283301",
     "https://storage.googleapis.com/microns-static-links/skel/swc/proofread/",
     "https://storage.googleapis.com/microns-static-links/skel/swc/dendrite/"),
    ("All 70,995 eligible retained cells, not every neuron/connection in the source volume.",
     "Only 1,689 trees have recorded axons; no missing axons are generated.",
     "Soma spike sources; area-sampled targets; E/I enabled by default; 0.5 nA at 10–60 ms.",
     "Exact per-SWC URLs/generations are in the download's SOURCE_CATALOG.jsonl."),
    reconstruction_version=1300, wiring_version=1300,
)


def description(dataset):
    """Use the same identity/source text in CLI help, logs and result metadata."""
    return "\n".join((dataset.title, f"Dataset ID: {dataset.key}",
                      *(f"Source: {url}" for url in dataset.sources), *dataset.notes))


def check_prepared(dataset, data, morphologies):
    """Reject the wrong release/cohort before allocating native NEURON cells."""
    summary = json.loads((data / "summary.json").read_text())
    geometry = json.loads((morphologies / "manifest.json").read_text())
    versions = (geometry.get("reconstruction_version"), geometry.get("wiring_version"))
    if versions != (dataset.reconstruction_version, dataset.wiring_version):
        raise ValueError(f"{dataset.key}: expected reconstruction/wiring versions "
                         f"{dataset.reconstruction_version}/{dataset.wiring_version}, got {versions}")
    if (summary.get("selected_neurons") != dataset.neurons
            or geometry.get("ready") != dataset.neurons
            or summary.get("internal_anatomical_synapses") != dataset.contacts):
        raise ValueError(f"{dataset.key}: expected {dataset.neurons} cells / {dataset.contacts} contacts; "
                         "use the generic python -m neuronsim runner for another cohort")
    if dataset.key == WHOLE_MINNIE.key and summary.get("cohort") != "all":
        raise ValueError("whole_minnie requires --cohort all preparation, not the proofread-only subset")
    # Compare exact string IDs. Converting these 18-digit identities to float
    # can silently connect the wrong cells even when cohort counts agree.
    with (data / "nodes.csv").open() as stream:
        nodes = list(csv.DictReader(stream))
    roots = {row["root_id"] for row in nodes}
    records = geometry["records"]
    if (len(nodes) != dataset.neurons or len(roots) != len(nodes)
            or len(records) != len(nodes) or roots != {row["root_id"] for row in records}):
        raise ValueError(f"{dataset.key}: node CSV and morphology identities differ")
    nuclei = {row["root_id"]: row["nucleus_id"] for row in records}
    if any(nuclei[row["root_id"]] != row["nucleus_id"] for row in nodes):
        raise ValueError(f"{dataset.key}: node/nucleus identity differs")
    return {**dataset.metadata(), "input_paths": {"data": str(data), "morphologies": str(morphologies)},
            "input_hashes": {name: sha256(data / name) for name in ("nodes.csv", "edges.csv", "summary.json")},
            "geometry_manifest_sha256": sha256(morphologies / "manifest.json"),
            "reported_neurons": len(nodes), "reported_contacts": dataset.contacts,
            "contact_count_policy": "Preflight checks the preparation summary; native construction/merge verifies actual contacts."}


def _prepared_main(dataset, argv):
    from .cli import parse_args, run_args

    args = parse_args(argv, description=description(dataset), defaults={
        "duration_ms": dataset.default_duration_ms, "threads": 1, "inhibition": True,
    })
    comm = None
    if args.mpi_shards:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
    # Only rank zero reads/hashes the large global inputs. Broadcast the small
    # identity record so 96 ranks do not each parse a full morphology manifest.
    if comm is None or comm.Get_rank() == 0:
        try:
            provenance = check_prepared(dataset, args.data, args.morphologies)
            provenance["run_settings"] = {
                "duration_ms": args.duration_ms, "inhibition": args.inhibition,
                "inhibitory_reversal_mv": args.inhibitory_reversal, "threads_per_process": args.threads,
                "paired_control": args.paired_control, "control_only": args.control_only,
                "validation_window_ms": args.validation_window_ms,
                "mpi_shards": str(args.mpi_shards) if args.mpi_shards else None,
            }
            if args.mpi_shards:
                shards = json.loads((args.mpi_shards / "manifest.json").read_text())
                if (shards["ranks"] != comm.Get_size()
                        or shards["neurons"] != dataset.neurons or shards["contacts"] != dataset.contacts
                        or shards["input_hashes"] != provenance["input_hashes"]
                        or shards["geometry_manifest_sha256"] != provenance["geometry_manifest_sha256"]):
                    raise ValueError("MPI shards were prepared from different inputs or for a different rank count")
            packet = (provenance, None)
        except Exception as error:
            packet = (None, f"{type(error).__name__}: {error}")
    else:
        packet = None
    if comm is not None:
        packet = comm.bcast(packet, root=0)
    provenance, error = packet
    if error:
        raise ValueError(error)
    args.dataset_provenance = provenance
    if comm is None or comm.Get_rank() == 0:
        args.output.mkdir(parents=True, exist_ok=True)
        (args.output / "dataset.json").write_text(json.dumps(provenance, indent=2) + "\n")
        print(description(dataset), flush=True)
    return run_args(args)


def launch(dataset, argv=None):
    """Four small entry files select metadata; two shared engines do the work."""
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv == ["--describe"]:
        print(json.dumps(dataset.metadata(), indent=2))
        return
    if dataset.kind == "archive":
        from .archive import main
        return main(argv, dataset=dataset)
    return _prepared_main(dataset, argv)
