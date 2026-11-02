# NEURON simulation implementation

NEURON simulation tools for reconstructed mouse cortical networks, contained in the `simulation/` folder of stimsim. The original shared scripts remain in `../neuron-python-share/`.

The `neuronsim` package imports per-cell SWCs, applies a checksum-pinned published ModelDB 241165 electrical profile, constructs native synaptic connections, and runs either native threads or distributed MPI processes. It includes transmission controls, replay checks, input provenance, recording and video tools.

**This is an adaptation of the original model.** Minnie supplies reconstructed anatomy, cell annotations and observed contact counts. The same published `cACint209` electrical profile is applied to every cell; weights, kinetics, delays and most synaptic attachment positions are assumptions. Missing axons are not generated. The exact author-specific mechanisms needed by the original scripts remain incomplete. This repository does not yet contain a fruit-fly or embodied-body implementation.

Start with [data provenance and exact source URLs](docs/DATA_PROVENANCE.md), the [code walkthrough](docs/CODE_WALKTHROUGH.md), and [import notes](docs/IMPORT_NOTES.md).

## Choose a dataset

Use these separate entry files from `simulation/`. Each prints its dataset identity/source URLs and saves `dataset.json` plus dataset metadata in its result summaries. `--describe` prints the source/version/assumption record without loading data or running NEURON.

| Entry file | Dataset and default selection |
| --- | --- |
| [`run_three_neurons.py`](run_three_neurons.py) | Exactly the original three L5 wrapper cells and three coordinate-mapped contacts |
| [`run_l23.py`](run_l23.py) | All 313 valid original L23 cells and 5,112 internal coordinate-mapped contacts |
| [`run_column.py`](run_column.py) | 1,352-cell column: v661 reconstructions linked by nucleus identity to v943 wiring; 124,341 contacts |
| [`run_whole_minnie.py`](run_whole_minnie.py) | Full eligible v1300 cohort: 70,995 cells / 13,309,933 contacts; predominantly dendrite-only trees |

See [commands, exact URLs and connection-method differences](docs/DATASET_ENTRY_POINTS.md). Dataset selection is centralized in `neuronsim/datasets.py`; archive experiments share `neuronsim/archive.py`, while prepared column/full-Minnie experiments share the existing native/MPI engine. Both reuse the same cell importer, published profiles, solver settings and whole-cell thread partitioner.

## Layout

| Path | Contents |
| --- | --- |
| `../neuron-python-share/` | Original files from this repository, retained unchanged |
| `neuronsim/` | Shared simulation package: parameters, SWC cells, contacts, protocols, recording, threading and MPI |
| `tools/` | Preparation, download, analysis, comparison and rendering entry points; see [tool index](tools/README.md) |
| `tests/` | Standard-library unittest checks, including native NEURON checks |
| `vendor/aberra2018/` | Published ModelDB 241165 HOC/MOD sources, licenses and checksum manifests |
| `vendor/m1_260015/` | Sources for the earlier, separate M1 experiments |
| `model/`, `mechanisms/` | Earlier working copies and mechanism-recovery notes |
| `docs/`, `examples/` | Walkthroughs, provenance, historical experiment reports and recording selections |

Datasets, credentials, environments, compiled mechanisms and generated outputs are excluded from new Git additions. Existing binaries already tracked under `../neuron-python-share/` are historical inputs, not the build used by this package.

## Install and compile

Use Linux with Python 3.11 or 3.12, a C++ compiler and `make`. The source simulation was validated with NEURON 9.0.1. From the `simulation/` directory:

```sh
cd simulation
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements-minnie.txt neuron==9.0.1
(cd vendor/aberra2018/mechanisms && ../../../.venv/bin/nrnivmodl)
.venv/bin/python -m neuronsim --help
```

`requirements-neuron.lock.txt` records an earlier environment, not a universal cross-platform lock. The historical M1 workflow has separate `requirements-mouse-m1.txt` dependencies. Rendering requires FFmpeg on PATH; `render_population_video.py` also accepts an installed `imageio-ffmpeg` fallback.

## Run prepared data

Provide `prepared/data/{nodes.csv,edges.csv,summary.json}` and `prepared/morphologies/{manifest.json,swc/}`. Raw downloads must be prepared and checked first; see [NAS dataset preparation](docs/MINNIE_NAS_DATA.md). A full download is about 31.5 GB and is not bundled.

```sh
.venv/bin/python -m neuronsim \
  --data /path/to/prepared/data \
  --morphologies /path/to/prepared/morphologies \
  --output outputs/experiment --threads 4 --inhibition \
  --duration-ms 1000 --validation-window-ms 100
```

`--inhibition` enables E/I polarity; omitting it requests the earlier all-excitatory configuration. This command simulates one second, then runs a 100 ms zero-transmission control and a 100 ms coupled replay. Add `--branch-roots /path/to/roots.json` to record selected cable segments in a single-process run. Full networks can require substantial runtime, RAM and output storage; test a small cohort first.

## MPI execution

Install an MPI implementation such as Open MPI and the optional Python dependency. Run all commands from `simulation/`; the MPI executable and `mpi4py` must use compatible MPI libraries.

```sh
.venv/bin/python -m pip install -r requirements-mpi.txt
.venv/bin/python -m neuronsim.sharding \
  --data /path/to/prepared/data --morphologies /path/to/prepared/morphologies \
  --output outputs/shards96 --ranks 96
mpiexec -n 96 --bind-to core .venv/bin/python -m neuronsim \
  --data /path/to/prepared/data --morphologies /path/to/prepared/morphologies \
  --mpi-shards outputs/shards96 --output outputs/mpi_run \
  --threads 1 --inhibition --duration-ms 1000 --validation-window-ms 100
.venv/bin/python -m neuronsim.merge_mpi \
  --run outputs/mpi_run --output outputs/merged
```

This is **96 processes with one native thread per process**. See [MPI design and prior reference comparisons](docs/MPI_CONSTRUCTION.md). MPI branch-voltage recording is currently unsupported.

## Analyze and render

For a run with the long protocol and both validation conditions:

```sh
.venv/bin/python tools/analyze_nas_trial.py \
  --help
.venv/bin/python tools/render_population_video.py \
  --run outputs/merged/long_protocol \
  --nodes /path/to/prepared/data/nodes.csv \
  --output outputs/video/activity.mp4 \
  --duration-ms 500 --playback-seconds 30 --workers 12
```

The population renderer shows saved soma voltage and spikes at soma positions. It does not invent branch or axon propagation. Analysis of historical full-cohort experiments also requires the provenance/input paths shown by its `--help`; those data are not included in a clone.

## Test

```sh
.venv/bin/python -m unittest discover -s tests -v
```

For an individual module, use `-p 'test_protocols.py'` with the same discovery command. Compile `mechanisms/` separately to enable the earlier xtra/pruning checks:

```sh
(cd mechanisms && ../.venv/bin/nrnivmodl)
```

Native tests skip explicitly when NEURON or required compiled mechanisms are unavailable. Tests use temporary synthetic fixtures and do not require the full Minnie dataset.

## Previous experiments and evidence

The imported reports describe prior runs from the separate NeuronSim workspace, including a 1,352-cell column and a 70,995-cell partial-morphology cohort. Their original NAS/Genoa locations and source versions are documented; local `outputs/...` links require those external artifacts. They do not imply that a new checkout has rerun the experiments. See [documentation index](docs/README.md) and [full-cohort limitations](docs/NAS_FULL_RESULTS.md).
