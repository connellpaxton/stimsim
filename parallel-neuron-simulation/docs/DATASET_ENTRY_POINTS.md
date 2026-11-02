# Dataset-specific simulation entry points

Run commands from `parallel-neuron-simulation/` after installing dependencies and compiling the published mechanisms as described in [the README](../README.md). No entry point downloads data automatically or silently selects a different release. Use a new output directory for each experiment.

## Select the experiment

| Python file | Required input | Full retained cohort | Default protocol |
| --- | --- | --- | --- |
| [`run_three_neurons.py`](../run_three_neurons.py) | Original archive `SWC_L5/` and `L5_all.csv` | Three exact wrapper IDs, three contacts | 25 ms; 1 nA into cell 1 at 0–1 ms; all contacts excitatory |
| [`run_l23.py`](../run_l23.py) | Original archive `SWC_L23/` and `L23_all.csv` | 349 supplied SWCs → 313 valid cells, 5,112 contacts | 1,000 ms; 1 nA into the highest-output selected cell at 0–1 ms; all contacts excitatory |
| [`run_column.py`](../run_column.py) | Prepared column tables and reconstructions | 1,352 cells, 124,341 contacts | 1,000 ms; 0.5 nA into up to 32 selected 5P-IT cells at 10–60 ms; E/I enabled |
| [`run_whole_minnie.py`](../run_whole_minnie.py) | Prepared v1300 `--cohort all` data | 70,995 cells, 13,309,933 contacts | Same current protocol as the column; E/I enabled |

All four use the published ModelDB 241165 `cACint209` adaptation, original-wrapper compartment rule, fixed 0.025 ms timestep and 34 °C. This does not recreate missing custom channels or supply cell-specific physiology. Weights, kinetics, delays and current injections remain assumptions.

Each file supports `--help` and standalone `--describe`. For example:

```sh
.venv/bin/python run_whole_minnie.py --describe
```

`--describe` reports the expected snapshot; it does not validate local files. Actual runs print the identity and URLs, save `dataset.json`, and embed the same provenance in `summary.json`. Prepared entry points check release versions, cohort counts, root/nucleus identities and input hashes before native construction. The full-Minnie entry rejects proofread-only preparation. These guards prevent common mix-ups; they do not establish biological completeness or authenticity independently of the underlying source manifests. Other cohorts remain available through the generic `python -m neuronsim` interface.

## Original three-neuron and L23 datasets

These SWCs came from the locally supplied `neuron_10-26-2026` archive. Its original SWC download URL and materialization version are **not established**. The [stimsim reference scripts and CSVs at the original commit](https://github.com/connellpaxton/stimsim/tree/5a99045f470db2900d1bf8e9ad5e2ebdf0c04be3/neuron-python-share) are available online; that repository is **not** a public source for the missing SWC archive. If a known archive URL becomes available, pass `--source-url https://...` to record it.

The three fixed root IDs are `864691135640672571`, `864691135993429185`, and `864691136025333561`. The runner requires those files and their three-contact cycle; another L5 trio is rejected.

```sh
archive=/path/to/neuron_10-26-2026
.venv/bin/python run_three_neurons.py --data "$archive" \
  --output outputs/three_neurons --threads 1
.venv/bin/python run_l23.py --data "$archive" \
  --output outputs/l23 --threads 32
```

L23 defaults to every valid cell. `--max-cells 3 --duration-ms 25` requests an explicit small subset, recorded as `subset: true`; it is not reported as a full L23 simulation. Invalid zero-radius SWCs are excluded, without repair. The full parent cohort is checked even when requesting a subset.

The shared archive engine attaches both ends of each contact to the nearest computational segment using the supplied nanometer coordinates. It preserves source row order and converts coordinates to SWC micrometers. Its presynaptic source may be an axon segment. This reproduces the coordinate-based assumptions, rather than replacing them with the prepared-model soma/GID routing.

## Reconstructed Minnie column

This is the previous 1,352-neuron column, not the separate 1,711-cell v1300 proofread trial. Sources:

- [v943 column annotations](https://storage.googleapis.com/mat_dbs/public/minnie65_phase3_v1/v943/allen_v1_column_types_slanted_ref_merged.csv.gz).
- [v943 aggregate connectivity](https://storage.googleapis.com/mat_dbs/public/minnie65_phase3_v1/v943/connections_with_nuclei.csv.gz).
- [v661 identity annotations](https://storage.googleapis.com/mat_dbs/public/minnie65_phase3_v1/v661/allen_v1_column_types_slanted_ref_merged.csv.gz).
- SWC prefix: `https://bossdb-open-data.s3.amazonaws.com/iarpa_microns/minnie/minnie65/skeletons/v661/skeletons/`. Each filename is `<v661_root_id>_<nucleus_id>.swc`; every exact source URL is retained in the prepared morphology manifest.

Preparation links v661 geometry to v943 wiring through stable nucleus IDs and soma-position checks. It does not assume segmentation root IDs are interchangeable across releases.

```sh
.venv/bin/python tools/prepare_minnie.py \
  --cache outputs/column_cache --output outputs/column/data
.venv/bin/python tools/prepare_minnie_skeletons.py \
  --data outputs/column/data --output outputs/column/morphologies
.venv/bin/python run_column.py \
  --data outputs/column/data --morphologies outputs/column/morphologies \
  --output outputs/column_run --threads 32 --validation-window-ms 100
```

## Full eligible Minnie65 v1300 cohort

The [official MICrONS static-repository documentation](https://tutorial.microns-explorer.org/static-repositories.html#neuron-skeletons-v1300) identifies the sources. The download uses these generation-pinned tables:

- [v1300 connections](https://storage.googleapis.com/mat_dbs/public/minnie65_phase3_v1/v1300/connections_with_nuclei.csv.gz?generation=1737139996006427).
- [v1300 cell metadata](https://storage.googleapis.com/mat_dbs/public/minnie65_phase3_v1/v1300/aibs_cell_info.csv.gz?generation=1752194440283301).
- SWC prefixes: `https://storage.googleapis.com/microns-static-links/skel/swc/proofread/` and `https://storage.googleapis.com/microns-static-links/skel/swc/dendrite/`. `SOURCE_CATALOG.jsonl` retains each exact URL, object generation and source checksum. [Full provenance](DATA_PROVENANCE.md) also lists the individual-contact table downloaded but not used for this run's placement.

```sh
raw=/path/to/minnie65/v1300
.venv/bin/python tools/download_minnie_simulation_data.py --output "$raw"
.venv/bin/python tools/index_minnie_metadata.py --data "$raw"
.venv/bin/python tools/prepare_nas_network.py --source "$raw" \
  --cohort all --workers 16 --output outputs/full_prepared
```

For the established 96-process execution:

```sh
.venv/bin/python -m neuronsim.sharding \
  --data outputs/full_prepared/data --morphologies outputs/full_prepared/morphologies \
  --output outputs/full_shards96 --ranks 96
mpiexec -n 96 --bind-to core .venv/bin/python run_whole_minnie.py \
  --data outputs/full_prepared/data --morphologies outputs/full_prepared/morphologies \
  --mpi-shards outputs/full_shards96 --threads 1 \
  --output outputs/full_run --validation-window-ms 100
.venv/bin/python -m neuronsim.merge_mpi --run outputs/full_run --output outputs/full_results
```

Install the optional MPI requirements first. This is 96 MPI processes with one native thread each. A single-process run omits `--mpi-shards`. Rank zero performs global identity preflight and broadcasts its provenance; every rank summary and the merged summary retain it. MPI merging rejects inconsistent dataset identities. Branch-voltage recording remains single-process only.

“Whole Minnie” here means the full **eligible retained cohort**, with only connections between retained cells. Just 1,689 trees contain recorded axons. No missing axons or external neurons are invented. Both prepared experiments use soma spike sources and area-sampled dendritic targets, not individual measured synaptic attachment sites.

## Shared code and validation

[`neuronsim/datasets.py`](../neuronsim/datasets.py) owns identities, URLs and prepared preflight. [`archive.py`](../neuronsim/archive.py) handles both coordinate-contact archives; `tools/run_l23.py` is a compatibility shim retaining its earlier defaults. The existing [`simulation.py`](../neuronsim/simulation.py) / [`mpi.py`](../neuronsim/mpi.py) handle both prepared cohorts. All reuse `cell.py`, published-model loading, `create_synapse()`, `configure_solver()` and `partition_cells()`; the four entry files contain no copied neuron equations or integration loops. Shared model constants live in `config.py`.

See [validation status](VALIDATION.md) for executed tests, exact small-run comparisons and the limits of the real prepared-input preflight. Historical full-scale results are separate from checks of these new entry files.
