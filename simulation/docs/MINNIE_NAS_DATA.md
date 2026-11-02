# Minnie simulation inputs on NAS

The public **Minnie65 v1300** shapes, connectivity and cell metadata are stored at:

```text
/mnt/nas-storage/NeuronSim/datasets/minnie65/v1300
```

The directory is shared with Genoa. Check `DOWNLOAD_STATUS.json` for the current transfer/completion state; a `.part` file is still being assembled or verified.

The initial download **completed on October 28, 2026**: all **71,240 selected source files** (31,475,900,787 bytes) matched their pinned source MD5 checksums and sizes. The final SWC filename/size inventory matched exactly, both wiring tables passed 1,000-row schema samples, and Genoa successfully opened the completed files and 144,120-row cell index. Completion evidence is in `TRANSFER_COMPLETION.json` on NAS.

## Dataset contents

| Path | Content |
| --- | --- |
| `swc/<root_id>.swc` | 71,234 unique original SWCs; proofread export preferred when both catalogs contain the same ID |
| `connectivity/connections_with_nuclei.csv.gz` | Aggregated cell-pair contact counts and pre/post nucleus IDs |
| `connectivity/synapses_pni_2_v1_filtered_view.csv.gz` | Individual contacts, their endpoint positions and cell/supervoxel IDs |
| `connectivity/*_header.csv` | Column names and SQL types for the headerless compressed tables |
| `metadata/aibs_cell_info.csv.gz` | Original cell identity, classification, visual-area and proofreading metadata |
| `metadata/cells.csv.gz` | Derived metadata with headers, stable IDs, soma positions in µm and SWC catalog membership |
| `metadata/cells.sqlite` | The same derived cell information indexed for queries; long IDs are TEXT |
| `metadata/*_table_metadata.json` | Coordinate resolution and table descriptions obtained through CAVE |
| `SOURCE_CATALOG.jsonl` | URLs, pinned object generations, expected sizes and source checksums |
| `VERIFIED_FILES.jsonl`, `SHA256SUMS` | Completed-download integrity records |

The selected source objects total **31,475,900,787 bytes** (31.48 GB), plus small manifests and derived indexes. This includes both connectivity representations. CSVs stay compressed and can be streamed; the original SWC bytes remain unchanged.

Only shapes, tabular wiring and metadata are selected. Image volumes, segmentation voxels, meshes and imaging movies are outside this download.

## Coordinate and identity conventions

SWCs use **micrometers**. Archived table coordinates use **4 × 4 × 40 nm units**, so conversion to micrometers is:

```python
x_um = x * 0.004
y_um = y * 0.004
z_um = z * 0.04
```

Cell-metadata points were exported as 3D EWKB geometry; `index_minnie_metadata.py` decodes them and applies this conversion. A sample of 1,000 downloaded SWC soma positions matched the corresponding unique-nucleus metadata exactly. IDs must remain strings or exact integers throughout CSV parsing. `size` in the synapse table describes detected cleft volume in voxel units, not electrical conductance.

Most static SWCs are dendrite-only; the proofread export includes axons. Neither export guarantees complete reconstruction. The original metadata includes roots associated with multiple nuclei and records without a broad cell-type label. The derived index retains these cases and exposes `root_nucleus_multiplicity` rather than silently treating them as unambiguous neurons. These inputs still need cohort and morphology selection before simulation. Electrical mechanisms, weights, delays and stimuli remain separate model choices.

The official static-export documentation says axons were deliberately removed from unproofread SWCs. Their absence from these files does not imply they are absent from the underlying imagery. The connection tables also contain endpoints outside the downloaded SWC catalog; simulation cohorts retain only connections with two usable endpoints.

`prepare_nas_network.py --cohort all --workers 16 --output prepared` attempts every uniquely identified E/I neuron with a downloaded SWC. The dendrite export labels its root as type 3; the preparation may convert that root to soma type 1 in a working copy only when it matches the unique-nucleus soma position within 0.01 µm and has a positive source radius. It records this normalization separately, preserves measured coordinates/radii, and reports other validation failures explicitly. `--cohort proofread` remains the default. Preparation workers are separate processes and do not change native simulation threading.

## Download and integrity

The downloader uses the public [documented static repositories](https://tutorial.microns-explorer.org/static-repositories.html#neuron-skeletons-v1300). Each source object is pinned to a GCS generation. Size and source MD5 are checked before a final filename is published; local SHA-256 hashes are retained. Large tables use parallel resumable ranges and whole-object verification. This verifies source identity, not biological completeness or positive radii.

```sh
.venv/bin/python tools/download_minnie_simulation_data.py
.venv/bin/python tools/index_minnie_metadata.py
```

The original datasets and scripts in the repository are not replaced. The download is outside Dropbox and Git. See the NAS directory's `README.md`, inventory, status and validation records for the exact files and transfer state.

See the [first v1300 simulation trial](NAS_V1300_TRIAL.md) for the proofread-cohort selection, its observed internal wiring and the commands that prepare these exports for the existing model.

Attribution and license: [MICrONS terms](https://www.microns-explorer.org/terms-and-conditions).
