# Dataset provenance and simulation validation

Recorded October 30, 2026. This note describes the completed **70,995-neuron Minnie65 v1300 experiment and its population video**. Anatomy and connectivity come from MICrONS; electrical behavior comes from a separate published model and explicitly assigned parameters. Successful execution does not establish biological accuracy.

## Public sources and local originals

The source repositories are listed in the [official MICrONS static-export documentation](https://tutorial.microns-explorer.org/static-repositories.html#neuron-skeletons-v1300). The original downloads completed October 28, 2026: **71,240 files, 31,475,900,787 bytes**, including **71,234 unique SWCs**, three compressed tables and their schemas. Raw imagery, segmentation volumes, surface meshes and functional imaging movies were excluded.

Original files and download evidence are stored at:

```text
/mnt/nas-storage/NeuronSim/datasets/minnie65/v1300/
```

The SWCs were downloaded from these documented object prefixes, using the exact per-file URLs in `SOURCE_CATALOG.jsonl`:

| Export | Source prefix | Contents |
| --- | --- | --- |
| Proofread | `https://storage.googleapis.com/microns-static-links/skel/swc/proofread/` | Reconstructed trees that can include axons |
| Dendrite | `https://storage.googleapis.com/microns-static-links/skel/swc/dendrite/` | Trees with axons excluded by the provider |

One source file per root ID was retained, preferring the proofread export when both existed. MICrONS describes these as v1300 skeleton exports; each object's immutable storage generation is also recorded. Neither export guarantees complete biological anatomy.

### Exact downloaded table URLs

| Local file | Generation-pinned download URL |
| --- | --- |
| `connectivity/connections_with_nuclei.csv.gz` | [connections_with_nuclei.csv.gz](https://storage.googleapis.com/mat_dbs/public/minnie65_phase3_v1/v1300/connections_with_nuclei.csv.gz?generation=1737139996006427) |
| `connectivity/connections_with_nuclei_header.csv` | [connections_with_nuclei_header.csv](https://storage.googleapis.com/mat_dbs/public/minnie65_phase3_v1/v1300/connections_with_nuclei_header.csv?generation=1737137741817434) |
| `connectivity/synapses_pni_2_v1_filtered_view.csv.gz` | [synapses_pni_2_v1_filtered_view.csv.gz](https://storage.googleapis.com/mat_dbs/public/minnie65_phase3_v1/v1300/synapses_pni_2_v1_filtered_view.csv.gz?generation=1745555441467935) |
| `connectivity/synapses_pni_2_v1_filtered_view_header.csv` | [synapses_pni_2_v1_filtered_view_header.csv](https://storage.googleapis.com/mat_dbs/public/minnie65_phase3_v1/v1300/synapses_pni_2_v1_filtered_view_header.csv?generation=1745551747386868) |
| `metadata/aibs_cell_info.csv.gz` | [aibs_cell_info.csv.gz](https://storage.googleapis.com/mat_dbs/public/minnie65_phase3_v1/v1300/aibs_cell_info.csv.gz?generation=1752194440283301) |
| `metadata/aibs_cell_info_header.csv` | [aibs_cell_info_header.csv](https://storage.googleapis.com/mat_dbs/public/minnie65_phase3_v1/v1300/aibs_cell_info_header.csv?generation=1752194418928634) |

### Example exact SWC URLs

These illustrate the two downloaded exports; the full catalog identifies every SWC. They are source examples, not a claim that both cells passed final simulation selection.

- Proofread: [864691134884807418.swc](https://storage.googleapis.com/microns-static-links/skel/swc/proofread/864691134884807418.swc?generation=1743799634586523)
- Dendrite: [864691134784325628.swc](https://storage.googleapis.com/microns-static-links/skel/swc/dendrite/864691134784325628.swc?generation=1743808312065723)

`SOURCE_CATALOG.jsonl` records every original URL, generation, expected size and provider checksum. `VERIFIED_FILES.jsonl` and `SHA256SUMS` record local verification. Generation parameters pin downloaded bytes; the URLs can reference files much larger than a browser should open casually.

## Electrical-model source

The electrical model is [ModelDB 241165: Aberra, Peterchev and Grill (2018)](https://modeldb.science/241165), associated with [PubMed 30127100](https://pubmed.ncbi.nlm.nih.gov/30127100/). Source code is pinned to [commit 93e01d4a927ab0cc5ca137329b86925777e63f02](https://github.com/ModelDBRepository/241165/tree/93e01d4a927ab0cc5ca137329b86925777e63f02).

We apply the published [`L4_LBC_cACint209_1/biophysics.hoc`](https://github.com/ModelDBRepository/241165/blob/93e01d4a927ab0cc5ca137329b86925777e63f02/cells/L4_LBC_cACint209_1/biophysics.hoc) regional profile to each imported Minnie tree. The [vendored source manifest](../vendor/aberra2018/SOURCE_MANIFEST.json) records individual model-file URLs and hashes. This transfers a shared electrical profile; it does not provide cell-specific Minnie physiology or reproduce the original supplied script's missing custom mechanisms exactly.

## Selection and changes before simulation

Of 71,234 SWCs, 71,000 had an unambiguous neuron identity and E/I classification. Five failed soma-position checks and were excluded, leaving **70,995 neurons**: 62,599 E and 8,396 I. Retaining only connections between these cells gives **13,309,933 observed contact counts across 8,886,062 directed pairs**. This is an induced subnetwork, not every neuron and connection in the biological volume.

Working copies received documented preparation: 5,414 zero-radius repairs, removal of 42,102 identical parent/child aliases, and root labels changed from dendrite to soma for 69,284 trees only after unique-nucleus position checks. Retained coordinates, positive radii and total cable length were preserved. Original downloads remain unchanged. Per-cell provenance records source/prepared hashes and repairs. See [preparation details](NAS_FULL_DATASET.md).

Only **1,689 retained trees contain a recorded axon**. No missing axons were generated. Presynaptic spikes originate at the soma and cross modeled connections using an assigned delay; they do not travel along a reconstructed missing axon.

## Measured anatomy versus assigned simulation settings

| Item | Origin or assumption |
| --- | --- |
| Skeleton coordinates, estimated radii, cell identities and classifications | Provider reconstructions/annotations, subject to source errors and incompleteness |
| Connected cell pairs and contact counts | Minnie v1300 aggregate connection table |
| Ion channels and regional densities | Shared published `cACint209` profile on every cell |
| Synaptic receptor | NEURON `Exp2Syn`, one object per retained contact count |
| Rise/decay, strength, delay | Assigned 0.5/3 ms, 0.01 µS, 1 ms |
| Reversal potential | Assigned 0 mV for E sources, −80 mV for I sources |
| Postsynaptic location | Dendritic-area sampling, seed 20261027; individual measured contact coordinates were downloaded but **not used in this run** |
| Stimulus | Assigned 0.5 nA into 32 selected somata from 10–60 ms |
| Numerical settings | Fixed 0.025 ms steps, 34 °C, initial −65 mV; compartment rule documented in the results |

The video uses simulated soma voltages and spikes, not experimental voltage recordings. It shows 0–500 ms at 60× slowdown. No branch voltage was recorded for this full-cohort run, and the renderer does not synthesize axonal travel.

## What the checks establish

- **Download integrity:** provider MD5 and byte counts matched the pinned source objects; local SHA-256 records support later comparison. Checksums establish byte identity, not anatomical correctness.
- **Input consistency:** root/nucleus identity, coordinates, topology and radii were checked; all 70,995 raw/prepared SWC hashes were verified after simulation.
- **Parallel implementation:** a smaller reference comparison matched placements, saved voltages and spikes exactly. The full cohort was not repeated serially.
- **Replay:** the full run's 100 ms replay matched all 71,065,995 saved soma-voltage samples exactly. All 102,652 spike events in `[0, 100)` matched exactly. Sixteen cutoff events from the longer run are reported separately; raw recordings were not altered.
- **Transmission control:** during the 100 ms zero-conductance control, only the directly injected cells spiked. Coupled transmission recruited unstimulated cells.

These checks were implemented by the simulation authoring agent and are available for independent review. They establish provenance and computational consistency. Full-network timestep/compartment convergence, cell-type-specific physiological fitting and agreement with experimental activity **have not been demonstrated**. Strong sustained activity should not be treated as proof of realistic cortical dynamics.

## Audit trail and reproduction

Start with the [exact command and runtime manifest](../outputs/mpi_full_trial/run_manifest.json), [validation results](../outputs/mpi_full_trial/results/validation.json), [full results](NAS_FULL_RESULTS.md), and [code walkthrough](CODE_WALKTHROUGH.md). The complete experiment, raw recordings, source snapshots, input links and checksum inventory are archived at:

```text
/mnt/nas-storage/NeuronSim/experiments/minnie65_v1300_full_mpi96_FzV1DzbT/
```

`FILE_MANIFEST.json` identifies archived files; `prepared/morphologies/manifest.json` identifies per-cell source/prepared geometry; `results/modeled_sites.csv.gz` records actual modeled synaptic placements. `video/first_500ms_population.json` records the source-recording and video hashes.

For independent review, trace a cell from its public URL through raw/prepared SWCs to its electrical compartments, trace a connected pair from source count to modeled synapses, then reproduce a small run and test finer temporal/spatial discretization. Those checks address different possible errors and should precede strong biological claims.

Earlier experiments used different inputs: the original archive was supplied locally; the earlier reconstructed column combined v661 SWCs with v943 annotations/wiring through nucleus identity. The separate v1926 completeness/connectivity census did not replace the inputs of this v1300 run. Consult each experiment's manifest rather than assuming every report uses the same dataset version.
