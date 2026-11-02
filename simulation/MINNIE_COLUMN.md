# Minnie reconstructed cortical column

For the subsequent experiment using the published model underlying the supplied Python script, see [Minnie with ModelDB 241165](MINNIE_241165.md). This report describes the earlier M1 electrical-profile implementation.

The model uses **each of 1,352 Minnie neurons' own reconstructed branches**, together with the public within-column wiring. Branch lengths, diameters, and topology determine the cable equations solved by NEURON. Electrical-density parameters are transferred from published M1 models; individual Minnie ion-channel conductances and synaptic strengths have not been measured or fitted for this model.

## Observed data

The inputs are the [MICrONS Minnie column annotations](https://www.microns-explorer.org/cortical-mm3) and [archived v943 connectivity tables](https://tutorial.microns-explorer.org/static-repositories.html). The classification and graph use the same materialization version; root IDs stay exact integers or strings. Download URLs and SHA-256 hashes are retained in the data summary. Five annotations without an unambiguous excitatory/inhibitory classification were excluded.

| Quantity | Count |
| --- | ---: |
| Classified neurons | 1,352 |
| Excitatory / inhibitory neurons | 1,188 / 164 |
| Observed internal directed cell pairs | 62,996 |
| Detected anatomical synapses between those cells | 124,341 |
| Simulated receptor objects | 203,738 |
| Retained reconstructed skeleton points | 5,603,234 |
| Imported NEURON sections | 156,468 |
| Electrical compartments | 736,076 |

The archive originally supplied by the user contains 1,097 excitatory-class morphologies: their soma coordinates match the public annotations. It omits the inhibitory population and several other annotated classes. This explains why its validated subset had 998 cells and only 19,929 internal contacts.

The full selected graph is sparse: 3.45% of possible directed cell pairs have a recorded connection. All 1,352 cells belong to one weakly connected component; 1,328 form a strongly connected component with directed paths between every pair. Twenty-two cells have no recorded output to another selected cell. Thus a column can communicate through chains of connections without every neuron contacting every other neuron. Detected wiring can still contain reconstruction errors and incomplete processes; connections outside the selected column are omitted.

## Individual branch geometry

Public [v661 SWC reconstructions](https://tutorial.microns-explorer.org/static-repositories.html) were matched to the v943 wiring by stable nucleus ID and soma coordinates agreeing within 0.01 µm. All 1,352 cells have reconstructions; there is no fallback to a representative shape. The reconstruction manifest records source URLs, raw and import SHA-256 hashes, identity checks, and every repair.

The import copies collapse 1,224 redundant parent/child vertices with identical coordinates, compartment labels, and radii. This preserves branch layout and cable length. There are 362 zero-radius samples in 96 cells; these inherit the nearest measured positive upstream radius. Retained coordinates and measured positive radii are unchanged, and raw downloads remain untouched. Three repaired samples inherit the soma-root radius, which adds local diameter uncertainty. No branches are invented. **512 cells have no recorded axon**; the remaining 840 include axonal reconstruction. Reconstructions can still be incomplete, and skeletons do not explicitly represent every dendritic spine.

NEURON Import3D preserves the recorded compartment labels and branching, using its standard SWC soma representation. This is a cable skeleton, not a membrane surface mesh. Spatial subdivision uses the standard odd-segment `d_lambda=0.1` rule at 100 Hz. Axial resistance and membrane capacitance therefore follow each imported cell's dimensions. All cells report outgoing spikes at the soma's upward 0 mV crossing; observed edges are delivered with modeled delays, without requiring a complete reconstructed axon.

## Electrical and synaptic assumptions

`tools/run_minnie_column.py` instantiates one cell for each observed root ID and preserves each recorded pair's anatomical contact count. It uses unchanged ModelDB 260015 M1 mechanism sources. `tools/export_m1_electrical_profiles.py` samples their input density parameters in a separate process. `tools/reconstructed_cell.py` transfers parameters by soma/axon/basal/apical class and physical path distance, using linear interpolation, median values at coincident reference distances, and endpoint clamping beyond the reference range. M1 shapes are not instantiated in the reconstructed network.

Broad laminar pyramidal classes map to corresponding M1 electrical profiles. Basket cells use PV parameters; Martinotti cells use SOM parameters. Bipolar, neurogliaform, and uncertain inhibitory classes also use generic SOM parameters. These cross-area and subtype mappings are approximations; `cell_models.csv` records every assignment.

Each anatomical contact creates the appropriate AMPA/NMDA or GABA receptor components. The aggregate table has no individual synapse coordinates, so attachment sections are sampled by reconstructed membrane area, with positions sampled within each section and a fixed seed. Excitatory contacts target dendrites; basket-cell contacts target the soma. Other inhibitory contacts use the modeled receptor/target rules. Delays are assumed to be 2 ms plus soma separation divided by 500 µm/ms. Nominal conductance is **0.001 µS (1 nS) per contact**, split between receptor components. These are raw conductances; the M1 EPSP-to-conductance normalization is not applied.

Temperature is 34°C and the fixed timestep is 0.025 ms. Thirty-two manually classified layer-5 IT cells receive 0.5 nA from 10–60 ms. This controlled input substitutes for unknown activity arriving from outside the column. It does not reproduce a recorded in-vivo stimulus condition. The compiler marks `catcb` as non-thread-safe, so the 105 cells using SOM electrical parameters remain on thread zero; other whole cells are distributed across native NEURON threads in one process.

## Completed reconstructed-network simulation

Both runs completed 100 ms on Genoa using 32 native threads and one MPI rank. The control reuses the same constructed cells, contact positions, and injected currents after reinitializing states and setting all transmission conductances to zero.

| Result | Coupled | Zero transmission |
| --- | ---: | ---: |
| Shared construction | 801.41 s | Reused |
| Initialization / setup | 18.39 s | 0.42 s |
| Solve, 100 ms simulated time | 84.28 s | 84.63 s |
| Total neuronal spikes | 152 | 111 |
| Spiking neurons | 54 | 32 |
| Spikes from uninjected neurons | 41 | 0 |
| Excitatory / inhibitory spikes | 121 / 31 | 111 / 0 |
| Resident RAM after simulation | 3.93 GiB | 3.94 GiB |

Twenty-two uninjected neurons fired only with transmission enabled: eight excitatory and fourteen inhibitory cells. Of the 1,320 uninjected neurons, 1,165 showed a maximum soma-voltage difference above 0.01 mV, and 1,141 exceeded 0.1 mV. The largest difference includes an action potential and should not be interpreted as an ordinary EPSP amplitude.

All 1,353,352 recorded soma-voltage samples per condition and all final compartment voltages were finite. Initial soma voltages, recording identities, and time axes matched. This verifies communication within the parameterized network; it does not calibrate the response against in-vivo recordings. Timings are single trials and memory values are resident snapshots, not measured peaks. Thread zero contains 102,980 segments versus about 20,400 per other thread because of the non-thread-safe mechanism, limiting parallel scaling.

Inspect the [simulation summary](outputs/NeuronSim-minnie-column.WCVsAeuu/outputs/reconstructed32/summary.json), [control summary](outputs/NeuronSim-minnie-column.WCVsAeuu/outputs/reconstructed32/zero_transmission_control/summary.json), [voltage comparison](outputs/NeuronSim-minnie-column.WCVsAeuu/outputs/reconstructed32/comparison.json), and [network response plot](outputs/NeuronSim-minnie-column.WCVsAeuu/outputs/reconstructed32/network_response.png).

## Incoming signals on an individual tree

`tools/probe_minnie_branches.py` applies the same isolated 1 nS AMPA event at two modeled apical sites on reconstructed layer-2/3 neuron `864691135782992179`. A separate unstimulated run supplies the baseline. These sites are fixed on the tree when discretization is refined.

| Input site | Soma path distance | Peak local response | Peak soma response |
| --- | ---: | ---: | ---: |
| Proximal apical branch | 60.0 µm | 1.61 mV | 1.01 mV |
| Distal apical branch | 239.6 µm | 15.73 mV | 0.82 mV |

The distal input produces a larger local depolarization but a smaller response at the soma. The two sites also differ in branch diameter and surrounding tree load; this is not a universal distance-only attenuation rule. Responses include active M1 electrical priors on the recorded geometry, rather than isolating geometry experimentally. Halving both `dt` (0.025 → 0.0125 ms) and `d_lambda` (0.1 → 0.05) changes the proximal/distal soma peaks by 0.99%/0.38%. The proximal local peak changes by about 7%; this check does not establish full-network convergence.

See [probe recordings and parameters](outputs/NeuronSim-minnie-column.WCVsAeuu/outputs/branch_probe/summary.json) and [the refinement comparison](outputs/NeuronSim-minnie-column.WCVsAeuu/outputs/branch_probe/convergence.json).

![Branch input and soma response](outputs/NeuronSim-minnie-column.WCVsAeuu/outputs/reconstructed32/branch_signal_transfer.png)

## Earlier representative-shape benchmark

The following results used representative M1 shapes (34,220 segments), before the individual reconstructions were incorporated. They remain a reference for the earlier implementation; they are not timings or responses of the reconstructed network.

| Result | Coupled, 32 threads | Zero transmission, 32 threads | Coupled, 96 threads |
| --- | ---: | ---: | ---: |
| Construction, including attachment-table export | 10.07 s | 10.07 s | 10.11 s |
| Initialization | 2.22 s | 2.23 s | 2.20 s |
| Solve, 100 ms simulated time | 7.94 s | 7.97 s | 6.50 s |
| Total neuronal spikes | 81 | 64 | 81 |
| Spikes from uninjected neurons | 16 | 0 | 16 |
| Excitatory / inhibitory spikes | 65 / 16 | 64 / 0 | 65 / 16 |
| Resident RAM after simulation | 835 MiB | 835 MiB | 849 MiB |

Nine uninjected inhibitory cells fired only when transmission was enabled. Of the 1,320 uninjected cells, 1,152 showed a voltage difference exceeding 0.01 mV relative to the zero-transmission control. This establishes synaptic communication within this parameterized model, not biological calibration of its responses.

The coupled 32- and 96-thread runs produced identical voltage samples and spike times. All 1,353,352 recorded soma-voltage samples and all final compartment voltages were finite. Inputs and modeled attachment tables matched across conditions. Sixteen repository tests passed. Timings are one trial per condition; memory is a snapshot, not a peak measurement.

## Inspect and reproduce

Local files are under [`outputs/NeuronSim-minnie-column.WCVsAeuu/`](outputs/NeuronSim-minnie-column.WCVsAeuu/):

- [Public-data summary and provenance](outputs/NeuronSim-minnie-column.WCVsAeuu/data/column/summary.json), [observed nodes](outputs/NeuronSim-minnie-column.WCVsAeuu/data/column/nodes.csv), and [observed cell-pair counts](outputs/NeuronSim-minnie-column.WCVsAeuu/data/column/edges.csv).
- [Reconstruction manifest](outputs/NeuronSim-minnie-column.WCVsAeuu/data/reconstructions/manifest.json), untouched `raw/` files, validated `swc/` import copies, and [exported electrical priors](outputs/NeuronSim-minnie-column.WCVsAeuu/data/reconstructions/electrical_profiles.json).
- [Recorded branch examples](outputs/NeuronSim-minnie-column.WCVsAeuu/outputs/reconstructed32/recorded_branches.png), [modeled branch signal transfer](outputs/NeuronSim-minnie-column.WCVsAeuu/outputs/reconstructed32/branch_signal_transfer.png), and [full reconstructed-network responses](outputs/NeuronSim-minnie-column.WCVsAeuu/outputs/reconstructed32/network_response.png).
- Earlier surrogate benchmark: [validated comparison](outputs/NeuronSim-minnie-column.WCVsAeuu/outputs/comparison_raw_us.json), per-condition recordings, cell-template assignments, and modeled attachment tables in the `*_us` folders.
- [Measured column positions and connectivity](outputs/NeuronSim-minnie-column.WCVsAeuu/outputs/minnie_anatomy.png) and [earlier surrogate recipient response and raster](outputs/NeuronSim-minnie-column.WCVsAeuu/outputs/minnie_response.png).

Remote workspace: `/home/shisoft/NeuronSim-minnie-column.WCVsAeuu`. Its cached public graph is approximately 3 GB compressed. Earlier prototype folders without the `_us` suffix used an unintended extra weight normalization and are superseded.

With the environment and compiled M1 library described in [the M1 source guide](vendor/m1_260015/README.md):

```sh
.venv-m1/bin/python tools/prepare_minnie.py --cache outputs/minnie_cache --output outputs/minnie_data
.venv-m1/bin/python tools/prepare_minnie_skeletons.py --data outputs/minnie_data --output outputs/minnie_skeletons
.venv-m1/bin/python tools/export_m1_electrical_profiles.py --output outputs/minnie_skeletons/electrical_profiles.json
.venv-m1/bin/python tools/run_minnie_column.py --data outputs/minnie_data --morphologies outputs/minnie_skeletons --output outputs/minnie_reconstructed --threads 32 --paired-control
```

Omit `--morphologies` to reproduce the earlier shape-surrogate mode. `--paired-control` reinitializes the same geometry and placements with all synaptic transmission weights set to zero. A 32-cell validation produced identical control voltages and spikes whether the control reused the coupled construction or ran in a separate process. Twenty repository tests pass, including branch geometry, distinct compartment parameters, lossless duplicate removal, zero-radius repair, and passive rest behavior.

On the retained Genoa workspace, the example branch probe and local plots can be reproduced with:

```sh
.venv/bin/python tools/probe_minnie_branches.py --data data/column --morphologies data/reconstructions --root-id 864691135782992179 --output outputs/branch_probe
.venv/bin/python tools/probe_minnie_branches.py --data data/column --morphologies data/reconstructions --root-id 864691135782992179 --output outputs/branch_probe_fine --d-lambda .05 --dt .0125
python tools/analyze_minnie_reconstructions.py --workspace outputs/NeuronSim-minnie-column.WCVsAeuu
```

The analysis command runs locally from this repository after copying the remote output folder. Downloaded data and simulation outputs are excluded from Git; the scripts and provenance reports remain reviewable source files.
