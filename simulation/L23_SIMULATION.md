# L23 cell and connectivity simulation

For the subsequent one-second run using the current `cACint209` profile and original-wrapper discretization, see [L23 with the current parallel model](docs/L23_CURRENT_PARALLEL.md). The results below retain the earlier L23 pyramidal-cell profile.

The downloaded L23 data includes **cell morphology and anatomical contact records**, not a complete set of fitted electrical models. `SWC_L23/*.swc` describes neuron geometry: 3D points, radii, structural types, and parent links. `L23_all.csv` identifies presynaptic and postsynaptic neuron IDs and contact coordinates.

## What is available

| Item | Count |
| --- | ---: |
| Supplied L23 morphologies | 349 |
| Contact records in `L23_all.csv` | 42,559 |
| Contacts whose two cells are in the supplied L23 folder | 6,366 |
| Morphologies passing tree and positive-radius checks | 313 |
| Contacts between those valid morphologies | 5,112 |
| Distinct directed pairs in the valid network | 4,611 |

The other table rows reference cells outside this supplied L23 set. Thirty-six morphologies contain 44 zero-radius nodes and were excluded without repairing or modifying the downloaded geometry. Affected IDs and reasons are saved in each run's summary.

## Model and assumptions

The adapter in `tools/published_swc.py` imports the downloaded SWCs into NEURON and applies the unchanged published [ModelDB 241165 L23 pyramidal-cell profile](https://github.com/ModelDBRepository/241165/blob/93e01d4a927ab0cc5ca137329b86925777e63f02/cells/L23_PC_cADpyr229_1/biophysics.hoc). Its files and hashes are recorded separately in `vendor/aberra2018/L23_SOURCE_MANIFEST.json`.

This is an active-channel demonstration with explicit assumptions:

- Every selected cell receives the same published L23 pyramidal-cell profile; the dataset does not provide per-cell channel fits or excitatory/inhibitory identity.
- SWC type 3 is treated as basal dendrite and type 4 as apical. The L23 files contain no type-4 labels, so apical compartments are not inferred from geometry.
- Contacts use the downloaded script's `Exp2Syn` defaults: reversal 0 mV, peak conductance 0.01 µS, rise/decay 0.5/3 ms, and 1 ms delay. All modeled synapses are therefore excitatory. These weights and electrical properties are assumptions, not measurements from the CSV.
- Presynaptic and postsynaptic sites map to the closest computational segment. Over the valid full network, median endpoint distance is 4.73 µm, 95th percentile 9.40 µm, and maximum 16.10 µm.
- Temperature is 34°C. Integration is fixed-step at 0.025 ms and defaults to one thread; soma traces are recorded every 0.1 ms. Compartment counts follow the published template's `1 + 2*int(L/40)` rule. No extracellular stimulation or artificial myelination is added.

Coordinate arrays are cached once per cell, and IDs use dictionary lookup. Static zero-conductance pruning is enabled, but this L23 profile has no removable instances of the two channels targeted by the earlier L4 optimization. Nonzero channels remain present.

These choices make the supplied geometry and connectivity runnable. They do not reproduce the missing custom mechanism files exactly or establish a calibrated biological L23 circuit.

## Completed full-network runs

Both runs simulated **313 cells, 5,112 synapses, 44,303 sections, and 143,291 computational segments for 25 ms**. One soma received a 1 ms current pulse, starting at time zero. All recordings were finite.

| Result | Original 1 nA pulse | Evoked 3 nA pulse |
| --- | ---: | ---: |
| Cells with soma spikes | 0 | 299 |
| Total soma spikes | 0 | 303 |
| Estimated synaptic events delivered by 25 ms | 0 | 5,612 |
| Cell loading and synapse setup | 24.3 s | 23.9 s |
| Integration time | 29.8 s | 29.6 s |
| Resident RAM before plotting | 793.2 MiB | 793.2 MiB |

Spikes are upward soma crossings of 0 mV. Delivery counts are inferred from recorded presynaptic threshold events plus the configured delay, counted separately for each contact. The stronger pulse was an explicit additional test, not a parameter inferred from the data. Extensive recruitment in the evoked run reflects this homogeneous excitatory demonstration and should not be interpreted as measured biological behavior.

The four-cell pilot had 13 contacts. With a 3 nA pulse, it produced one stimulated-cell spike and spikes in two recipient cells. A matched zero-conductance synapse control retained the stimulated-cell spike but removed both recipient spikes, confirming functional synaptic transmission in the implementation. Nine automated tests also passed, including regional parameters, coordinate conversion, multi-cell ownership, cell selection, and the existing numerical/connection checks.

## Reproduce and inspect

```sh
# Connected four-cell pilot, original stimulus
.venv/bin/python tools/run_l23.py

# All valid L23 cells, original stimulus
.venv/bin/python tools/run_l23.py --max-cells 0 --output outputs/l23_valid_full

# All valid cells, explicitly stronger evoked stimulus
.venv/bin/python tools/run_l23.py --max-cells 0 --amplitude 3 --output outputs/l23_valid_evoked

# Matched small-network control with no synaptic conductance
.venv/bin/python tools/run_l23.py --max-cells 4 --amplitude 3 --synaptic-weight 0 --output outputs/l23_no_transmission_control
```

`--max-cells 0` selects every valid morphology. The smaller pilot deterministically chooses a neighborhood with many recorded contacts; the stimulated cell is `864691134941559651`.

Each output folder contains a run summary, voltage recordings, soma spike times, and precise contact-to-segment mappings. The evoked full run also has `soma_voltage.png` and `spike_raster.png` for inspection. Both the raw source geometry and all published profile files remain unchanged.

## Multithreading

Use `--threads N` to run the NEURON solver with N threads in one process:

```sh
.venv/bin/python tools/run_l23.py --threads 4 --max-cells 0 --amplitude 3 --output outputs/l23_threads4_balanced
```

Whole cells are assigned largest-first to the thread with the fewest assigned compartments. NumPy/BLAS thread counts remain one; NEURON's solver thread count is controlled separately. The global time recorder is associated with the stimulus point process so NEURON can identify its recording thread. Thread IDs are read through a HOC helper because NEURON 9's direct Python binding coerces `sec_in_thread()` results to Boolean values.

Summaries record the hostname, NEURON version, actual segments per thread, wall time, and process CPU time. `integration_seconds` includes initialization; `initialization_seconds` and `solve_seconds` separate its two parts. Compare thread counts on the same machine with matching parameters, and compare voltage traces and spike times as well as runtime: parallel execution can change numerical results.

The verified local four-thread, 313-cell evoked run assigned 35,966 / 35,775 / 35,774 / 35,776 segments to its threads. Initialization took 21.35 s and solving took 7.74 s, totaling 29.09 s versus 29.60 s in the saved single-thread run. This short run therefore showed little overall improvement. All recorded soma voltages and spike times matched that baseline exactly (303 soma spikes, 5,612 estimated deliveries). Results and comparison are in `outputs/l23_threads4_balanced_verified/`. The earlier two-thread four-cell pilot retained identical spike times but differed by up to 0.0861 mV, so exact agreement is not guaranteed for other partitions.

### Genoa launch from a terminal

Genoa is reachable as `shisoft@192.168.10.239` (`Shisoft-Genoa`), with an AMD EPYC 9654 processor, 96 cores / 192 hardware threads, and 755 GiB RAM. From a terminal with SSH access, run:

```sh
bash tools/run_genoa.sh
```

The launcher packages only L23 input data and required source files, creates a new `~/NeuronSim-genoa.XXXXXXXX` folder on `shisoft@192.168.10.239`, installs pinned dependencies in its own virtual environment, recompiles mechanisms, and benchmarks 1/2/4/8/16/32 solver threads (up to available CPUs). It retrieves recordings, logs, hardware information, and comparisons under local `outputs/NeuronSim-genoa.XXXXXXXX/outputs/`. It requires an existing trusted SSH host key, key authentication, Python 3.11–3.13 with venv support, g++, make, and package-index access on Genoa. It does not install system packages.

`bash tools/run_genoa.sh --prepare-only` prepares the local archive without networking. Remote installation, compilation, and all six benchmarks through 32 threads completed successfully in `/home/shisoft/NeuronSim-genoa.Kd664Nlj`, matching every recorded voltage and spike time. Additional 96-thread runs retained exact spike times but showed voltage-recording differences between repeats. See [GENOA_PERFORMANCE.md](GENOA_PERFORMANCE.md) for measured timings, numerical limitations, and retrieved results. Each setting through 32 ran once; 96 ran twice, so performance timings remain preliminary.
