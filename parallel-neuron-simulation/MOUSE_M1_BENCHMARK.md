# Mouse M1 network benchmark on Genoa

The selected network completed **100 ms of simulated time** on Genoa using 32 and 96 native NEURON solver threads, one process per run. It contains **903 biophysical neurons, 630 artificial input cells, and 1,400,571 synaptic receptor objects**.

## What the data represents

This uses the published [mouse primary motor-cortex reconstruction, ModelDB 260015](https://modeldb.science/260015), from Dura-Bernal and colleagues (2023). The full published model is described as having over 10,000 neurons and 30 million synapses, constrained by experimental measurements. The supplied morphology, fitted electrical models, and connectivity rules are retained in `vendor/m1_260015/` at commit `45e7d5c4ada0cd32da1a601b3f761fe001a9de39`.

**These are generated model connections, not 1.4 million individually measured synapses.** The selected population scale is 0.09. Some cell models are morphologically detailed; others are reduced compartment models. AMPA and NMDA receptor components can occupy the same modeled connection location.

| Network quantity | Count |
| --- | ---: |
| Biophysical neurons | 903 |
| Artificial input cells | 630 |
| Computational compartments | 113,797 |
| Synaptic receptor objects / synaptic NetCons | 1,400,571 |
| Distinct `(pre, post, section, location)` keys | 724,394 |
| Directed cell pairs, including external inputs | 340,153 |
| AMPA / NMDA receptor objects | 647,189 each |
| GABAA / slow GABAA / GABAB receptor objects | 90,517 / 8,171 / 7,505 |

The previously proposed larger population scale of 0.12 was stopped at the user's request. The completed experiment uses the selected 1.4M-receptor version.

## Protocol and implementation

Runs used NEURON 9.0.1, NetPyNE 1.1.1, Python 3.12.3, 34°C, and fixed 0.025 ms steps on `Shisoft-Genoa` (AMD EPYC 9654, 96 cores / 192 hardware threads). The publication's optional 0.5 nA IT5B current input was enabled; it remained active throughout this 100 ms trial. Published connectivity seeds and parameters were retained. LFP/dipole calculations, publication plotting, and the later conductance-modification protocol were excluded.

All mechanism sources compiled unchanged. Forty-six SOM cells using the compiler-marked non-thread-safe `catcb` channel stayed on thread zero. Other whole cells were balanced across threads by compartment count. NumPy/BLAS threads stayed at one. Voltage recording used explicit cell ownership at every fixed step and exported every fourth step.

## Measured performance

| Measurement | 32 threads | 96 threads |
| --- | ---: | ---: |
| Network construction and counting | 147.27 s | 147.71 s |
| Initialization | 15.48 s | 15.49 s |
| Solve, 100 ms biological time | 24.29 s | 17.11 s |
| Construction + initialization + solve | 187.05 s | 180.31 s |
| Resident RAM after simulation | 3.70 GiB | 3.71 GiB |
| Recorded neuronal spike events | 605 | 605 |
| Recorded artificial-input events | 28 | 28 |

96 threads accelerated solving **1.42×** over 32. Construction dominated total runtime. These timings exclude exporting connection tables and recordings; memory values are snapshots, not peak measurements. Each full-size setting ran once.

## Verification and files

Both runs had the same connectivity SHA-256 (`c523d502369b81230ccc9f55d56f12901a42c7238d9bdcf5c38319bf3f244176`). All **32,032 recorded soma-voltage samples**, recorded cell IDs, time arrays, and recorded spike times/IDs matched exactly. All recorded voltages were finite. The smaller 94-cell pilot also matched exactly on one and four threads; 12 repository tests passed.

The exported connection table was independently checked: 1,400,571 rows, 724,394 location keys, matching receptor counts, valid source/target IDs, finite positive weights, and valid locations/delays. No exported weights were zero.

Local results: [`outputs/NeuronSim-mouse-million.rDKJmxUP/outputs/`](outputs/NeuronSim-mouse-million.rDKJmxUP/outputs/).

- [Comparison and timings](outputs/NeuronSim-mouse-million.rDKJmxUP/outputs/receptors14m_comparison.json)
- [Generated connections, CSV.gz](outputs/NeuronSim-mouse-million.rDKJmxUP/outputs/receptors14m_threads32/connections.csv.gz), [cell metadata](outputs/NeuronSim-mouse-million.rDKJmxUP/outputs/receptors14m_threads32/cells.csv), and [export verification](outputs/NeuronSim-mouse-million.rDKJmxUP/outputs/export_validation.json)
- [Activity plot](outputs/NeuronSim-mouse-million.rDKJmxUP/outputs/activity.png); raw recordings and per-run summaries are in the corresponding thread folders.

Remote workspace: `/home/shisoft/NeuronSim-mouse-million.rDKJmxUP`. Environment versions are pinned in `requirements-mouse-m1.lock.txt`; setup is described in [the source guide](vendor/m1_260015/README.md).

```sh
.venv-m1/bin/python tools/run_mouse_m1.py --scale .09 --threads 96 --duration 100 --evoked --minimum-synapses 1000000 --output outputs/m1_96
```
