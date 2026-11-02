# L23 data with the current parallel model

The supplied `SWC_L23/` geometry and `L23_all.csv` contacts completed **1,000 ms of model time** on Genoa with **32 native NEURON threads in one process**. This run retains the current published ModelDB 241165 `cACint209` profile and original-wrapper discretization used in the recent three-neuron comparison. It does not use the different pyramidal-cell profile from the earlier L23 benchmark.

## Inputs and protocol

| Item | Value |
| --- | ---: |
| Supplied L23 SWCs | 349 |
| Valid SWCs simulated | 313 |
| Excluded SWCs with zero-radius points | 36 |
| Contact rows in the full CSV | 42,559 |
| Contacts between simulated cells | 5,112 |
| Directed connected pairs | 4,611 |
| NEURON sections / electrical compartments | 44,303 / 185,035 |
| Duration / timestep / recording interval | 1,000 / 0.025 / 0.1 ms |

No SWCs were repaired or updated to another materialization. The exact source IDs and contact coordinates were retained. Rows referencing unavailable or excluded cells were omitted. All 350 input files, including excluded SWCs, matched local hashes before and after execution.

The most-connected selected source, `864691134941559651`, received **1 nA from 0–1 ms** at its soma. The electrical profile is shared by every cell at 34°C, with initialization at −65 mV. Each contact uses the wrapper's nearest presynaptic and postsynaptic computational segments, threshold 0 mV, positive 0.01 µS conductance, rise/decay 0.5/3 ms, delay 1 ms, and reversal 0 mV. Median contact-endpoint mapping distance was 3.84 µm; maximum was 17.01 µm.

These are the approved published electrical parameters and the original script's all-excitatory synaptic assumptions. The unavailable author-specific channels have not been recreated. No inhibitory population or external input population was added.

## Measured results

| Result | Value |
| --- | ---: |
| Cell loading and synapse setup | 59.98 s |
| Initialization | 2.82 s |
| Native integration, one model second | **105.29 s** |
| Entire run process, including exports and plots | **174.84 s** |
| Resident RAM before plotting | 1.30 GiB |
| Soma spikes / cells with soma spikes | **484 / 307** |
| Saved soma-voltage values | 3,130,313 |

The 32 thread loads ranged from 5,483 to 5,876 compartments and their actual ownership was verified. The full-process time excludes the separate serial validation run. Timings are single measurements.

There were 334 somatic spikes in the first 100 ms and another 150 afterward; the last occurred at 999.45 ms. Activity remained strongly depolarized: 301 somata ended above −40 mV. These outputs do not establish realistic L23 activity under the homogeneous all-excitatory assumptions. The estimated 680,374 synaptic deliveries count local presynaptic threshold events separately for each contact; they are not 680,374 neuronal soma spikes.

## Numerical validation

A fresh **one-thread, 25 ms** run used the same cells, contacts and stimulus. Its **78,563 soma-voltage samples and every soma spike time matched the first 25 ms of the 32-thread run exactly**, with maximum voltage difference **0 mV**. All saved values and final compartment voltages were finite. The full second was not rerun serially.

The runner now records voltage vectors every fixed step with point-process anchors on each owning thread, then downsamples by four. The current `pc.psolve` solver advances the native threads. Cell construction occurs before enabling the requested thread count. The earlier profile remains the runner's default; `--original-family` selects this experiment.

## Inspect and reproduce

- [Network voltage and spike overview](../outputs/NeuronSim-l23-current.9iguVDCM/outputs/parallel32_1s/network_response.png)
- [Recorded arrays](../outputs/NeuronSim-l23-current.9iguVDCM/outputs/parallel32_1s/recordings.npz)
- [Run summary and excluded IDs](../outputs/NeuronSim-l23-current.9iguVDCM/outputs/parallel32_1s/summary.json)
- [Serial/parallel validation](../outputs/NeuronSim-l23-current.9iguVDCM/outputs/validation.json)
- [Activity analysis](../outputs/NeuronSim-l23-current.9iguVDCM/outputs/parallel32_1s/activity_analysis.json)
- [Runner](../tools/run_l23.py)

On Genoa, the isolated workspace is `/home/shisoft/NeuronSim-l23-current.9iguVDCM`. Its retained data, vendor library and environment link to the previously prepared local resources. The saved code and input manifests preserve provenance.

```sh
.venv/bin/python tools/run_l23.py --data data --original-family \
  --max-cells 0 --duration 1000 --amplitude 1 --threads 32 \
  --output outputs/parallel32_1s
```
