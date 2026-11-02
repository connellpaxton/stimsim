# Three-neuron wrapper comparison

The exact dataset selected in `wrapper-synaptic-connections.py` now runs through the current infrastructure. It uses these original `SWC_L5/` files, without replacing their geometry or IDs with newer materializations:

1. `864691135640672571.swc` — injected cell
2. `864691135993429185.swc`
3. `864691136025333561.swc`

Only three rows of `L5_all.csv` have both endpoints in this set: **1 → 2 → 3 → 1**, one contact per edge. The stimulus is **1 nA from 0–1 ms**, with **25 ms** simulated at **dt = 0.025 ms**, starting at **−65 mV**. There are 2,708 electrical compartments. Each condition runs in one process with one native thread.

## Comparison and results

Both implementations use the same published ModelDB 241165 `cACint209` profile at 34°C. The original author-specific channel variants remain unavailable, so this compares our approved published-model adaptation, not the unavailable original electrical model. The wrapper does not explicitly set temperature.

| Condition | Neuron 1 spikes | Neuron 2 spikes | Neuron 3 spikes | Maximum axon voltage difference from reference |
| --- | --- | --- | --- | ---: |
| Wrapper reference, nearest contact segments | 0.900 ms | 6.025 ms | None | Reference |
| Current solver, same contact segments | 0.900 ms | 6.025 ms | None | **0 mV** |
| Current column defaults, soma sources and area-sampled targets | 0.900 ms | None | None | 108.702 mV |
| Current defaults with transmission disabled | 0.900 ms | None | None | 114.568 mV |

The matched-location run reproduces every recorded soma/axon sample, synaptic conductance sample, spike time, and final compartment voltage exactly. It uses the current cached geometry, zero-conductance pruning, GID routing and fixed-step `pc.psolve` integration. The reference executes the repository working wrapper with a published-cell adapter and cached coordinate lookup replacing inactive `xtra` coordinate storage.

Changing the connection assumptions changes firing: neuron 2 peaks at −53.938 mV under current defaults, rather than producing the reference spike. Neuron 3 consequently receives no onward synaptic response. The 108.702 mV maximum difference includes a spike being present in one trace and absent in the other; it is not solver error. Source and target assumptions change together here, so their individual contributions are not isolated. The area-sampled result is one deterministic placement realization.

## Inspect and reproduce

- [Voltage comparison plot](../outputs/wrapper_three_neuron/comparison.png)
- [Complete results and input hashes](../outputs/wrapper_three_neuron/comparison.json)
- [The three original contact rows](../outputs/wrapper_three_neuron/selected_contacts.csv)
- [Generated experiment report](../outputs/wrapper_three_neuron/README.md)
- [Runner and comparison checks](../tools/compare_wrapper_three.py)

Each condition directory contains `recordings.npz`, `voltage.csv`, and `summary.json`. All original source hashes were checked after execution. The existing 26 tests pass. Timing definitions differ between the wrapper and engine runs, so these runs are not a speed benchmark.

```sh
.venv/bin/python tools/compare_wrapper_three.py
```

The tool verifies the original wrapper's cell selection, CSV filename, stimulus, timestep, initialization and stop time before executing. Existing column simulation inputs and source scripts are unchanged.

## One-second runtime measurement

Extending only the stop time to 1,000 ms took **24.983 s** of native integration for the wrapper adaptation and **25.160 s** for the current matched-location engine. Setup/execution including recording export was about **25.9 s** for each, excluding initial module imports. These are single sequential measurements on the same machine with one thread, rather than a demonstrated speed difference.

All 40,001 time samples, soma/axon/source traces, synaptic conductances/currents, final compartment voltages and spike times matched exactly. The first 25 ms also exactly matched the earlier experiment. See the [one-second timing report](../outputs/wrapper_three_neuron_1s/README.md) and [validation data](../outputs/wrapper_three_neuron_1s/timing_comparison.json). The source wrapper is unchanged; timing instrumentation and the stop-time override are applied only in memory.
