# Minnie column with ModelDB 241165 and inhibitory outputs

This experiment uses the **1,352-neuron Minnie column's recorded branches, identities, anatomical cell classes, and within-column wiring**, with published ModelDB 241165 electrical parameters. Outputs from Minnie's inhibitory neurons now have an inhibitory reversal potential.

## Data and model provenance

| Input | Source |
| --- | --- |
| 1,188 excitatory and 164 inhibitory neurons | Minnie v943 manual annotations |
| Individual reconstructed branches | Public Minnie v661 SWCs, matched by nucleus identity and soma position |
| 62,996 directed pairs and 124,341 contacts | Minnie v943 within-column connectivity |
| Intrinsic channel mechanisms and regional densities | Published ModelDB 241165 `cACint209` profile |
| Synaptic strengths, kinetics, attachment positions, delays, and current injections | Explicit modeling assumptions |

Thus all neuron-specific structural inputs are from Minnie. The structural dataset does not provide the calibrated conductances or synaptic dynamics used by this simulation. The source versions, limited SWC repairs, and missing axonal reconstructions are documented in [the geometry report](MINNIE_COLUMN.md). The published electrical profile and discretization are documented in [the original-family experiment](MINNIE_241165.md).

One published L4 basket-cell electrical profile is still applied to every tree. Minnie E/I labels determine output polarity; they do not supply distinct fitted intrinsic electrical models for every class. This is an anatomically constrained model rather than a validated physiological reconstruction.

## Inhibition configuration

Every contact remains a deterministic NEURON `Exp2Syn` with positive 0.01 µS (10 nS) conductance, 0.5/3 ms rise/decay, and 1 ms delay. The only change from the previous all-excitatory run is the reversal potential of contacts from inhibitory source neurons:

| Presynaptic Minnie class | Neurons | Contacts | Reversal |
| --- | ---: | ---: | ---: |
| Excitatory | 1,188 | 43,773 | 0 mV |
| Inhibitory | 164 | 80,568 | −80 mV |

The −80 mV assumption agrees with the [published archive's GABA-A reversal default](vendor/aberra2018/mechanisms/ProbGABAAB_EMS.mod); the retained `Exp2Syn` kinetics are the supplied Python script's defaults. Conductance weights stay positive: inhibition comes from the reversal potential. Having both E/I polarities does not by itself establish a balanced circuit.

Recorded geometry, sampled contact locations, conductances, delays, initial voltages, temperature, discretization, and injected currents match the previous run. The short protocol injects one selected neuron with 1 nA from 0–1 ms and runs for 25 ms. The longer protocol injects the same 32 Minnie neurons with 0.5 nA from 10–60 ms and runs for 100 ms. Each has a zero-transmission control. Temperature is 34°C, initialization −65 mV, and the fixed timestep 0.025 ms.

## Completed Genoa results

All four conditions completed on 32 native NEURON threads in one process. The 525,498-compartment network took 210.35 s to construct; construction was shared across conditions.

| Protocol | Previous all-excitatory spikes | Minnie E/I spikes | Zero-transmission spikes |
| --- | ---: | ---: | ---: |
| Short pulse, 25 ms | 1,380 | 54 | 1 |
| Longer drive, 100 ms | 2,681 | 781 | 50 |

The short E/I run recruited 52 uninjected neurons; the longer run recruited 331. No uninjected neuron fired in either zero-transmission control. Inhibitory neurons produced 25 and 253 spikes respectively. E/I runs had 53 and 363 spiking neurons, versus all 1,352 in the corresponding all-excitatory runs.

The short/long E/I solves took **7.26/28.86 s**, using **3.26/3.30 GiB** resident RAM after simulation. Their controls took 7.20/28.95 s. At the end of the short/long E/I runs, 2/14 soma voltages remained above −40 mV, compared with 1,199/600 under the all-excitatory settings. This configuration substantially reduces recruitment and prolonged depolarization; it has not been fitted to measured Minnie electrical activity.

All recorded soma samples and final compartment voltages were finite. The comparison verified that exactly **80,568 contact reversals** changed and every exported contact identity, location, strength, and delay matched the reference. Channel profile, geometry, initial voltage, numerical settings, and stimuli also matched. Both zero-transmission voltage and spike recordings were identical to the corresponding all-excitatory controls. Runtime values are single trials; memory values are snapshots, not peak measurements.

- [Verified polarity comparison](outputs/NeuronSim-minnie-column.WCVsAeuu/outputs/aberra241165_ei32/polarity_comparison.json)
- [Short-protocol inhibition comparison](outputs/NeuronSim-minnie-column.WCVsAeuu/outputs/aberra241165_ei32/original_stimulus/inhibition_comparison.png)
- [Long-protocol inhibition comparison](outputs/NeuronSim-minnie-column.WCVsAeuu/outputs/aberra241165_ei32/matched_protocol/inhibition_comparison.png)
- [E/I voltage responses and spike raster](outputs/NeuronSim-minnie-column.WCVsAeuu/outputs/aberra241165_ei32/matched_protocol/network_response.png)

## 32 versus 96 native threads

The same model was rerun on Genoa with 96 native NEURON threads in one process. Each longer condition simulates **100 ms (0.1 seconds) of model time**; solver times below are wall-clock seconds and exclude setup.

| Condition | Model time | 32-thread solve | 96-thread solve | Speedup |
| --- | ---: | ---: | ---: | ---: |
| Short E/I stimulus | 25 ms | 7.26 s | 4.76 s | 1.52× |
| Short zero transmission | 25 ms | 7.20 s | 4.58 s | 1.57× |
| Longer E/I stimulus | 100 ms | 28.86 s | 17.85 s | 1.62× |
| Longer zero transmission | 100 ms | 28.95 s | 17.95 s | 1.61× |

The longer coupled solve used 38.1% less wall time with 96 threads. Resident RAM was 3.31 GiB versus 3.30 GiB with 32 threads. Construction remained approximately 210.53 s because model loading is serial. The 96-thread partition contained 5,430–5,526 compartments per thread.

Every recorded soma-voltage value and every spike time/neuron ID matched exactly across thread counts in all four conditions; maximum voltage difference was 0 mV. The contact-table contents, assignments, model parameters, and stimuli also matched. The coupled short/long runs therefore retain 54/781 spikes. These are single trials, and memory measurements are resident snapshots rather than peaks.

Results are retained in [`aberra241165_ei96/`](outputs/NeuronSim-minnie-column.WCVsAeuu/outputs/aberra241165_ei96/), with a [validated numerical and performance comparison](outputs/NeuronSim-minnie-column.WCVsAeuu/outputs/aberra241165_ei96/thread_comparison.json). Reproduce the run using the command below with `--threads 96 --output outputs/aberra241165_ei96`. Compare the copied recordings locally with:

```sh
.venv/bin/python tools/compare_minnie_threads.py --baseline outputs/NeuronSim-minnie-column.WCVsAeuu/outputs/aberra241165_ei32 --candidate outputs/NeuronSim-minnie-column.WCVsAeuu/outputs/aberra241165_ei96 --output outputs/NeuronSim-minnie-column.WCVsAeuu/outputs/aberra241165_ei96/thread_comparison.json
```

## Inspect and reproduce

The retained Genoa workspace is `/home/shisoft/NeuronSim-minnie-column.WCVsAeuu`. The prepared graph and individual skeletons are reused without resampling the dataset:

```sh
.venv/bin/python tools/run_minnie_241165.py --data data/column --morphologies data/reconstructions --output outputs/aberra241165_ei32 --threads 32 --inhibition --paired-control --matched-protocol --stimulus-reference outputs/reconstructed32/summary.json
```

The source classification is validated before selecting polarity. `--inhibitory-reversal` can change the assumed negative reversal potential; omitting `--inhibition` retains the earlier all-excitatory behavior.

Local results are under [`outputs/NeuronSim-minnie-column.WCVsAeuu/outputs/aberra241165_ei32/`](outputs/NeuronSim-minnie-column.WCVsAeuu/outputs/aberra241165_ei32/). Each condition contains soma recordings, spikes, and its parameter/provenance summary. Contact locations and actual reversal assignments are exported in `modeled_sites.csv.gz`.

After copying the run folder locally, regenerate comparisons and plots:

```sh
.venv/bin/python tools/analyze_minnie_241165.py --workspace outputs/NeuronSim-minnie-column.WCVsAeuu --run aberra241165_ei32 --reference-run aberra241165_32
.venv/bin/python -m unittest discover -s tests -v
```

The comparison checks every contact against the all-excitatory reference, allowing only inhibitory reversal changes. It also checks profile, geometry, numerical settings, stimuli, and equality of zero-transmission recordings. Twenty-three tests pass, including a native NEURON event test that verifies hyperpolarization with positive inhibitory conductance and rejects unknown E/I classes.
