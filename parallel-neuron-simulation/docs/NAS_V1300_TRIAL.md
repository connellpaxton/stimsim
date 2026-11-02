# First simulation from the NAS Minnie dataset

Historical 1,711-cell trial completed October 29, 2026, before the stimsim import. Current-checkout checks are listed in [validation status](VALIDATION.md); evidence files are indexed under [external artifact locations](ARTIFACT_LOCATIONS.md#proofread-v1300-trial).

This trial uses the same published ModelDB 241165 `cACint209` adaptation as the previous column experiment, with **Minnie65 v1300 geometry and wiring from the completed NAS download**. Its cohort is the proofread-export subset across Minnie65, rather than the earlier curated 1,352-cell column or the entire downloaded catalog.

## Selected anatomy

| Input | Count |
| --- | ---: |
| Unique-nucleus, classified proofread neurons | 1,711 |
| Excitatory / inhibitory neurons | 1,492 / 219 |
| Directed connected cell pairs | 85,158 |
| Observed internal synaptic contacts | 156,945 |
| Original / prepared SWC points | 16,649,550 / 16,648,165 |
| Trees containing recorded axon | 1,689 |

All 1,711 candidates passed validation. Root IDs and nucleus IDs agree between wiring and metadata; every soma matches its metadata position within 0.01 µm. The preparation checked source MD5 checksums, finite coordinates, rooted tree structure and supported compartment types.

Metadata assigns 1,509 selected cells to V1, 145 to RL and 57 to AL. Only connections whose two endpoints belong to this cohort are included; external neurons and their background input are absent.

The same preprocessing as the column workflow repaired 1,316 zero-radius points in 224 trees using the nearest positive upstream radius and collapsed 1,385 identical parent/child points in 919 trees. Retained coordinates, positive radii and total cable length were preserved. Original NAS files were not edited. These exported reconstructions are not guaranteed biologically complete.

The broader metadata-eligible catalog contained **71,000 neurons and 13,311,432 internal contacts** across 8,887,257 directed pairs before morphology checks. That larger cohort was outside this trial. A subsequent preparation excluded five cells and produced the **70,995-neuron / 13,309,933-contact** network that [completed the full experiment](NAS_FULL_RESULTS.md).

## Model and protocol

The experiment uses 96 native NEURON threads in one process on Genoa. Each neuron retains its own recorded branches but receives the shared published L4 basket-cell electrical profile. The timestep is 0.025 ms, temperature 34°C and initial voltage −65 mV.

Each observed contact becomes one `Exp2Syn`: rise/decay 0.5/3 ms, 0.01 µS conductance, 1 ms delay, reversal 0 mV for excitatory sources and −80 mV for inhibitory sources. As in the column model, sources detect somatic spikes and targets are sampled by dendritic area using seed 20261027. The downloaded individual-contact positions are retained on NAS but are not used for placement in this run.

The 32 highest-output eligible `5P-IT` excitatory cells receive 0.5 nA somatic current from 10–60 ms. The main simulation lasts 1,000 ms. Two subsequent 100 ms conditions use the same constructed network and current stimulus: a zero-transmission control, followed by a coupled replay for reproducibility checking.

Anatomy and contact counts are data; conductances, delays, kinetics and intrinsic physiology are model assumptions. The unavailable custom mechanisms from the original supplied script remain unresolved.

## Measured results

The trial completed successfully on October 29, 2026 UTC (October 28 locally), using NEURON 9.0.1 and Python 3.12.3 on `Shisoft-Genoa`.

| Measurement | Result |
| --- | ---: |
| Input preparation, measured separately | 204.84 s |
| Network construction | 528.85 s |
| First initialization | 20.04 s |
| Integration of one model second | **387.77 s (6 min 28 s)** |
| Entire simulation process, including both short checks and exports | **1,036.30 s (17 min 16 s)** |
| Sections / electrical compartments | 265,110 / 1,105,118 |
| Contacts from excitatory / inhibitory sources | 67,718 / 89,227 |
| Soma spikes / spiking cells over one second | **6,994 / 1,613** |
| Spikes from directly unstimulated cells | 6,783, across 1,581 cells |

GNU `time` reported a peak resident set of 8,147,636,224 bytes (**7.59 GiB**). The separate `/proc` snapshot after the main solve was 8,258,125,824 bytes (**7.69 GiB**). These counters give slightly different readings; approximately **7.6–7.7 GiB** describes the observed memory use. They are saved separately rather than treating the snapshot as a peak.

Activity spread beyond the 32 injected cells. The first soma spike occurred at 11.4 ms and the last at 239.1 ms; only 25 spikes occurred after 200 ms, with no further spikes through 1,000 ms. All final somata were below −40 mV. This trial produced a transient response, rather than sustained firing for the full second.

The **100 ms zero-transmission control** produced 53 spikes in the 32 injected cells and **zero spikes in unstimulated cells**. The coupled replay produced 4,081 spikes, including 3,937 in 1,534 unstimulated cells. Synaptic transmission therefore changes the response and recruits unstimulated cells in this modeled cohort. These controls do not establish physiological realism.

Every saved voltage, time axis, cell ID and spike in the 100 ms replay matched the main run's prefix **exactly: maximum voltage difference 0 mV**. All recorded voltages and final compartment voltages were finite; actual compartment ownership matched the requested 96-thread partition. Post-run hashes verified all original and prepared SWCs, input tables, geometry manifest and 17 simulation source files. A separate comparison confirmed the electrical profile, solver settings, synaptic parameters and placement seed match the earlier 96-thread column run. This trial was not rerun serially, and timestep convergence was not tested here.

## Saved artifacts

- [Population response and control plot](ARTIFACT_LOCATIONS.md#proofread-v1300-trial)
- [Main run summary](ARTIFACT_LOCATIONS.md#proofread-v1300-trial)
- [Soma voltage and spike arrays](ARTIFACT_LOCATIONS.md#proofread-v1300-trial)
- [Control and replay validation](ARTIFACT_LOCATIONS.md#proofread-v1300-trial)
- [Comparison with the previous column parameters](ARTIFACT_LOCATIONS.md#proofread-v1300-trial)
- [Process resource measurement](ARTIFACT_LOCATIONS.md#proofread-v1300-trial)

The full prepared shapes, tables, source snapshot, published model, logs and results are on NAS at:

```text
/mnt/nas-storage/NeuronSim/experiments/minnie65_v1300_proofread_1711_uwUH4zIS
```

The original workstation's selected recordings and manifests are under `outputs/nas_v1300_trial/experiment/`; this is not an output directory bundled in stimsim. The original download remains at `/mnt/nas-storage/NeuronSim/datasets/minnie65/v1300`.

## Reproduce and inspect

On Genoa, the isolated workspace is `/home/shisoft/NeuronSim-nas-v1300.uwUH4zIS`. It reuses the existing NEURON environment and pinned ModelDB vendor library through symlinks. Prepared inputs and experiment outputs are separate from the original NAS download.

```sh
.venv/bin/python tools/prepare_nas_network.py \
  --source /mnt/nas-storage/NeuronSim/datasets/minnie65/v1300 \
  --cohort proofread --output prepared
mkdir -p results
/usr/bin/time -v -o results/process_resources.txt \
  .venv/bin/python -u -m neuronsim \
  --data prepared/data --morphologies prepared/morphologies \
  --output results --threads 96 --inhibition \
  --duration-ms 1000 --validation-window-ms 100
```

Run from `parallel-neuron-simulation/`. This reproduces the trial's native protocol through the generic runner; the named column/full-Minnie entries intentionally reject this different 1,711-cell cohort. Preparation requires pandas. Use a new output directory; the prepared `data/` contains `nodes.csv`, `edges.csv` and `summary.json`, while `morphologies/manifest.json` records per-cell provenance and repairs.

`tools/analyze_nas_trial.py` is an archive-specific analyzer. It also requires the archived source snapshot and `source_hashes.json`, original/prepared inputs, cell assignments and `results/process_resources.txt`. A fresh CLI run alone does not create that full context; the historical analysis command must not be treated as a general one-command analysis of any run.
