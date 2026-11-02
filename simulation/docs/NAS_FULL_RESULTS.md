# Full Minnie v1300 partial-morphology simulation results

The **70,995-neuron, 13,309,933-contact** network completed one second of model time on Genoa using **96 MPI ranks with one native thread each**. Controls, replay validation and NAS archive checks completed successfully on October 29, 2026 UTC.

This is the acknowledged partial-morphology experiment. Only 1,689 trees contain a recorded axon; missing axons were not generated. The run retains the published ModelDB 241165 `cACint209` profile on every cell and the previous column's compartment-sizing rule. See [preparation and exclusions](NAS_FULL_DATASET.md) and [MPI construction](MPI_CONSTRUCTION.md).

See [dataset provenance and validation](DATA_PROVENANCE.md) for exact source download URLs, transformations, assigned electrical parameters and what the checks establish.

## Network and protocol

| Quantity | Value |
| --- | ---: |
| Excitatory / inhibitory cells | 62,599 / 8,396 |
| Sections / electrical compartments | 3,934,427 / 10,856,919 |
| Directed connected pairs | 8,886,062 |
| Contacts from E / I sources | 5,596,484 / 7,713,449 |
| Duration / fixed timestep / saved interval | 1,000 / 0.025 / 0.1 ms |

The 32 selected high-output `5P-IT` cells received 0.5 nA at their soma from 10–60 ms. Contacts use `Exp2Syn`, rise/decay 0.5/3 ms, 0.01 µS weight, 1 ms delay and E/I reversals of 0/−80 mV. Presynaptic events originate at the soma; targets are sampled by dendritic area with seed 20261027. These electrical parameters and placements are assumptions, not individual-contact measurements. Raw geometry and wiring both come from v1300.

## Performance

| Measurement | Time |
| --- | ---: |
| Maximum rank cell-construction time | 82.52 s |
| Maximum rank total network-construction time | 91.23 s |
| First initialization | 5.00 s |
| Main integration, one model second | **4,175.29 s — 69 min 35 s** |
| Zero-transmission / coupled-replay integration | 410.34 / 396.48 s |
| Complete MPI execution, including controls and rank exports | **5,097 s — 84 min 57 s** |

These are single measurements. MPI execution excludes input preparation/sharding, the separate output merge, analysis, archive copying, and the earlier abandoned serial-construction attempt.

The maximum sampled aggregate **PSS was 115.00 GiB**, and summed RSS reached **118.82 GiB**. There were 240 complete-job samples, nominally every 15 seconds plus sampling overhead. These are observed sampled maxima, not exact peaks. Early samples include two small launcher processes; subsequent samples cover exactly 96 ranks. GNU `time` around `mpiexec` does not report an aggregate MPI memory peak. The 121.53 GiB sum of post-main rank RSS snapshots was measured at different instants and should not be substituted for a simultaneous measurement.

## Activity and transmission control

The main run recorded **864,908 soma spikes in 65,328 cells**, including **864,552 spikes in 65,296 directly unstimulated cells**. E-labeled cells produced 735,206 spikes and I-labeled cells 129,702. The first event was at 11.4 ms and the last at 999.975 ms; 672,296 spikes occurred after 200 ms. Activity continued after the current pulse and through the end of the run.

During the first 100 ms, the coupled replay produced **102,652 spikes**, including 102,597 in 49,362 unstimulated cells. The zero-transmission control produced **36 spikes in the 32 injected cells**, with **zero unstimulated-cell spikes**. Synaptic transmission therefore recruits the larger modeled network. At the end of the main run, 17,671 somata were above −40 mV. These responses have not been established as realistic mouse cortical dynamics under the shared electrical profile and incomplete anatomy.

## Numerical and provenance validation

All saved voltages and every final compartment voltage were finite. Every contact and cell was accounted for during the MPI merge, and actual thread ownership was verified. The 100 ms replay matched all **71,065,995 soma-voltage samples exactly**, maximum difference **0 mV**, including the sample at 100 ms. All **102,652 spike times and cell IDs in `[0, 100)` ms** also matched exactly.

The initial closed-interval spike comparison exposed an endpoint issue: the longer run records 16 events at `100.00000000011417` ms as it continues; the replay stops at 100 ms. Those cutoff events are now reported separately. Raw recordings were not modified. Three regression tests verify that cutoff events are reported, missing interior events fail validation, and altered voltage at the cutoff fails validation. The full suite passes **37 tests**. The original analysis, failed comparison log and original source hashes are preserved.

Post-run checks verified all **70,995 raw and prepared SWC hashes**, prepared tables and model provenance. Every original source file except the corrected analysis script remained unchanged; both analysis versions and hashes are retained. The electrical parameters match the previous column reference, with MPI execution explicitly identified as the architecture change. The full experiment was not repeated serially and timestep convergence was not tested here.

## Artifacts

- [First 500 ms activity video](../outputs/mpi_full_trial/video/first_500ms_population.mp4) / [HTML player](../outputs/mpi_full_trial/video/first_500ms_population.html): 1080p, 30 fps, 30 seconds of playback (60× slowdown). All 70,995 recorded soma voltages appear at measured positions in two projections; white rings mark spikes in the preceding 2 ms. This shows population recruitment, not voltage traveling along branches or missing axons. No simulation rerun was needed.
- [Population rates and recruitment](../outputs/mpi_full_trial/results/population_rates.png)
- [Voltage and control overview](../outputs/mpi_full_trial/results/network_response.png)
- [Main summary](../outputs/mpi_full_trial/results/long_protocol/summary.json)
- [Validation and endpoint accounting](../outputs/mpi_full_trial/results/validation.json)
- [Memory methodology and observations](../outputs/mpi_full_trial/results/aggregate_memory.json)
- [Original command and source-update provenance](../outputs/mpi_full_trial/run_manifest.json)
- [Completion record](../outputs/mpi_full_trial/EXPERIMENT_COMPLETE.json)

The full archive, including merged recordings and rank outputs, is:

```text
/mnt/nas-storage/NeuronSim/experiments/minnie65_v1300_full_mpi96_FzV1DzbT
```

`results/long_protocol/recordings.npz` contains the full soma/spike arrays; `results/modeled_sites.csv.gz` contains all modeled placements. The archive's `prepared` link references the checksum-verified prepared inputs in the preceding full-cohort workspace. Local `outputs/mpi_full_trial/` contains summaries, figures and provenance; it does not duplicate the large NPZ files. The Genoa workspace remains `/home/shisoft/NeuronSim-mpi.FzV1DzbT`.

The derived video, preview frames and provenance JSON are archived under `video/`. Reproduce with `tools/render_population_video.py` (NumPy, Matplotlib and FFmpeg, or `imageio-ffmpeg`):

```bash
python tools/render_population_video.py \
  --run /path/to/archive/results/long_protocol \
  --nodes /path/to/archive/prepared/data/nodes.csv \
  --output /path/to/new-video/activity.mp4 \
  --duration-ms 500 --playback-seconds 30 --workers 12
```

The renderer samples the nearest saved voltage (maximum timing error 0.050 ms) and converts display values to float32 (maximum change 0.00000382 mV). Source recordings remain unchanged. The provenance JSON identifies the source and encoded video by SHA-256.
