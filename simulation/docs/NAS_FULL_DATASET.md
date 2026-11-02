# Full eligible Minnie65 v1300 cohort

The original serial-construction run was replaced by the [validated 96-process MPI workflow](MPI_CONSTRUCTION.md). All cells and contacts were constructed in 91.23 seconds, and the one-second experiment completed successfully. See the [full results](NAS_FULL_RESULTS.md). The original workspace below is retained for preparation provenance and partial serial-construction logs.

The original full-cohort test pipeline was launched on Genoa on October 29, 2026 UTC. Preparation produced **70,995 neurons and 13,309,933 observed internal contacts** across 8,886,062 directed cell pairs. The completed MPI replacement uses the same published ModelDB 241165 profile, E/I synapses and one-second stimulus protocol as the smaller trial, with 96 MPI ranks and 100 ms control/replay checks.

## What “full” includes

The download has 71,234 SWCs. Exactly 71,000 correspond to an unambiguous, uniquely identified neuron with an E/I classification. Five failed soma-position checks and were excluded; this removed 1,499 contacts involving those cells. The retained cells comprise 62,599 excitatory and 8,396 inhibitory neurons.

Only contacts whose two endpoints have a retained SWC are included. The source wiring tables contain other endpoints without an available shape. This cohort is therefore an induced subnetwork of the observed volume, not the entire biological circuit.

There are 1,711 proofread-source and 69,284 dendrite-source trees. Only 1,689 contain a recorded axon. MICrONS [documents that unproofread static exports deliberately exclude axons](https://tutorial.microns-explorer.org/static-repositories.html#neuron-skeletons-v1300); this does not imply axons are absent from the underlying imagery. The simulation generates presynaptic events at the soma and uses the assumed 1 ms delay, so missing axonal branches do not prevent modeled communication. They do prevent interpreting it as measured propagation along a fully reconstructed axon.

## Preparation checks and adaptations

Every retained tree passed source-checksum, root/nucleus position, finite-geometry, topology and radius checks. The prepared cohort has 200,591,690 SWC points. The same documented cleanup repaired 5,414 zero-radius points and collapsed 42,102 identical parent/child aliases while preserving cable length, retained coordinates and positive radii.

The dendrite export labels its root as dendrite type 3. For 69,284 cells, that root was changed to soma type 1 in a working copy only after its position matched the unique-nucleus soma metadata within 0.01 µm and its supplied radius was positive. No axons or surrogate shapes were generated. Original NAS files remain unchanged.

The aggregate wiring uses `-1` where a nucleus annotation is absent: 36,563 retained pairs / 48,254 contacts have at least one such annotation. These rows retain their original fields and use the independently verified v1300 root identity. Every supplied, nonmissing nucleus ID matched metadata. Actual annotated conflicts cause preparation to fail.

The first validation attempts exposed the root-label convention and missing annotation sentinel; their logs remain on Genoa. The successful wiring pass reused checksum-verified prepared morphologies with `--connections-only`. The 33 repository tests pass, including source-integrity, deterministic parallel preparation and guarded root-normalization checks. A normalized dendrite tree also imported and ran 25 ms with the published NEURON profile and finite voltages.

## Monitor and inspect

- Genoa workspace: `/home/shisoft/NeuronSim-nas-full.0UyboUvY`
- NAS experiment: `/mnt/nas-storage/NeuronSim/experiments/minnie65_v1300_full_0UyboUvY`
- Local preparation summary: [summary.json](../outputs/nas_full_trial/prepared/data/summary.json)
- Exact commands: [run_manifest.json](../outputs/nas_full_trial/run_manifest.json)
- Background pipeline: [run_full_trial.py](../outputs/nas_full_trial/run_full_trial.py)

`SIMULATION_STATUS.json` records the pipeline stage on Genoa and NAS. Genoa's `simulation.log` records cell construction and integration progress. The pipeline saves finite-output checks, exact replay comparison, comparison with the previous column's electrical parameters, resource measurements and an archive inventory. It writes `EXPERIMENT_COMPLETE.json` only after successful simulation, validation and archive verification; otherwise its status becomes `failed` with the error preserved. Final full-cohort runtime and memory use are not yet measured.

Prepared shapes/tables are archived before simulation. Results will appear under NAS `results/`; the original download remains under `datasets/minnie65/v1300`. The pipeline is independent of the chat session.
