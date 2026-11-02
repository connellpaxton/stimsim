# Parallel cell construction and simulation

The full Minnie65 v1300 experiment now builds cells and synapses using **96 MPI processes, one native NEURON thread per process**. It retains the published ModelDB 241165 parameters, original compartment-sizing rule, exact prepared SWCs and deterministic synapse placements. The previous serial-construction run was stopped after the replacement built successfully and passed reference comparisons. Its inputs and logs are preserved.

## Measured construction times

| Cohort | Serial network construction | MPI cell construction | MPI network construction |
| --- | ---: | ---: | ---: |
| 1,711 cells / 156,945 contacts | 528.85 s | 9.21 s | **9.34 s** |
| 70,995 cells / 13,309,933 contacts | Stopped during construction | **82.52 s** | **91.23 s** |

Times exclude input sharding and MPI launcher startup. The first row measures approximately **57× faster network construction** on the same anatomy and wiring. Numerical integration has a separate cost and now also runs across MPI ranks.

The full cohort contains **10,856,919 electrical compartments**. Rank loads range from 111,833 to 114,012 compartments; each rank owns 738–740 cells. The one-second full simulation completed successfully; see the [full results](NAS_FULL_RESULTS.md) for runtime, memory and validation.

## How ownership works

[`sharding.py`](../neuronsim/sharding.py) assigns whole neurons to ranks, greedily balancing SWC point counts. Global GIDs retain the original sorted-root-ID order. Each incoming synapse is assigned to the rank owning its postsynaptic neuron. Rank input files contain each contact exactly once.

[`mpi.py`](../neuronsim/mpi.py) initializes MPI before constructing `ParallelContext`, then calls the shared cell and synapse builders on each local subset. `pc.gid_connect()` connects local targets to local or remote spike-source GIDs. Live NEURON objects remain in their owning process. This follows NEURON's [documented distributed-network interface](https://nrn.readthedocs.io/en/8.2.2/python/modelspec/programmatic/network/parcon.html).

Every contact retains its offset in the original input stream. The builder advances the same PCG64 stream over skipped contacts, consuming the original choice/position draws for local contacts. This preserves the random placements across process counts. [`merge_mpi.py`](../neuronsim/merge_mpi.py) restores global cell/contact order, verifies coverage, and writes the existing population NPZ/JSON format.

## Numerical validation

A two-rank three-cell test reproduced all placements and spike times; its maximum voltage difference from the earlier local reference was 1.12 × 10⁻¹¹ mV.

The 96-rank, 1,711-cell test reproduced **every one of 156,945 modeled contact rows** exactly. Each of three 100 ms conditions—coupled, zero-transmission and coupled replay—matched **1,712,711 soma-voltage samples and all spike times/IDs exactly**, with maximum voltage difference **0 mV**. The coupled condition had 4,081 spikes; the control had 53. [Comparison evidence](../outputs/mpi_validation/proofread96/results/reference_comparison.json).

The 34 repository tests pass. The new sharding regression verifies destination ownership, unique cell/contact coverage, global contact offsets, stimulus selection and unchanged source fields. The previous one-process CLI remains available when `--mpi-shards` is omitted.

## Run and inspect

Install an MPI implementation and the optional Python dependencies, then prepare shards and launch:

```sh
.venv/bin/python -m pip install -r requirements-mpi.txt
.venv/bin/python -m neuronsim.sharding \
  --data prepared/data --morphologies prepared/morphologies \
  --output shards96 --ranks 96
mpiexec -n 96 --bind-to core .venv/bin/python -m neuronsim \
  --data prepared/data --morphologies prepared/morphologies \
  --mpi-shards shards96 --output mpi_run --threads 1 --inhibition \
  --duration-ms 1000 --validation-window-ms 100
.venv/bin/python -m neuronsim.merge_mpi --run mpi_run --output results
```

`--threads` is the number of native workers **per MPI rank**. MPI branch-voltage recording is currently unsupported. The existing single-process branch-recording workflow remains available.

Genoa's isolated workspace is `/home/shisoft/NeuronSim-mpi.FzV1DzbT`. The dedicated MPI environment reuses the existing NEURON installation through a `.pth` file and installs `mpi4py` separately. The original run's environment and source snapshot were not changed.

- [Full construction timings and rank loads](../outputs/mpi_full_trial/construction.json)
- [Full command and runtime manifest](../outputs/mpi_full_trial/run_manifest.json)
- [Background pipeline](../outputs/mpi_full_trial/run_mpi_full_trial.py)
- NAS experiment: `/mnt/nas-storage/NeuronSim/experiments/minnie65_v1300_full_mpi96_FzV1DzbT`

The pipeline's live `SIMULATION_STATUS.json` and `simulation.log` are on Genoa; status is mirrored to NAS. Rank recordings are merged and checked before a completion marker is written. Memory telemetry sums rank RSS and PSS; early 98-process samples include the two small launcher processes. These are sampled memory observations rather than exact aggregate peaks. The telemetry controller was corrected and resumed without restarting the running MPI simulation.
