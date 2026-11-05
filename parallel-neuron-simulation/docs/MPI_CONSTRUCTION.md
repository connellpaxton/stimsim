# Parallel cell construction and simulation

The historical full Minnie65 v1300 experiment completed on October 29, 2026 using **96 MPI processes, one native NEURON thread per process**. It retained the published ModelDB 241165 parameters, original compartment-sizing rule, exact prepared SWCs and deterministic synapse placements. The previous serial-construction run was stopped after the replacement built successfully and passed reference comparisons. Its inputs and logs were preserved. See [current-checkout validation limits](VALIDATION.md) and [external evidence locations](ARTIFACT_LOCATIONS.md).

## Measured construction times

| Cohort | Serial network construction | MPI cell construction | MPI network construction |
| --- | ---: | ---: | ---: |
| 1,711 cells / 156,945 contacts | 528.85 s | 9.21 s | **9.34 s** |
| 70,995 cells / 13,309,933 contacts | Stopped during construction | **82.52 s** | **91.23 s** |

Times exclude input sharding and MPI launcher startup. The first row measures approximately **57× faster network construction** on the same anatomy and wiring. Numerical integration had a separate cost and also ran across MPI ranks.

The full cohort contains **10,856,919 electrical compartments**. Rank loads range from 111,833 to 114,012 compartments; each rank owns 738–740 cells. The one-second full simulation completed successfully; see the [full results](NAS_FULL_RESULTS.md) for runtime, memory and validation.

## How ownership works

[`sharding.py`](../neuronsim/sharding.py) assigns whole neurons to ranks, greedily balancing SWC point counts. Global GIDs retain the original sorted-root-ID order. Each incoming synapse is assigned to the rank owning its postsynaptic neuron. Rank input files contain each contact exactly once.

[`mpi.py`](../neuronsim/mpi.py) initializes MPI before constructing `ParallelContext`, then calls the shared cell and synapse builders on each local subset. `pc.gid_connect()` connects local targets to local or remote spike-source GIDs. Live NEURON objects remain in their owning process. This follows NEURON's [documented distributed-network interface](https://www.neuronsimulator.org/en/latest/progref/modelspec/programmatic/network/parcon.html).

Every contact retains its offset in the original input stream. The builder advances the same PCG64 stream over skipped contacts, consuming the original choice/position draws for local contacts. This preserves the random placements across process counts. [`merge_mpi.py`](../neuronsim/merge_mpi.py) restores global cell/contact order, verifies coverage, and writes the existing population NPZ/JSON format.

## Numerical validation

A two-rank three-cell test reproduced all placements and spike times; its maximum voltage difference from the earlier local reference was 1.12 × 10⁻¹¹ mV.

The 96-rank, 1,711-cell test reproduced **every one of 156,945 modeled contact rows** exactly. Each of three 100 ms conditions—coupled, zero-transmission and coupled replay—matched **1,712,711 soma-voltage samples and all spike times/IDs exactly**, with maximum voltage difference **0 mV**. The coupled condition had 4,081 spikes; the control had 53. [Comparison evidence](ARTIFACT_LOCATIONS.md#mpi-reference-comparison).

At that MPI-development stage, 34 tests passed. The sharding regression verifies destination ownership, unique cell/contact coverage, global contact offsets, stimulus selection and unchanged source fields. The one-process CLI remains available when `--mpi-shards` is omitted. Current tests additionally cover dataset labeling and mismatch guards; their executed scope is recorded in [validation status](VALIDATION.md).

## Run and inspect

From `parallel-neuron-simulation/`, install a compatible MPI implementation and optional dependencies, then prepare shards and launch the full eligible cohort:

Here `prepared/` must contain the v1300 `--cohort all` inputs. For another prepared cohort, use the generic `python -m neuronsim` runner instead of the size/version-guarded `run_whole_minnie.py` entry.

```sh
.venv/bin/python -m pip install -r requirements-mpi.txt
.venv/bin/python -m neuronsim.sharding \
  --data prepared/data --morphologies prepared/morphologies \
  --output shards96 --ranks 96
mpiexec -n 96 --bind-to core .venv/bin/python run_whole_minnie.py \
  --data prepared/data --morphologies prepared/morphologies \
  --mpi-shards shards96 --output mpi_run --threads 1 --inhibition \
  --duration-ms 1000 --validation-window-ms 100
.venv/bin/python -m neuronsim.merge_mpi --run mpi_run --output results
```

`--threads` is the number of native workers **per MPI rank**. MPI branch-voltage recording is currently unsupported. The existing single-process branch-recording workflow remains available.

The original Genoa workspace was `/home/shisoft/NeuronSim-mpi.FzV1DzbT`. Its dedicated MPI environment reused the existing NEURON installation through a `.pth` file and installed `mpi4py` separately. These paths describe the historical run, not a portable environment setup or active job.

- [Full construction timings and rank loads](ARTIFACT_LOCATIONS.md#full-minnie-v1300)
- [Full command and runtime manifest](ARTIFACT_LOCATIONS.md#full-minnie-v1300)
- [Archived completion pipeline](ARTIFACT_LOCATIONS.md#full-minnie-v1300)
- NAS experiment: `/mnt/nas-storage/NeuronSim/experiments/minnie65_v1300_full_mpi96_FzV1DzbT`

The archived `SIMULATION_STATUS.json` records completion; `simulation.log` contains execution progress. Rank recordings were merged and checked before `EXPERIMENT_COMPLETE.json` was written. Memory telemetry summed rank RSS and PSS; early 98-process samples included two small launcher processes. These are sampled observations rather than exact aggregate peaks. During that experiment, the telemetry controller was corrected and resumed without restarting the simulation. No live-job state is implied by this report.
