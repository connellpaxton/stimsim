# Genoa multithreading results

On October 27, 2026, the full valid L23 demonstration ran successfully on **Shisoft-Genoa** (`shisoft@192.168.10.239`): AMD EPYC 9654, 96 physical cores / 192 hardware threads, and 755 GiB RAM. An isolated environment used Python 3.12.3, NEURON 9.0.1, and GCC 13.3.0. Pinned dependencies passed `pip check`, and the published mechanisms were compiled on Genoa.

## Workload and measured times

Every run used 313 cells, 5,112 synapses, 143,291 compartments, 25 ms of simulated time, a fixed 0.025 ms step, and a 3 nA / 1 ms pulse. Each ran in one process with the requested NEURON solver threads; cells were assigned by compartment count. NumPy/BLAS threads remained one.

| Threads | Initialization | Solve | Initialization + solve | Solve speedup |
| ---: | ---: | ---: | ---: | ---: |
| 1 | 1.63 s | 49.25 s | 50.88 s | 1.00× |
| 2 | 29.78 s | 24.83 s | 54.61 s | 1.98× |
| 4 | 28.84 s | 12.45 s | 41.30 s | 3.95× |
| 8 | 28.30 s | 6.25 s | 34.56 s | 7.87× |
| 16 | 27.83 s | 3.38 s | 31.21 s | 14.56× |
| 32 | 28.08 s | 1.90 s | 29.98 s | 25.92× |
| 96 | 28.62 s | 1.36 s | 29.98 s | 36.12× |

**96 threads gave the fastest solve, but almost identical total time to 32 threads.** Its solve accelerated 36.1× over one thread and 1.39× over 32 threads. Initialization plus solving accelerated 1.70× over one thread. Loading and network construction took another 46–47 seconds per run: including that work, elapsed time was 96.94 s on one thread, 76.53 s on 32, and 76.48 s on 96. These measurements exclude exporting recordings and plotting. Initialization dominates this short simulation. The table uses one run per setting; an additional 96-thread run took 1.35 s solving and 30.00 s including initialization. Other counts above 32 were not tested.

## Validation and memory

The six runs through 32 threads reproduced identical recorded time arrays, cell IDs, all 78,250 soma-voltage samples, and exact soma spike times. Both 96-thread runs retained all 303 soma spike times and an estimated 5,612 synaptic deliveries, but their recorded voltages differed from the baseline and from each other. The first 96-thread run had a maximum sample difference of 11.93 mV and RMS difference of 0.0917 mV from the baseline; 308 samples differed by more than 0.1 mV. Its repeat had a maximum baseline difference of 11.47 mV and RMS difference of 0.0786 mV. The cause of the voltage-recording variability has not been established, so 96-thread traces should not be treated as equivalent to the baseline.

Thread partitions were populated and their compartment totals checked. The 96-thread partition assigned 1,389–1,705 compartments per thread. Resident RAM before plotting ranged from 791.7 to 804.5 MiB; this is not a peak-memory measurement.

The model retains the electrical assumptions documented in [L23_SIMULATION.md](L23_SIMULATION.md), including a shared published L23 profile and excitatory contact defaults. These results establish execution and agreement for this workload.

## Files and reproduction

Remote workspace: `/home/shisoft/NeuronSim-genoa.Kd664Nlj`.

Retrieved data: [benchmark CSV](outputs/NeuronSim-genoa.Kd664Nlj/outputs/benchmark.csv), [validated detailed results](outputs/NeuronSim-genoa.Kd664Nlj/outputs/validated_benchmark.json), and per-thread recordings, summaries, logs, and plots in the same folder.

From the local repository, `bash tools/run_genoa.sh` creates a fresh remote workspace and repeats the six settings through 32 threads. The additional trials are saved under `threads_96/` and `threads_96_repeat/`, each with a comparison JSON. To repeat 96 threads in the existing remote workspace:

```sh
.venv/bin/python tools/run_l23.py --data data --max-cells 0 --amplitude 3 --duration 25 --threads 96 --output outputs/threads_96_new
```
