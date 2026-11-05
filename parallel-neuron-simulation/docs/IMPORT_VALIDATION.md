# Initial import validation (historical)

This records checks performed November 2, 2026 at initial import commit `fd3080d`. The branch was then named `feat/neuron-simulation` and the folder `simulation/`. Both are now named `parallel-neuron-simulation`. [Validation status](VALIDATION.md) records the later dataset refactor and folder-rename checks; the measurements below remain the original import results.

## Environment and compilation

A fresh virtual environment was installed for this checkout: Python 3.11, NEURON 9.0.1, NumPy 2.4.6, Matplotlib 3.11.2 and pandas 3.1.0. `pip check` reported no broken requirements. Both the published `vendor/aberra2018/mechanisms/` and earlier recovered `mechanisms/` sources compiled successfully using `nrnivmodl`; compiled output is ignored by Git.

At the import, the simulation, sharding, merge, analysis and population-renderer CLI help commands loaded successfully. Source manifests verified all 29 published baseline entries, three additional L23 entries and 111 older M1 entries. Copied implementation files matched their import-manifest hashes at that stage; the documentation index was adapted explicitly. Later refactors are tracked in Git. `IMPORT_MANIFEST.json` is an initial-source snapshot, not a checksum lock for today's implementation.

## Tests

The test command, using today's folder name, is run from `parallel-neuron-simulation/`:

```sh
NEURON_MODULE_OPTIONS=-nogui OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  .venv/bin/python -m unittest discover -s tests -v
```

**37 tests passed, with no skips**, including native NEURON tests, channel/profile checks, synaptic polarity, source-integrity checks, deterministic morphology preparation, MPI input-sharding coverage, branch recordings, original CSV parsing and replay-cutoff regressions.

## End-to-end native simulation

A temporary **synthetic two-cell, five-contact fixture** exercised the real CLI, published profile, native threads, synaptic transmission, protocol reinitialization and recording export. It is not Minnie anatomy and is not biological validation. The main condition lasted 100 ms, followed by an 80 ms zero-transmission control and an 80 ms coupled replay.

| Check | Result |
| --- | --- |
| One versus two native threads | All five contact placements, saved soma voltages and spike times/IDs identical; maximum voltage difference 0 mV |
| Main coupled run | 14 spikes, including seven in the unstimulated cell |
| Zero-transmission control | Eight spikes in the injected cell; none in the unstimulated cell |
| Coupled replay | Inclusive saved voltages and interior spike events matched the main prefix exactly |
| Numerical output | Recorded and final compartment voltages finite |
| Initial relocation check | The CLI simulation passed after moving into the code subdirectory, then named `simulation/`, and matched the one-thread reference exactly |

After the subsequent folder rename, those local smoke fixtures/results are under ignored `parallel-neuron-simulation/outputs/import_smoke/`; environment/build/test logs are local and excluded from Git. Original tracked files under `neuron-python-share/` were verified unchanged. Staged-path checks excluded generated data, environments, caches and compiled binaries. Whitespace checks passed for the shared package and newly authored documentation/configuration; checksum-pinned upstream and historical files retained their original formatting.

## Limits

An actual distributed MPI simulation was not rerun in this checkout: the local host has no MPI launcher. MPI input-sharding tests passed; prior distributed-reference comparisons remain documented separately in [MPI construction](MPI_CONSTRUCTION.md). The full Minnie dataset, older M1 experiment, full-cohort performance benchmark, large-run videos and embodiment were not rerun as part of this import. Historical results are not presented as new validation.
