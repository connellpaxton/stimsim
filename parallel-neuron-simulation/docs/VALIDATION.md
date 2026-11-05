# Validation status

Reviewed November 2, 2026 on branch `parallel-neuron-simulation`. This page distinguishes checks of the imported/refactored code from earlier large experiments in the separate NeuronSim workspace. The folder is `parallel-neuron-simulation/`; the Python package remains `neuronsim`.

## Checks of this checkout

The latest full suite passed **49 tests with no skips** after the folder rename (`1754908`). The dataset refactor is `f8e1077`; the subsequent critical-path annotations changed comments/documentation only. The documentation review does not change executable code.

| Check | Executed scope and result |
| --- | --- |
| Unit/native checks | Protocols, source integrity, SWC validation/cleanup, published channels, polarity, thread ownership, recording, sharding and dataset-identity guards |
| Original three-neuron entry | All three real cells, 25 ms, one native thread; two soma spikes. Saved soma samples and spike times/IDs matched the earlier coordinate-matched reference exactly; maximum voltage difference 0 mV. |
| L23 entry | Three real L23 cells, 25 ms, two native threads; saved arrays and soma spikes matched the previous runner exactly. This was a subset, not all 313 cells. |
| Prepared-data engine | Synthetic two-cell/five-contact network, 100 ms main run plus 80 ms control/replay. Native recordings remained identical; dataset identity and actual settings appeared in each condition's summary. |
| Real column/full-Minnie inputs | Both passed version, root/nucleus, size and input-hash preflight. Their full networks were not integrated with the new entry files. |
| Folder rename | Python source bytes were unchanged; all four dataset descriptions and the 49-test suite passed from the renamed folder. |

The local logs are in the ignored development virtual environment. Refactor comparisons and preflight records are under ignored `outputs/dataset_entries/`; they are not bundled in Git. [Import validation](IMPORT_VALIDATION.md) separately records the earlier 37-test snapshot.

## Limits and earlier results

The full L23, column and 70,995-cell Minnie experiments were **not rerun for the dataset refactor**. Actual distributed MPI execution was not rerun in this checkout; the local host has no MPI launcher. Sharding and mixed-dataset merge rejection were tested, while prior distributed numerical comparisons are recorded in [MPI construction](MPI_CONSTRUCTION.md).

Earlier full-scale results remain historical evidence, including the [completed one-second full Minnie run](NAS_FULL_RESULTS.md). Their external recordings/source snapshots are indexed in [artifact locations](ARTIFACT_LOCATIONS.md). Passing software checks establishes computational consistency, not biological accuracy or timestep/compartment convergence.

## Re-run local tests

From `parallel-neuron-simulation/`, after installing dependencies and compiling mechanisms:

```sh
NEURON_MODULE_OPTIONS=-nogui OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  .venv/bin/python -m unittest discover -s tests -v
```

Tests explicitly skip native cases when their dependencies/builds are missing; inspect the summary rather than treating a partially skipped suite as this recorded result.
