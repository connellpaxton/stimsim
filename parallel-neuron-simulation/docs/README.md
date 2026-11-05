# Documentation index

The implementation is in `parallel-neuron-simulation/`. Run commands from that folder. Datasets, recorded outputs and historical run workspaces are external to this checkout.

## Current code and usage

- [Dataset entry points](DATASET_ENTRY_POINTS.md): separate three-neuron, L23, column and full-eligible-Minnie launchers, input versions, source URLs and commands.
- [Critical simulation path](CRITICAL_SIMULATION_PATH.md): highlighted blocks from SWC import and channels to the native `pc.psolve()` call.
- [Code walkthrough](CODE_WALKTHROUGH.md): construction, connections, threading, recording and solver orchestration.
- [Setup and requirements](../README.md).
- [Validation status](VALIDATION.md): executed checks of this checkout, small real-data comparisons and full-run/MPI limits.
- [Minnie data on NAS](MINNIE_NAS_DATA.md): the v1300 download layout, coordinate conventions and preparation tools.
- [Dataset provenance](DATA_PROVENANCE.md): exact published URLs, transformations and electrical assumptions for the historical full-Minnie experiment.
- [Recording selection example](../examples/minnie/README.md).
- [Import notes](IMPORT_NOTES.md) / [initial import checks](IMPORT_VALIDATION.md): migration provenance and the earlier validation snapshot.

## Historical completed experiments

These reports preserve measurements from the separate NeuronSim workspace. They do not imply that the current entry files reran every full dataset. [Artifact locations](ARTIFACT_LOCATIONS.md) indexes the external evidence; an ordinary clone does not contain those recordings/videos.

| Report | Scope |
| --- | --- |
| [Full Minnie results](NAS_FULL_RESULTS.md) | Completed one-second 70,995-cell run, runtime, sampled memory, control and replay-cutoff accounting |
| [MPI construction](MPI_CONSTRUCTION.md) | Prior 96-process benchmarks/reference comparisons and the current distributed execution interface |
| [Full-cohort preparation](NAS_FULL_DATASET.md) | v1300 selection, exclusions, root normalization and the replaced serial attempt |
| [First v1300 trial](NAS_V1300_TRIAL.md) | Separate 1,711-cell proofread cohort and completed one-second result |
| [Three-neuron comparison](WRAPPER_THREE_NEURON.md) | Historical four-condition wrapper comparison and the effect of changing contact locations |
| [L23 parallel result](L23_CURRENT_PARALLEL.md) | Historical full 313-cell, one-second experiment on 32 native threads |
| [Column branch video](../MINNIE_PROPAGATION_VIDEO.md) | Earlier column stimulus, branch recordings, video and validation |
| [Column E/I results](../MINNIE_241165_EI.md) | Earlier column conditions and 32/96-native-thread agreement |
| [Published-reference adaptation](../MINNIE_241165.md) | Earlier all-excitatory column reference |
| [ModelDB 241165 review](../MODELDB_241165_REVIEW.md) | Upstream equations, regional parameters and missing custom variants |

## Earlier model alternatives and source analysis

These retained reports use configurations distinct from the current four dataset entries. Their measurements, commands and environment paths describe the experiments at the time.

| Report | Scope |
| --- | --- |
| [Minnie with M1 profiles](../MINNIE_COLUMN.md) | Earlier M1 physiology and shape-surrogate comparisons |
| [Mouse M1 benchmark](../MOUSE_M1_BENCHMARK.md) | ModelDB 260015 and generated receptor/connectivity experiment |
| [Earlier L23 simulation](../L23_SIMULATION.md) | Older L23 electrical profile and subset |
| [Earlier Genoa performance](../GENOA_PERFORMANCE.md) | Older thread tests and recording limitations |
| [Single-thread optimization](../SINGLE_THREAD_PERFORMANCE.md) | Static channel pruning and solver comparisons |
| [Four-cell memory experiment](../MEMORY_REQUIREMENTS.md) | Four copies of the bundled published cell |
| [Original archive analysis](../ANALYSIS.md) | Morphologies, contacts and passive demonstration |
| [Working original scripts](../model/README.md) | Earlier CSV/runtime changes; custom mechanisms remain incomplete |
| [Mechanism search](../mechanisms/SOURCE_SEARCH.md) | Retrieved candidates and unresolved definitions |
