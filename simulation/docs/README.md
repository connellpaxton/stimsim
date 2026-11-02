# Documentation index

This code was imported from a separate NeuronSim workspace. Read the [import notes](IMPORT_NOTES.md) for source hashes, external-data requirements and the distinction between historical results and checks run in this checkout.

## Read the code

- [Code walkthrough](CODE_WALKTHROUGH.md): execution order, assumptions, file formats, and functions to inspect.
- [Repository setup and commands](../README.md).
- [Real Minnie recording selection](../examples/minnie/README.md).
- [Minnie simulation data on NAS](MINNIE_NAS_DATA.md): file layout, source versions, coordinate conventions and download/integrity tools.

## Current Minnie / ModelDB 241165 experiment

Start with [dataset provenance and validation](DATA_PROVENANCE.md) for exact source URLs, data transformations, model assumptions and limits of the evidence.

| Report | What it explains |
| --- | --- |
| [Completed full Minnie run](NAS_FULL_RESULTS.md) | One-second 70,995-cell result, runtime, aggregate memory, transmission control and replay boundary validation |
| [Parallel MPI construction](MPI_CONSTRUCTION.md) | 96-process construction, exact reference validation, measured speedup and full-run status |
| [Full NAS cohort preparation and run status](NAS_FULL_DATASET.md) | All eligible v1300 shapes, explicit exclusions, missing anatomy and the background full-cohort test |
| [First NAS v1300 trial](NAS_V1300_TRIAL.md) | Version-consistent proofread cohort, preparation, one-second protocol and short validation controls |
| [Three-neuron wrapper comparison](WRAPPER_THREE_NEURON.md) | Exact original SWCs/contacts, matched-location solver agreement, and deviations from column routing defaults |
| [L23 with the current parallel model](L23_CURRENT_PARALLEL.md) | One-second, 313-cell coordinate-based experiment on 32 Genoa threads |
| [One-second propagation video](../MINNIE_PROPAGATION_VIDEO.md) | Video, stimulus, branch recording, validation, runtime, and artifacts |
| [Minnie with inhibitory outputs](../MINNIE_241165_EI.md) | Provenance, E/I contacts, controls, and 32/96-thread agreement |
| [Published-reference adaptation](../MINNIE_241165.md) | Regional electrical profile and earlier all-excitatory reference |
| [ModelDB 241165 review](../MODELDB_241165_REVIEW.md) | Upstream equations, regional parameters, and missing custom variants |

## Earlier experiments and provenance

These retained reports describe configurations distinct from the current Minnie E/I model.

| Report | Scope |
| --- | --- |
| [Minnie with M1 profiles](../MINNIE_COLUMN.md) | Earlier M1 electrical parameters and shape-surrogate comparisons |
| [Mouse M1 benchmark](../MOUSE_M1_BENCHMARK.md) | ModelDB 260015 and generated 1.4-million-receptor network |
| [L23 simulation](../L23_SIMULATION.md) | Downloaded L23 subset with published electrical parameters |
| [L23 Genoa performance](../GENOA_PERFORMANCE.md) | Earlier thread tests and their recording limitations |
| [Single-thread optimization](../SINGLE_THREAD_PERFORMANCE.md) | Static zero-channel pruning and solver comparisons |
| [Memory measurement](../MEMORY_REQUIREMENTS.md) | Four copies of the published example cell |
| [Original archive analysis](../ANALYSIS.md) | Morphologies, contacts, source inventory, and passive demonstration |
| [Working original scripts](../model/README.md) | Runtime/CSV changes and unresolved author-specific mechanisms |
| [Mechanism source search](../mechanisms/SOURCE_SEARCH.md) | Retrieved MOD candidates and unresolved definitions |

Generated files remain in `outputs/` and on Genoa at `/home/shisoft/NeuronSim-minnie-column.WCVsAeuu`. Existing report locations and artifact paths are preserved.
