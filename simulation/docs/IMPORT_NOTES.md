# Simulation code import

This branch adds the simulation implementation developed in the separate NeuronSim workspace to `connellpaxton/stimsim`. The existing `neuron-python-share/` files are preserved unchanged. The imported `model/` directory is a separate earlier working copy, not a replacement for those originals.

## Imported scope

- The `neuronsim` package, including native-thread and MPI execution, source checks, deterministic contact placement, protocols, recording, sharding and merging.
- Preparation/download tools, L23 and three-cell comparisons, analysis and soma/branch video renderers.
- Tests, documentation, recording selections and pinned upstream model sources with their original license and checksum files.
- Earlier M1 and mechanism-recovery work, labeled separately in the documentation index.

The current model is the published ModelDB 241165 adaptation. This import does not add a fruit-fly model, connectome download or embodied simulator. It does not recover the missing custom mechanisms of the original supplied script.

## Source identity

[`IMPORT_MANIFEST.json`](IMPORT_MANIFEST.json) records the source workspace's base commit, its dirty/untracked state, the target base commit and SHA-256 hashes of all 266 copied files. It captures working files beyond the old source commit; referring only to that commit would omit the shared package and MPI implementation. Copied paths are relative to `simulation/`. Its README and ignore rules, the documentation index, and the small root README link were adapted separately. All added implementation code lives under `simulation/`.

For upstream public data and model URLs, see [data provenance](DATA_PROVENANCE.md) and the manifests in `vendor/`. A matching hash establishes byte identity; it does not establish scientific validity.

## External inputs and historical reports

The import excludes environments, credentials, compiled mechanism builds, raw cell datasets and generated recordings/videos. Historical reports retain the original experiment paths. Their `outputs/...` links are references to external artifacts, not files bundled in a fresh clone. NAS and Genoa paths are specific to the original environment and must be replaced when reproducing elsewhere.

Use explicit data/output arguments when running the older tools: some retain the original author's archive location as their default. The main `python -m neuronsim` interface requires explicit prepared data, morphology and output paths. The NAS downloader also has an original-machine default; set `--output` for another location.

The original working-copy manifests under `model/` and historical run source hashes describe earlier snapshots. The import manifest describes the copied source at migration time. Original scientific evidence should be read with its corresponding source version and settings.

## Validation

See [`IMPORT_VALIDATION.md`](IMPORT_VALIDATION.md) for checks run in this checkout. Prior 70,995-neuron timing/replay results are historical measurements, not a claim that this import reran that large experiment. Native mechanisms must be compiled locally; the original Windows DLL is not used by the imported runner.

The branch and any commit remain local until the user explicitly authorizes a push.
