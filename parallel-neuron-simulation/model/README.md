# Supplied NEURON model: working copy

These three Python scripts were copied from the downloaded `neuron_10-26-2026` archive. `source-manifest.json` records hashes of the original files and the initial working copies. The download itself was not edited.

Changes are limited to runtime setup and CSV parsing:

- The wrapper resolves sample SWCs in `SWC_L5/` and connectivity in the data directory. Set `NEURONSIM_DATA` to use another archive location.
- The unused pandas import was removed from the wrapper. `cell.py` uses standard-library `csv.DictReader` instead of pandas for connectivity. IDs remain strings and position coordinates are converted to floats before segment lookup.
- Animation defaults to off for headless execution. Set `NEURONSIM_ANIMATION=1` when running the wrapper directly with a display available.
- `plotting3D.py` is unchanged apart from newline normalization.

The channel insertions, channel parameters, discretization, stimulation, synapse settings, and integration settings are unchanged. The optional-parameter bug documented in the analysis is also unchanged; the wrapper's constructor call avoids it.

Run the actual wrapper through the launcher, which checks mechanisms and exports its recorded vectors and figure:

```sh
.venv/bin/python tools/run_original.py --check
.venv/bin/python tools/run_original.py --mechanisms /path/to/compiled/mechanisms
```

To invoke the wrapper directly:

```sh
NEURON_MODULE_OPTIONS=-nogui MPLBACKEND=Agg .venv/bin/python model/wrapper-synaptic-connections.py
```

The model still needs the author's custom ion-channel mechanisms. A separately retrieved and tested `xtra.mod` candidate is in `../mechanisms/`; see its README for source and compatibility limits. No ion-channel mechanisms have been substituted.
