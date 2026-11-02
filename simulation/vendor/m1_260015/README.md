# Published mouse motor-cortex model

`sim/` contains unchanged files from [ModelDB 260015](https://modeldb.science/260015), the mouse M1 model of Dura-Bernal and colleagues (2023). Its full network is described as containing over 10,000 neurons and 30 million synapses constrained by experimental data. This is a computational reconstruction with generated connectivity, not an electron-microscopy contact table.

Source repository: `ModelDBRepository/260015`, commit `45e7d5c4ada0cd32da1a601b3f761fe001a9de39`. `SOURCE_MANIFEST.json` records all downloaded paths, source URLs, sizes, and SHA-256 hashes. Original attribution and copyright notices remain in the source files.

The repository wrapper is `tools/run_mouse_m1.py`. It verifies every source hash, loads the published cell models and connection rules, scales populations, and disables LFP/dipole calculations and publication plotting. It explicitly scales both population density and fixed input counts because NetPyNE 1.1.1 does not consistently apply its general scale factor to both forms.

The wrapper builds on thread zero and uses whole-cell partitions for NEURON's native solver threads. Cells containing mechanisms that the compiler marks as non-vectorized remain on thread zero; for the supplied active populations this affects SOM cells using `catcb`. Channel source and equations are unchanged.

`--evoked` enables the source configuration's optional `IClamp1`: 0.5 nA into the IT5B population, starting at zero. Benchmark durations stay below the publication's mechanism-change protocol at 500 ms. Voltage traces are recorded every fixed step with explicit cell ownership and exported at every fourth step; spike events are sorted by time and cell ID.

Synaptic receptor-object counts, NetCon counts, and distinct `(pre_gid, post_gid, section, location)` keys are reported separately. The latter are modeled connection locations, not measured anatomical counts. `--save-connections` exports the generated receptor-level connection table and cell metadata.

On a separate Python 3.12 environment:

```sh
python3.12 -m venv .venv-m1
.venv-m1/bin/python -m pip install -r requirements-mouse-m1.lock.txt
(cd vendor/m1_260015/sim && ../../../.venv-m1/bin/nrnivmodl mod)
.venv-m1/bin/python tools/run_mouse_m1.py --scale .01 --threads 4 --duration 100
```

Use a separate process from the Aberra model: the two libraries contain overlapping mechanism names.
