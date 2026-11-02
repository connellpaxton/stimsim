# Four-neuron memory requirements

The measured reference is **four independent copies of the published `L4_LBC_cACint209_1` cell**, with 1,808 sections and 5,352 computational segments in total. This is not a measurement of four unspecified SWC morphologies or a connected network. The earlier timed simulation used one cell; the downloaded network wrapper selects three.

## Measured simulation memory

The optimized single-thread, fixed-step process used **77.6 MiB of resident RAM** after initializing and running four cells for 25 ms. This includes Python, NumPy, NEURON, mechanisms, cells, and preallocated recording buffers large enough for three seconds of soma recordings. Synapses, plotting, and full-compartment recording were disabled.

| Stage | Resident RAM |
| --- | ---: |
| Runtime loaded | 44.1 MiB |
| Four cells constructed | 73.3 MiB |
| Three-second soma buffers allocated | 74.4 MiB |
| After 25 ms probe | 77.6 MiB |

The soma buffers themselves need only **1.14 MiB**: four voltage vectors and one shared time vector, sampled every 0.1 ms. A practical allowance of **128–256 MiB** gives headroom for this headless configuration. Physical resident memory is the relevant measurement here, rather than the process's larger reserved virtual address space.

The full three-second four-cell run was not executed for this measurement. Core model storage is independent of simulation duration; recording buffers were preallocated for the target duration. Detailed measurements are in `outputs/memory/published_4_cells.json`.

## Full morphology recording

Recording every segment for the supplied 3D animator costs much more. Calculated array payloads for four cells are:

| Duration | Recording interval | Raw NEURON recordings | Animator conversion peak payload |
| --- | ---: | ---: | ---: |
| 25 ms | 0.1 ms | 10.3 MiB | 20.5 MiB |
| 25 ms | 0.025 ms | 40.9 MiB | 81.8 MiB |
| 3,000 ms | 0.1 ms | 1.20 GiB | 2.39 GiB |
| 3,000 ms | 0.025 ms | 4.79 GiB | 9.57 GiB |

These large-buffer figures are calculations, not measured allocations. Add the runtime's resident memory, vector spare capacity, geometry, plotting, and object overhead. The conversion peak comes from the original animator retaining double-precision NEURON vectors while creating per-vector float32 copies and then a stacked float32 matrix.

Storage scales with total segment count and recorded samples, rather than neuron count alone:

```text
double-precision voltage bytes ≈ segments × samples × 8
samples ≈ duration_ms / recording_interval_ms + 1
```

## Retained optimizations and reproduction

The normal published-example launcher now removes only static zero-conductance `Nap_Et2` and `K_Tst` instances by default, retaining fixed timestep `0.025 ms` and cache-efficient layout. A complete rerun still matched all 30,000 original soma samples bit for bit. Use `--no-prune-zero` to retain the unoptimized channel instances. Adaptive CVode remains an explicit benchmark option.

```sh
.venv/bin/python tools/run_published_example.py --step 1
.venv/bin/python tools/memory_published.py --cells 4 --duration 3000
```
