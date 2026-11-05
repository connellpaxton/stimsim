# Single-thread performance review

The published `L4_LBC_cACint209_1` example has useful optimization opportunities, primarily in the compiled simulation rather than Python orchestration. Four sequential benchmarks used one NEURON thread, one process, the same 452 sections and 1,338 computational segments, and the same 3-second step-1 current-clamp protocol. Upstream sources remain unchanged.

## Measured results

| Configuration | Integration time | Speedup | Maximum recorded soma difference from fixed-step baseline |
| --- | ---: | ---: | ---: |
| Published fixed step, `dt=0.025 ms` | 52.32 s | 1.00× | 0 |
| Fixed step, remove static zero-conductance channels | 37.30 s | 1.40× | 0 mV, exactly identical samples |
| Adaptive CVode, `atol=rtol=1e-6` | 13.57 s | 3.86× | 0.005244 mV |
| Adaptive CVode, `atol=rtol=1e-8` | 16.54 s | 3.16× | 0.005342 mV |

These are single measurements per configuration, not statistical confidence intervals. Timings cover `h.run()`, including initialization and recording. They exclude imports, compilation, setup, and export. The measured fixed-step baseline exactly matches the earlier published-example recording.

Setup took 0.24–0.26 seconds. Array copying and binary export took about 0.0014 seconds. The baseline spent 52.32 of 52.60 measured seconds in `h.run()`. Native integration time was measured explicitly because direct HocObject calls are not visible to Python's cProfile hooks. The benchmark now wraps that call in a Python function for inclusive profiling.

## Most conservative improvement

The published configuration sets `Nap_Et2` conductance to zero in 1,337 segments and `K_Tst` conductance to zero in 1,001 segments. Their gating equations still execute. The `prune-zero` experiment removes these mechanisms only from sections where their conductance is exactly zero in every segment.

This saved about 15 seconds, or 29%, while preserving the timestep, spatial discretization, and all 30,000 recorded soma-voltage samples exactly. The soma's nonzero persistent sodium conductance remains present. This optimization is limited to the static protocol: it is inappropriate if later code activates those channels, records their gating states, or references them through pointers. The experiment does not claim equality of states that were removed.

## Largest measured improvement

The source explicitly disables adaptive integration in `run.py`. Enabling [NEURON's CVode integrator](https://nrn.readthedocs.io/en/8.2.7/hoc/simctrl/cvode.html) lets it choose timesteps according to error tolerances. Both adaptive experiments retain all original mechanisms and set the calcium absolute-tolerance scale to `1e-3`.

The two adaptive traces differ by at most **0.000155 mV**, with RMS difference **0.0000331 mV**. Their differences from the fixed-step trace are not certified error bounds: the fixed-step trace also has discretization error. The tested protocol is subthreshold, with no positive soma crossings. Spike timing, recruitment thresholds, and other cells need separate convergence checks before applying the same solver settings.

The bundled example disables synapses. The archive's `ProbAMPANMDA_EMS` and `ProbGABAAB_EMS` mechanisms report that they cannot be used with CVode. This result therefore does not establish compatibility for every published stimulation or synaptic configuration.

## Other findings

- Cache-efficient memory layout is already enabled, and the installed NEURON build is a release build. Enabling those again offers no improvement here.
- For repeated sweeps, reuse an instantiated cell and update current amplitudes, reinitializing all state for each run. This saves setup time, which is small relative to fixed-step integration in this example.
- NumPy array-copy and text-export cleanup can reduce allocations and output size, but cannot explain or materially eliminate the 52-second solver cost.
- Several MOD rate functions repeatedly evaluate exponentials and temperature factors. Hoisting invariant factors, rate tables, or CoreNEURON are further candidates, but were not benchmarked. Rate tables change numerical approximations and require validation.
- Separately, the downloaded network's `find_synapse_segment()` rebuilds segment-coordinate arrays for every contact. Cache per-cell coordinates and cell-ID lookups, and consider a spatial index for larger networks. Recording every segment for animation can also dominate memory. Neither path ran in this published-example benchmark.

## Reproduce

```sh
.venv/bin/python tools/benchmark_published.py --solver baseline
.venv/bin/python tools/benchmark_published.py --solver prune-zero
.venv/bin/python tools/benchmark_published.py --solver cvode
.venv/bin/python tools/benchmark_published.py --solver cvode-tight
```

Run cases sequentially to avoid CPU contention. JSON timings, voltage arrays, and profiling notes are saved in `outputs/benchmarks/`. The normal published-example runner now enables the validated zero-conductance removal by default while retaining fixed timesteps. Pass `--no-prune-zero` to retain the unoptimized channel instances. Adaptive integration remains an explicit benchmark option. The shared implementation is in `tools/simulation_optimizations.py`.
