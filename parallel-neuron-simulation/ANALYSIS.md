# Analysis of the neuron_10-26-2026 archive

The archive contains real neuron reconstructions and contact coordinates plus a prototype NEURON simulation. Its morphology and connectivity data can be explored now. Python runtime dependencies have been resolved; reproducing the original active simulation still requires the author's custom mechanisms.

**Published-reference update:** the linked [ModelDB 241165 source was inspected and its bundled Python example executed successfully](MODELDB_241165_REVIEW.md). It implements regional behavior through shared mechanisms with per-section parameter assignments. This provides a documented route to adapting the downloaded script without finding every renamed MOD file. The downloaded three-cell network itself has not yet been ported or run successfully.

## Inventory and data quality

| Layer folder | SWC morphologies | Pass tree and positive-radius checks | Contact-table records |
| --- | ---: | ---: | ---: |
| L23 | 349 | 313 | 42,559 |
| L4 | 266 | 243 | 25,446 |
| L5 | 147 | 130 | 11,530 |
| L6 | 335 | 312 | 11,971 |
| Total | 1,097 | 998 | 91,506 |

- The morphologies contain 4,008,435 SWC nodes. All have one root, unique node IDs, valid parents, finite coordinates, and connected tree topology; no zero-length parent edges were found.
- There are **114 zero-radius nodes in 99 morphologies**, with no negative radii. Review the source reconstruction or choose and document a repair policy before using these cells in simulation. This analysis does not alter them.
- Contact tables contain 91,506 distinct rows, with no exact row duplication across tables. Both endpoint morphologies are supplied for 24,842 records; both pass geometry checks for **19,929** records. These checks do not establish biological completeness or validate membrane parameters.
- Table names do not imply that both cells belong to that layer. For example, `L5_all.csv` contains 3,732 contacts with both morphologies supplied, but only 753 have both endpoints in `SWC_L5`.
- Neuron root IDs are retained as strings to avoid loss of precision in identifiers of this size.

## What the supplied code does

`cell.py` imports SWC morphology, inserts compartment-specific active channels and extracellular mechanisms, chooses segment counts, and maps contact coordinates to the nearest computational segment. `SynapseConstructor` creates voltage-triggered `NetCon` connections to `Exp2Syn` synapses. `plotting3D.py` records segment voltage and animates it on the morphology.

The wrapper selects three layer-5 cells and connects them in a directed cycle:

```text
864691135640672571 -> 864691135993429185
864691135993429185 -> 864691136025333561
864691136025333561 -> 864691135640672571
```

Each edge has one contact in `L5_all.csv`. Under the code's nm-to-µm conversion, sample contact endpoints are 0.15–1.47 µm from their nearest SWC nodes. This supports coordinate alignment for the sample; it is not a whole-archive synapse-placement validation.

The wrapper injects 1 nA for 1 ms into the first soma, then runs for 25 ms at `dt=0.025 ms`. Its default synapses are all excitatory (`e=0 mV`, weight `0.01 µS`); the tables do not supply inhibitory/excitatory labels or calibrated synaptic weights.

## Original-run blockers and code findings

The initial attempt stopped at `ModuleNotFoundError: pandas`. An isolated Python 3.11 `.venv` now contains NEURON 9.0.1, NumPy 2.4.6, and Matplotlib 3.11.2, installed from existing local package caches because shell network access could not resolve PyPI. `pip check` passed and NEURON ran a built-in-channel runtime check for 3 ms (121 samples).

A working copy of the supplied scripts is in `model/`. It fixes morphology paths, disables animation by default, and replaces pandas CSV reading with `csv.DictReader`, preserving ID strings, numeric coordinates, and the original contact-selection logic. The wrapper's simulation statements and all `Cell` methods are unchanged, verified by comparing their Python syntax trees. A regression test confirms connection direction and parameters while preserving large IDs. The download's files remain unchanged; hashes are in `model/source-manifest.json`.

**The supplied wrapper working copy was then executed directly.** It got through imports and loaded the first morphology, then stopped at `sec.insert('xtra')` with `ValueError: argument not a density mechanism name`. This is a missing model-definition error, not a Python import error. The full traceback is in `outputs/original_attempt.log`; no original active-model voltage output has been produced.

`Cell` additionally references **28 custom density mechanisms**, including `xtra`, multiple calcium and sodium channels, and soma/apical variants. No `.mod` sources or compiled mechanism libraries are included. Installing NEURON alone cannot supply those custom definitions: NEURON requires the corresponding [compiled NMODL mechanisms](https://neuron.yale.edu/neuron/docs/using-nmodl-files). Obtain the matching files and intended NEURON version from the author rather than substituting similarly named models.

A subsequent source search located an `xtra.mod` with the expected `es`, `stim`, `ex`, and coordinate interface in [ModelDB 241165](https://github.com/ModelDBRepository/241165/blob/93e01d4a927ab0cc5ca137329b86925777e63f02/mechanisms/xtra.mod). It is saved unchanged in `mechanisms/`, compiled with NEURON 9.0.1, and passed an interface and voltage-pointer test. Exact provenance relative to the downloaded model is unconfirmed. With it loaded, the supplied wrapper advances past `xtra` and fails at `sec.insert('Ca')`; 27 other referenced mechanisms remain missing. This has not produced a completed active-model simulation.

The next source search recovered published candidates matching 11 of those 27 mechanism names, plus a related `NaTs2_t` base channel, from the pinned Aberra TMS simulation source. All 12 compile and insert successfully, with finite voltages in an isolated 0.1 ms check. All three passive reversal potentials also match a published cACint209 regional parameter profile. The remaining ten `_soma`, five `_apic`, and one `_nor` definitions were not located. This is evidence of shared model lineage, not proof of original parameter equivalence. The supplied wrapper now stops at `Ca_soma`; the complete 27-name results and caveats are in [mechanisms/SOURCE_SEARCH.md](mechanisms/SOURCE_SEARCH.md).

Other findings:

- The wrapper uses bare SWC filenames although files are inside `SWC_L5/`. The added launcher resolves full paths without editing the download.
- `SynapseConstructor.update_params(pre_weight=...)` raises `TypeError` because it compares optional time constants before resolving their existing values. This was reproduced in isolation. Initial construction supplies both constants and avoids this failure.
- Synapse insertion scans all rows and repeatedly rebuilds cell-coordinate arrays. Cache a cell-ID lookup and spatial index before scaling to thousands of cells.
- Channel parameters, temperatures, biological provenance, and validation targets are not documented sufficiently to treat output as a calibrated experiment.

## What actually ran

A new, explicitly **passive** cable demonstration ran the three original morphologies independently: 5,567, 6,811, and 6,451 compartments. Each cell received its own 1 nA, 1 ms soma pulse beginning at 1 ms, over a 25 ms simulation with `dt=0.025 ms`.

| Cell ID | Peak soma voltage (mV) | Peak farthest-node voltage (mV) |
| --- | ---: | ---: |
| 864691135640672571 | -36.84 | -64.74 |
| 864691135993429185 | -41.59 | -64.37 |
| 864691136025333561 | -38.33 | -64.67 |

This demonstrates current integration, passive voltage spread, and relaxation in the supplied geometry. It does not establish a spike threshold or reproduce network activity. Each non-root SWC node represents its parent edge as a cylindrical compartment; the root is a sphere. The model uses backward Euler and a tree matrix solver, `Ra=100 Ω·cm`, `cm=1 µF/cm²`, uniform leak `6.3e-5 S/cm²`, and rest `-65 mV`. It retains recorded radii and floors edges shorter than 0.01 µm; none of the sample edges needed that floor. The original model's channel physics and compartment discretization differ.

Five tests passed: the tree solver against a dense solve, a soma against the analytic RC response with time-step refinement, an unstimulated branched tree at rest, rejection of malformed trees, and preservation of the original model's CSV-driven connection direction and parameters. Passive output data are finite. The explorer's JavaScript passed syntax checking; browser visual verification was blocked by the local-file URL policy.

## Useful next steps

1. Explore reconstructed cell shapes and contact placement using `outputs/explorer.html`; compare passive responses using the saved CSV/SVG results.
2. Recover the author's exact mechanism sources and configuration, fix the optional-parameter update bug, and run the original three-cell cycle with the added launcher.
3. Resolve zero-radius geometry and audit all contact-to-morphology distances. Add provenance, coordinate units, cell identities, and reconstruction completeness metadata.
4. Once the model is validated, compare connected versus disconnected responses, sweep current amplitude and synaptic weight, and expand to selected within-layer or cross-layer subnetworks.

See [README.md](README.md) for commands. Full measurements and affected filenames are in `outputs/analysis.json` and `outputs/morphologies.csv`.
