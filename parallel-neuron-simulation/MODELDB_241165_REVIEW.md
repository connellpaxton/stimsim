# Review of ModelDB 241165

The linked [ModelDB archive](https://modeldb.science/241165) is the published code for [Aberra, Peterchev and Grill (2018), “Biophysically realistic neuron models for simulation of cortical stimulation”](https://pubmed.ncbi.nlm.nih.gov/30127100/). It provides a complete published reference for the underlying mechanisms and regional configuration. The extra names used by the downloaded Python script are not required by the published implementation.

## How it resolves the mechanism question

The archive has 16 unique MOD filenames, including the shared ion-channel/calcium mechanisms, two synapse models, and `xtra`. Inspecting all 366 MOD paths, including cell-specific copies, found no separate `_soma`, `_apic`, or `_nor` MOD files.

Instead, the HOC biophysics files insert shared mechanisms into section lists and assign parameters separately for each region. Examples from the published `L4_LBC_cACint209_1` profile:

| Downloaded script name | Published counterpart | Published regional setting |
| --- | --- | --- |
| `Ca_soma` | `Ca` on somatic sections | `gCabar=0.000028 S/cm²` |
| `NaTs2_t_soma` | `NaTs2_t` on somatic sections | `gNaTs2_tbar=0.395881 S/cm²` |
| `NaTs2_t_apic` | `NaTs2_t` on apical sections | `gNaTs2_tbar=0.001373 S/cm²` |
| `CaDynamics_E2_soma` | `CaDynamics_E2` on somatic sections | `gamma=0.000814`, `decay=967.678789 ms` |
| `NaTa_t_nor` | `NaTa_t` on node sections | The published myelin routine doubles the initial-segment sodium conductance |

The last row identifies a documented nodal implementation; equivalence to the download's unknown `_nor` file is inferred, not established. The paper's [myelination methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC6239949/) describe this doubled nodal sodium conductance.

The published L4 basket-cell profile also exactly matches the downloaded script's three passive reversal potentials: axon `-64.601696 mV`, soma `-69.781406 mV`, and basal/apical dendrites `-63.118492 mV`. This supports shared model lineage. It does not prove that the missing renamed files have identical channel kinetics or defaults, or that this profile is calibrated for every supplied SWC morphology.

The practical repair is to use published mechanisms with explicit regional assignments, including spatially varying parameters where specified. Treating 16 renamed files as the only route to a runnable model was too restrictive. The public reference supplies a documented implementation to adapt; an exact reproduction of the author's missing modifications remains unverified.

## Verified run

A representative subset of the official archive is saved in `vendor/aberra2018`, pinned to commit `93e01d4a927ab0cc5ca137329b86925777e63f02`. All 29 downloaded files match their published Git blob hashes. No upstream source files were modified. All 16 shared mechanisms compiled with NEURON 9.0.1.

The bundled `cells/L4_LBC_cACint209_1/run.py` was executed through its own `init_simulation()` and `run_step(1, plot_traces=False)` functions:

- Duration: 3,000 ms; integration step: 0.025 ms; temperature: 34°C.
- Holding current: -0.031738 nA; step current: 0.0630936 nA from 700 to 2,700 ms.
- Output: 30,000 finite soma-voltage samples at 0.1 ms intervals.
- Response: subthreshold; no upward soma-voltage crossings of 0 mV.

This is the archive's bundled Blue Brain baseline example. It does not execute the adult scaling/myelination pipeline or reproduce a paper figure. It also does not run the downloaded three-cell SWC network.

Results are in `outputs/published_model/`: the original DAT recording, a CSV export, a PNG plot, and `summary_step1.json`.

## Reproduce

```sh
(cd vendor/aberra2018/mechanisms && ../../../.venv/bin/nrnivmodl)
.venv/bin/python tools/run_published_example.py --step 1
```

The launcher verifies source checksums, loads the official mechanism library, invokes the unchanged upstream Python functions, and exports the results. Upstream attribution and licensing are preserved in the vendored files and manifest. The working copy in `model/` still requires an explicit port to the published regional configuration before its three-cell simulation can run.

The launcher now wraps cell creation to apply the validated static zero-conductance optimization by default. The source files remain unchanged. Add `--no-prune-zero` to reproduce the unoptimized setup. A full optimized step-1 rerun produced the same recorded soma trace exactly; the integration method remains fixed-step.
