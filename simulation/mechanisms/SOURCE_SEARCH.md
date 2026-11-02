# Mechanism source search

**Published-reference follow-up:** [ModelDB 241165 was checked directly](../MODELDB_241165_REVIEW.md), its mechanisms compiled, and its bundled Python example ran successfully. Its regional behavior is implemented using shared mechanisms plus section-specific parameter assignments; separate renamed MOD files are not required. The inventory below remains a record of exact-name matches, not a count of independent missing channel equations.

Of the 27 channel mechanisms missing after retrieving `xtra`, **11 have published sources with matching mechanism names**. Exact definitions for **16 variants were not located. All retrieved files remain candidates relative to the download: matching names do not establish identical author-specific parameters.

## Retrieved sources

Twelve unchanged channel files were retrieved from [Aman Aberra’s published TMS simulation repository](https://github.com/Aman-A/TMSsim_Aberra2019/tree/8ec575a1f1042abc64ed05d74e5f7e6672375d5d/nrn/mechanisms). Eleven match names required by `cell.py`; the extra `NaTs2_t.mod` is a related base channel for two unresolved variants. Sources, Git blob hashes, and SHA-256 checksums are in `reference/download-manifest.json`. Attribution and upstream licensing notes are preserved in `reference/UPSTREAM-LICENSE`.

All 12 retrieved channels compiled with NEURON 9.0.1, inserted into a section, and produced finite voltages in a 0.1 ms isolated runtime check. The existing `xtra` mechanism also remains in the compiled library. These are software checks, not a calibration against the downloaded model.

## Every requested name

| Required name | Search outcome | Published source or related base |
| --- | --- | --- |
| `Ca` | Matching SUFFIX found | [Ca.mod](https://github.com/Aman-A/TMSsim_Aberra2019/blob/8ec575a1f1042abc64ed05d74e5f7e6672375d5d/nrn/mechanisms/Ca.mod) |
| `CaDynamics_E2` | Matching SUFFIX found | [CaDynamics_E2.mod](https://github.com/Aman-A/TMSsim_Aberra2019/blob/8ec575a1f1042abc64ed05d74e5f7e6672375d5d/nrn/mechanisms/CaDynamics_E2.mod) |
| `CaDynamics_E2_soma` | Exact variant not located | [CaDynamics_E2.mod](https://github.com/Aman-A/TMSsim_Aberra2019/blob/8ec575a1f1042abc64ed05d74e5f7e6672375d5d/nrn/mechanisms/CaDynamics_E2.mod) (related base only) |
| `Ca_LVAst` | Matching SUFFIX found | [Ca_LVAst.mod](https://github.com/Aman-A/TMSsim_Aberra2019/blob/8ec575a1f1042abc64ed05d74e5f7e6672375d5d/nrn/mechanisms/Ca_LVAst.mod) |
| `Ca_LVAst_soma` | Exact variant not located | [Ca_LVAst.mod](https://github.com/Aman-A/TMSsim_Aberra2019/blob/8ec575a1f1042abc64ed05d74e5f7e6672375d5d/nrn/mechanisms/Ca_LVAst.mod) (related base only) |
| `Ca_soma` | Exact variant not located | [Ca.mod](https://github.com/Aman-A/TMSsim_Aberra2019/blob/8ec575a1f1042abc64ed05d74e5f7e6672375d5d/nrn/mechanisms/Ca.mod) (related base only) |
| `Ih` | Matching SUFFIX found | [Ih.mod](https://github.com/Aman-A/TMSsim_Aberra2019/blob/8ec575a1f1042abc64ed05d74e5f7e6672375d5d/nrn/mechanisms/Ih.mod) |
| `Im` | Matching SUFFIX found | [Im.mod](https://github.com/Aman-A/TMSsim_Aberra2019/blob/8ec575a1f1042abc64ed05d74e5f7e6672375d5d/nrn/mechanisms/Im.mod) |
| `Im_apic` | Exact variant not located | [Im.mod](https://github.com/Aman-A/TMSsim_Aberra2019/blob/8ec575a1f1042abc64ed05d74e5f7e6672375d5d/nrn/mechanisms/Im.mod) (related base only) |
| `Im_soma` | Exact variant not located | [Im.mod](https://github.com/Aman-A/TMSsim_Aberra2019/blob/8ec575a1f1042abc64ed05d74e5f7e6672375d5d/nrn/mechanisms/Im.mod) (related base only) |
| `K_Pst` | Matching SUFFIX found | [K_Pst.mod](https://github.com/Aman-A/TMSsim_Aberra2019/blob/8ec575a1f1042abc64ed05d74e5f7e6672375d5d/nrn/mechanisms/K_Pst.mod) |
| `K_Pst_apic` | Exact variant not located | [K_Pst.mod](https://github.com/Aman-A/TMSsim_Aberra2019/blob/8ec575a1f1042abc64ed05d74e5f7e6672375d5d/nrn/mechanisms/K_Pst.mod) (related base only) |
| `K_Pst_soma` | Exact variant not located | [K_Pst.mod](https://github.com/Aman-A/TMSsim_Aberra2019/blob/8ec575a1f1042abc64ed05d74e5f7e6672375d5d/nrn/mechanisms/K_Pst.mod) (related base only) |
| `K_Tst` | Matching SUFFIX found | [K_Tst.mod](https://github.com/Aman-A/TMSsim_Aberra2019/blob/8ec575a1f1042abc64ed05d74e5f7e6672375d5d/nrn/mechanisms/K_Tst.mod) |
| `K_Tst_apic` | Exact variant not located | [K_Tst.mod](https://github.com/Aman-A/TMSsim_Aberra2019/blob/8ec575a1f1042abc64ed05d74e5f7e6672375d5d/nrn/mechanisms/K_Tst.mod) (related base only) |
| `K_Tst_soma` | Exact variant not located | [K_Tst.mod](https://github.com/Aman-A/TMSsim_Aberra2019/blob/8ec575a1f1042abc64ed05d74e5f7e6672375d5d/nrn/mechanisms/K_Tst.mod) (related base only) |
| `NaTa_t` | Matching SUFFIX found | [NaTa_t.mod](https://github.com/Aman-A/TMSsim_Aberra2019/blob/8ec575a1f1042abc64ed05d74e5f7e6672375d5d/nrn/mechanisms/NaTa_t.mod) |
| `NaTa_t_nor` | Exact variant not located | [NaTa_t.mod](https://github.com/Aman-A/TMSsim_Aberra2019/blob/8ec575a1f1042abc64ed05d74e5f7e6672375d5d/nrn/mechanisms/NaTa_t.mod) (related base only) |
| `NaTs2_t_apic` | Exact variant not located | [NaTs2_t.mod](https://github.com/Aman-A/TMSsim_Aberra2019/blob/8ec575a1f1042abc64ed05d74e5f7e6672375d5d/nrn/mechanisms/NaTs2_t.mod) (related base only) |
| `NaTs2_t_soma` | Exact variant not located | [NaTs2_t.mod](https://github.com/Aman-A/TMSsim_Aberra2019/blob/8ec575a1f1042abc64ed05d74e5f7e6672375d5d/nrn/mechanisms/NaTs2_t.mod) (related base only) |
| `Nap_Et2` | Matching SUFFIX found | [Nap_Et2.mod](https://github.com/Aman-A/TMSsim_Aberra2019/blob/8ec575a1f1042abc64ed05d74e5f7e6672375d5d/nrn/mechanisms/Nap_Et2.mod) |
| `Nap_Et2_soma` | Exact variant not located | [Nap_Et2.mod](https://github.com/Aman-A/TMSsim_Aberra2019/blob/8ec575a1f1042abc64ed05d74e5f7e6672375d5d/nrn/mechanisms/Nap_Et2.mod) (related base only) |
| `SK_E2` | Matching SUFFIX found | [SK_E2.mod](https://github.com/Aman-A/TMSsim_Aberra2019/blob/8ec575a1f1042abc64ed05d74e5f7e6672375d5d/nrn/mechanisms/SK_E2.mod) |
| `SK_E2_soma` | Exact variant not located | [SK_E2.mod](https://github.com/Aman-A/TMSsim_Aberra2019/blob/8ec575a1f1042abc64ed05d74e5f7e6672375d5d/nrn/mechanisms/SK_E2.mod) (related base only) |
| `SKv3_1` | Matching SUFFIX found | [SKv3_1.mod](https://github.com/Aman-A/TMSsim_Aberra2019/blob/8ec575a1f1042abc64ed05d74e5f7e6672375d5d/nrn/mechanisms/SKv3_1.mod) |
| `SKv3_1_apic` | Exact variant not located | [SKv3_1.mod](https://github.com/Aman-A/TMSsim_Aberra2019/blob/8ec575a1f1042abc64ed05d74e5f7e6672375d5d/nrn/mechanisms/SKv3_1.mod) (related base only) |
| `SKv3_1_soma` | Exact variant not located | [SKv3_1.mod](https://github.com/Aman-A/TMSsim_Aberra2019/blob/8ec575a1f1042abc64ed05d74e5f7e6672375d5d/nrn/mechanisms/SKv3_1.mod) (related base only) |

## Evidence of model lineage

The downloaded code’s passive reversal potentials match all three values in the published [cACint209 regional biophysics file](https://github.com/Aman-A/TMSsim_Aberra2019/blob/8ec575a1f1042abc64ed05d74e5f7e6672375d5d/nrn/cells/L4_LBC_cACint209_1/biophysics.hoc): axon `-64.601696 mV`, soma `-69.781406 mV`, and basal/apical dendrites `-63.118492 mV`. A copy is saved in `reference/cACint209_biophysics.hoc` for review, without executing it.

This suggests a shared model family. It does not prove that the missing regional mechanism files use identical kinetics or defaults. In the published model, conductances and calcium-handling parameters are assigned separately for each region. For example, `NaTa_t.mod` defaults to `gNaTa_tbar=0.00001 S/cm²`, while this regional profile assigns `3.993125 S/cm²` to axons. The downloaded Python code does not assign active-channel conductances explicitly; its missing mechanism files may encode different defaults.

## Unresolved definitions

The 16 unresolved names comprise ten `_soma` variants, five `_apic` variants, and `NaTa_t_nor`. Exact SUFFIX searches were completed for each; query results are retained in `reference/variant-search.json`. Additional name searches and web searches found no exact definition. Public search coverage is incomplete, so this is not proof that the files do not exist online.

A lowercase `ca_soma` search hit was rejected: it is an unrelated calcium-accumulation mechanism, not the case-sensitive `Ca_soma` calcium-current mechanism requested here. Hits such as `gSK_E2bar_SK_E2_soma` can describe a parameter on a soma section without defining a mechanism named `SK_E2_soma`.

No regional mechanisms were manufactured by renaming base files. Reconstructing a runnable model using the published regional profile is possible as a separately documented adaptation, but would require explicit choices about regional conductances, calcium handling, ion reversal potentials, temperature, and the meaning of `_nor`.

## Current run status

The supplied wrapper now gets past the retrieved axon mechanisms and stops at `sec.insert('Ca_soma')`. The traceback is saved as `outputs/original_attempt_with_published_channels.log`. The complete original active-model simulation has not run.

```sh
.venv/bin/python tools/run_original.py --mechanisms mechanisms --check
```
