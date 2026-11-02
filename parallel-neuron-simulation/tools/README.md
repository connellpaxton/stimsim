# Tool entry points

The current simulation implementation lives in [`../neuronsim/`](../neuronsim/). Start with [the code walkthrough](../docs/CODE_WALKTHROUGH.md).

For new runs, select one of the four [dataset entry files](../docs/DATASET_ENTRY_POINTS.md) in `parallel-neuron-simulation/`. `tools/run_l23.py` remains a legacy compatibility launcher; the new `../run_l23.py` defaults to all valid L23 cells and the current electrical profile.

## Current Minnie / ModelDB 241165 workflow

| Tool | Purpose |
| --- | --- |
| `prepare_minnie.py` | Extract typed Minnie identities and observed within-column contact counts |
| `prepare_minnie_skeletons.py` | Fetch/validate per-neuron SWCs and preserve preprocessing provenance |
| `download_minnie_simulation_data.py` | Download the v1300 SWCs, wiring tables and cell metadata to NAS with source checksums and resumable ranges |
| `index_minnie_metadata.py` | Create queryable cell metadata with exact IDs, decoded coordinates and links to the SWC catalog |
| `prepare_nas_network.py` | Select and validate v1300 proofread E/I cells; stream observed internal contact counts into simulation inputs |
| `analyze_nas_trial.py` | Validate a proofread-cohort run against its zero-transmission control and exact replay; plot population activity |
| `compare_mpi_reference.py` | Check every merged MPI synapse placement, voltage sample and spike against an earlier run |
| `run_minnie_241165.py` | Compatible launcher for `python -m neuronsim` |
| `analyze_minnie_241165.py` | Compare coupled/control and E/I/all-excitatory conditions |
| `compare_minnie_threads.py` | Compare recordings and resources across native thread counts |
| `compare_wrapper_three.py` | Run the exact three-cell wrapper dataset; compare coordinate-matched and column-default connections |
| `run_l23.py --original-family` | Run the supplied L23 SWCs/contact coordinates with the current cACint209 profile and native threads |
| `validate_minnie_video.py` | Check the long recording against the earlier benchmark |
| `render_minnie_video.py` | Create MP4/player from recorded soma and segment voltages |
| `render_population_video.py` | Create a spatial soma-voltage/spike MP4 from existing recordings and node positions; supports parallel rendering without branch recordings |

`published_swc.py`, `branch_recording.py`, and `simulation_optimizations.py` retain earlier import names and forward to the shared package.

## Earlier experiments

| Tools | Configuration |
| --- | --- |
| `run_minnie_column.py`, `reconstructed_cell.py`, `export_m1_electrical_profiles.py`, `probe_minnie_branches.py`, `analyze_minnie_reconstructions.py` | Earlier Minnie geometry with transferred ModelDB 260015 M1 electrical profiles |
| `run_mouse_m1.py` | Published M1 network with generated connectivity |
| `run_l23.py`, `run_genoa.sh` | Downloaded L23 subset and earlier remote benchmark |
| `run_published_example.py`, `benchmark_published.py`, `memory_published.py` | Published single-cell baseline, optimization, and memory measurements |
| `run_original.py` | Original supplied wrapper; exact custom mechanisms remain incomplete |
| `analyze.py`, `passive_demo.py`, `explorer.html` | Archive inspection and passive cable demonstration |

`minnie_network.py` retains the earlier M1 template/receptor mappings; its generic graph helper is shared from `neuronsim.graph`. It does not configure the current ModelDB 241165 model.
