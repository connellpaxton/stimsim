# External experiment artifacts

Reviewed November 2, 2026. Historical recordings, videos and run manifests are outside this Git repository. Report links point here instead of implying that a fresh checkout contains `outputs/` artifacts. The listed roots are machine-specific; access requires the original workstation/NAS. Each listed file was confirmed present during this review. Files were not copied or modified.

For current checkout validation, see [validation status](VALIDATION.md). For published download URLs, see [data provenance](DATA_PROVENANCE.md).

## Full Minnie v1300

Artifact root: `/mnt/nas-storage/NeuronSim/experiments/minnie65_v1300_full_mpi96_FzV1DzbT`. Paths below are relative to this root.

- `EXPERIMENT_COMPLETE.json`
- `construction.json`
- `results/aggregate_memory.json`
- `results/long_protocol/summary.json`
- `results/network_response.png`
- `results/population_rates.png`
- `results/validation.json`
- `run_manifest.json`
- `run_mpi_full_trial.py`
- `video/first_500ms_population.html`
- `video/first_500ms_population.mp4`

## Proofread v1300 trial

Artifact root: `/mnt/nas-storage/NeuronSim/experiments/minnie65_v1300_proofread_1711_uwUH4zIS`. Paths below are relative to this root.

- `results/long_protocol/recordings.npz`
- `results/long_protocol/summary.json`
- `results/model_parameter_comparison.json`
- `results/network_response.png`
- `results/process_resources.txt`
- `results/validation.json`

## Full-cohort preparation

Artifact root: `/home/shisoft/Dropbox/Code/NeuronSim/outputs/nas_full_trial`. Paths below are relative to this root.

- `prepared/data/summary.json`
- `run_full_trial.py`
- `run_manifest.json`

## MPI reference comparison

Artifact root: `/home/shisoft/Dropbox/Code/NeuronSim/outputs/mpi_validation`. Paths below are relative to this root.

- `proofread96/results/reference_comparison.json`

## Column regression

Artifact root: `/home/shisoft/Dropbox/Code/NeuronSim/outputs/NeuronSim-minnie-column.WCVsAeuu`. Paths below are relative to this root.

- `outputs/cleanup_full100/refactor_regression.json`
- `outputs/cleanup_full100/validation.json`

## L23 Genoa

Artifact root: `/home/shisoft/Dropbox/Code/NeuronSim/outputs/NeuronSim-l23-current.9iguVDCM`. Paths below are relative to this root.

- `outputs/parallel32_1s/activity_analysis.json`
- `outputs/parallel32_1s/network_response.png`
- `outputs/parallel32_1s/recordings.npz`
- `outputs/parallel32_1s/summary.json`
- `outputs/validation.json`

## Three-neuron comparison

Artifact root: `/home/shisoft/Dropbox/Code/NeuronSim/outputs/wrapper_three_neuron`. Paths below are relative to this root.

- `README.md`
- `comparison.json`
- `comparison.png`
- `selected_contacts.csv`

## Three-neuron one-second

Artifact root: `/home/shisoft/Dropbox/Code/NeuronSim/outputs/wrapper_three_neuron_1s`. Paths below are relative to this root.

- `README.md`
- `timing_comparison.json`

## Original workspace mirrors

The original workstation output tree is `/home/shisoft/Dropbox/Code/NeuronSim/outputs/`. For example, `mpi_full_trial/` contains selected full-run summaries, figures and the video; the NAS full-run root above additionally contains the large merged recordings and per-rank output. Historical `outputs/...` paths in prose refer to that original workspace unless explicitly labeled as new stimsim checks.

Genoa workspace paths in individual reports identify the original execution environments. They are historical provenance, not evidence of an active job. No background simulation is being monitored by these documents.
