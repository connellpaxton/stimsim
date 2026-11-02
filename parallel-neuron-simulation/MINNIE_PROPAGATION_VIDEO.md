# One-second Minnie propagation video

The completed experiment simulates **1,000 ms (one second)** of the current Minnie E/I column using published ModelDB 241165 parameters. This follows the revised request for one second rather than one minute.

- [Watch the MP4](outputs/NeuronSim-minnie-column.WCVsAeuu/outputs/aberra241165_ei96_1s/video/propagation.mp4)
- [Open the local video player](outputs/NeuronSim-minnie-column.WCVsAeuu/outputs/aberra241165_ei96_1s/video/propagation.html)
- [Animated preview of the first 90 ms](outputs/NeuronSim-minnie-column.WCVsAeuu/outputs/aberra241165_ei96_1s/video/propagation_preview.gif)

## Experiment and result

The model retains the [same Minnie anatomy, connectivity, and E/I configuration](MINNIE_241165_EI.md): 1,352 neurons, 124,341 contacts, and 525,498 electrical compartments. Thirty-two excitatory neurons receive **0.5 nA from 10–60 ms**, with no further injected current during the second. NEURON advances 40,000 fixed steps of 0.025 ms on 96 native threads; the quiet part of the recording is also simulated.

There were **938 spikes across 402 neurons**: 632 excitatory-cell spikes and 306 inhibitory-cell spikes. The first 100 ms contained 781 spikes and the next 100 ms contained 157; the final spike occurred at **134.725 ms**. The remaining interval shows voltage recovery after the initial stimulus-driven response. This demonstrates the behavior of the chosen model, not measured in-vivo activity.

Construction took **209.64 s**, solving took **174.14 s**, and resident memory after simulation was **4.78 GiB**, including the additional recording buffers. Video rendering took approximately 131 s locally. Timings are single trials and memory is a resident snapshot rather than a measured peak.

## Reading the video

The left panel shows all neurons at an oblique projection of their actual soma positions. Circles/triangles distinguish anatomical E/I classes, gold rings identify the injected cells, and white halos identify recorded spikes within the preceding 2 ms. Colors show recorded soma voltage.

The three detailed panels show actual imported cable centerlines colored by **individually recorded computational-segment voltage**, rather than coloring an entire tree with its soma voltage:

| Panel | Minnie root ID | Class | Recorded segments |
| --- | --- | --- | ---: |
| 1: injected source | `864691135395143157` | 5P-IT, excitatory | 844 |
| 2: connected recipient | `864691136444943747` | 5P-PT, excitatory | 468 |
| 3: connected recipient | `864691135561619681` | BC, inhibitory | 773 |

The source has two recorded contacts onto recipient 2 and seven onto recipient 3. These cells were selected from the earlier benchmark for their connectivity and early responses. Their voltages reflect inputs from the full network; the video does not isolate a single presynaptic cause for each response. No traveling arrows or uncomputed activity are added.

The MP4 contains **1,401 frames at 30 fps**, encoded as H.264/yuv420p at 1600×1000, for **46.7 seconds of playback**. The first 200 ms advances by 0.2 model ms per frame; the rest advances by 2 model ms per frame. The on-screen clock and playback scale identify this change. The GIF previews only the initial 90 model ms.

All soma voltages are saved every 0.1 ms as float64. The selected 2,085 branch-segment recordings are visualization copies saved every 0.1 ms as float32. Branch polylines preserve the imported 3D bends and segment boundaries; line widths are a display convention.

## Verification and limitations

The first 100 ms reproduces every soma voltage and spike time/neuron ID from the prior 96-thread benchmark exactly, with a maximum voltage difference of **0 mV**. All selected branch soma recordings match their population recordings after float32 conversion. All recorded voltages and final compartment voltages are finite. See [the validation report](outputs/NeuronSim-minnie-column.WCVsAeuu/outputs/aberra241165_ei96_1s/validation.json).

Twenty-four tests pass, including preservation of a curved cable's geometry, branch/soma recording alignment, regional densities, and inhibitory-event polarity. The full MP4 decodes without errors, and representative active/recovery frames were visually checked.

Minnie supplies the structural data and E/I labels. The same published L4 basket-cell profile is applied to all neurons; channel parameters are not fitted to individual Minnie cells. Synaptic strengths, kinetics, delays, and attachment positions remain assumptions. The video depicts the resulting simulation rather than a validated physiological reconstruction.

## Reproduce

The Genoa workspace is `/home/shisoft/NeuronSim-minnie-column.WCVsAeuu`. The selected roots are in `data/branch_video_roots.json`.

```sh
.venv/bin/python tools/run_minnie_241165.py --data data/column --morphologies data/reconstructions --output outputs/aberra241165_ei96_1s --threads 96 --inhibition --duration-ms 1000 --branch-roots data/branch_video_roots.json --stimulus-reference outputs/reconstructed32/summary.json
```

After copying recordings locally, use the repository environment and system FFmpeg:

```sh
.venv/bin/python tools/validate_minnie_video.py --workspace outputs/NeuronSim-minnie-column.WCVsAeuu
.venv/bin/python tools/render_minnie_video.py --workspace outputs/NeuronSim-minnie-column.WCVsAeuu
```

The renderer generates the MP4, sample frames, playback metadata, and an HTML player. Automatic opening through the browser tool was blocked because that tool disallows local `file:` URLs; the files can be opened directly by the user.
