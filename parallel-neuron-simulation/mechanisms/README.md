# Retrieved mechanism sources

## Channel source search

See [SOURCE_SEARCH.md](SOURCE_SEARCH.md) for the full 27-name inventory. Eleven required channel names now have published source candidates, and `NaTs2_t.mod` was also retrieved as a related base channel. The other 16 exact regional/node variants were not located. The sources are pinned to an upstream commit and preserved unchanged; matching names do not verify the original author's parameter values.

All 12 retrieved channels compile and pass an isolated insertion/finite-voltage check. With the complete retrieved library loaded, the supplied wrapper now reaches the first missing soma variant, `Ca_soma`. No regional variants have been created by renaming base mechanisms. The full original simulation remains incomplete.

## xtra.mod

Retrieved from the [ModelDB 241165 source repository](https://github.com/ModelDBRepository/241165/blob/93e01d4a927ab0cc5ca137329b86925777e63f02/mechanisms/xtra.mod), pinned to commit `93e01d4a927ab0cc5ca137329b86925777e63f02`. The source credits Ted Carnevale and a 2018 modification by Aman Aberra. The file is preserved byte for byte; its Git blob SHA is `6fc62b8310d7c821cc602d55c4802f34429f9368`.

This version supplies the interface expected by the downloaded `cell.py`: `SUFFIX xtra`, per-segment `x`, `y`, `z` and `es`, global `stim`, and pointer `ex`. It sets the pointed-to voltage to `stim * es` during initialization and before voltage equations are assembled. The file contains additional `type` and `order` fields that the Python script does not use.

It compiled successfully with NEURON 9.0.1. An isolated section test verified insertion, both coordinate-access styles, the `ex` pointer connection to `e_extracellular`, amplitude scaling, and resetting the external voltage when the stimulus becomes zero. These checks establish interface compatibility, not that this is the original author's exact version. No ion-channel models were substituted.

Build and check from the repository root:

```sh
(cd mechanisms && ../.venv/bin/nrnivmodl)
.venv/bin/python -m unittest discover -s tests -v
.venv/bin/python tools/run_original.py --mechanisms mechanisms --check
```

The compiler reports that POINTER use is not thread safe; the validation uses a single-threaded simulation. The initial `xtra`-only attempt moved the wrapper to a missing `Ca` mechanism; the subsequent channel-source search above moves it to `Ca_soma`. See `outputs/xtra_validation.json`, `outputs/original_attempt_with_xtra.log`, and `outputs/original_attempt_with_published_channels.log`.
