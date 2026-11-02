# Minnie recording example

`branch_roots.json` selects the three real Minnie neurons used in the propagation video:

| Root ID in v943 | Role | Class |
| --- | --- | --- |
| `864691135395143157` | Injected source | 5P-IT, excitatory |
| `864691136444943747` | Connected recipient | 5P-PT, excitatory |
| `864691135561619681` | Connected recipient | BC, inhibitory |

The prepared v943 graph contains two contacts from the source to the excitatory recipient and seven to the inhibitory recipient. All three still receive input from the full network. The file selects recordings; it does not change connections or membrane properties.

Use it with `--branch-roots examples/minnie/branch_roots.json`. The default deterministic stimulus selection matches the saved video's 32 injected cells on this data snapshot.

The input graph and SWCs are downloaded/generated data under `outputs/`; they are not committed. See [the setup and code walkthrough](../../docs/CODE_WALKTHROUGH.md).
