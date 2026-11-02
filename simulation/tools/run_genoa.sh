#!/usr/bin/env bash
# Run from a terminal with SSH access. --prepare-only creates the payload locally.
set -euo pipefail
repo=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
host=${1:-shisoft@192.168.10.239}
data=${NEURONSIM_DATA:-/home/shisoft/Downloads/neuron_10-26-2026-20261027T002512Z-1-001/neuron_10-26-2026}
archive="$repo/outputs/genoa_payload.tar.gz"
mkdir -p "$repo/outputs"
python3 - "$repo" "$data" "$archive" <<'PY'
from pathlib import Path
import sys, tarfile
root, data, archive = map(Path, sys.argv[1:])
files = [root / 'requirements-neuron.lock.txt']
files += [root / 'tools' / name for name in ['run_l23.py', 'published_swc.py', 'simulation_optimizations.py', 'analyze.py']]
files += sorted((root / 'neuronsim').glob('*.py'))
published = root / 'vendor/aberra2018'
files += [published / 'L23_SOURCE_MANIFEST.json', published / 'SOURCE_MANIFEST.json']
files += sorted((published / 'mechanisms').glob('*.mod'))
files += sorted((published / 'cells/L23_PC_cADpyr229_1').glob('*'))
swcs = sorted((data / 'SWC_L23').glob('*.swc'))
if not swcs or not (data / 'L23_all.csv').is_file():
    raise SystemExit('L23 data not found; set NEURONSIM_DATA to its folder')
with tarfile.open(archive, 'w:gz') as package:
    for path in files:
        if path.is_file(): package.add(path, arcname=str(path.relative_to(root)))
    for path in swcs + [data / 'L23_all.csv']:
        package.add(path, arcname='data/' + str(path.relative_to(data)))
print(f'Prepared {archive} ({archive.stat().st_size / 2**20:.1f} MiB)')
PY
if [[ "$host" == '--prepare-only' ]]; then exit 0; fi
if [[ ! "$host" =~ ^[a-zA-Z0-9_.@:-]+$ || "$host" == -* ]]; then
    echo 'Invalid SSH destination' >&2; exit 2
fi
ssh_options=(-F /dev/null -o BatchMode=yes -o StrictHostKeyChecking=yes -o ConnectTimeout=10)
if [[ -f "$HOME/.ssh/config" ]]; then ssh_options[1]="$HOME/.ssh/config"; fi
remote_dir=$(ssh "${ssh_options[@]}" "$host" 'mktemp -d "$HOME/NeuronSim-genoa.XXXXXXXX"')
if [[ ! "$remote_dir" =~ ^/[a-zA-Z0-9_./-]+$ ]]; then
    echo 'Unexpected remote directory name; stopping' >&2; exit 2
fi
echo "Remote directory: $remote_dir"
ssh "${ssh_options[@]}" "$host" "tar xzf - -C '$remote_dir'" < "$archive"
ssh "${ssh_options[@]}" "$host" "bash -s -- '$remote_dir'" <<'REMOTE'
set -euo pipefail
cd -- "$1"
mkdir -p outputs
hostname > outputs/hardware.txt
lscpu >> outputs/hardware.txt
free -h >> outputs/hardware.txt
python_cmd=''
for candidate in python3.11 python3.12 python3.13 python3; do
    if command -v "$candidate" >/dev/null && "$candidate" -c 'import sys; sys.exit(not ((3,11) <= sys.version_info[:2] < (3,14)))'; then
        python_cmd=$candidate; break
    fi
done
if [[ -z "$python_cmd" ]]; then
    echo 'Python 3.11–3.13 with venv support is required' >&2; exit 1
fi
command -v g++ >/dev/null
command -v make >/dev/null
"$python_cmd" -m venv .venv
.venv/bin/python -m pip install -r requirements-neuron.lock.txt
.venv/bin/python -m pip check
(cd vendor/aberra2018/mechanisms && ../../../.venv/bin/nrnivmodl)
.venv/bin/python - <<'PY'
import csv, json, os, platform, subprocess, sys
from pathlib import Path
import numpy as np
cpus = len(os.sched_getaffinity(0))
rows, baseline, baseline_spikes = [], None, None
for threads in [n for n in [1, 2, 4, 8, 16, 32] if n <= cpus]:
    output = Path('outputs') / f'threads_{threads}'
    command = [sys.executable, 'tools/run_l23.py', '--data', 'data', '--max-cells', '0', '--amplitude', '3', '--duration', '25', '--threads', str(threads), '--output', str(output)]
    with Path(f'outputs/threads_{threads}.log').open('w') as log:
        subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
    summary = json.loads((output / 'summary.json').read_text())
    with np.load(output / 'recordings.npz') as recorded:
        current = {key: recorded[key].copy() for key in recorded.files}
    spikes = json.loads((output / 'soma_spikes.json').read_text())
    if baseline is None: baseline, baseline_spikes = current, spikes
    if not np.array_equal(baseline['cell_ids'], current['cell_ids']) or not np.array_equal(baseline['time_ms'], current['time_ms']):
        raise RuntimeError('Incompatible recordings across thread counts')
    delta = current['soma_mv'] - baseline['soma_mv']
    row = {key: summary[key] for key in ['threads', 'integration_seconds', 'initialization_seconds', 'solve_seconds', 'integration_cpu_seconds', 'resident_ram_bytes_before_plotting', 'soma_spike_count', 'estimated_delivered_synaptic_events']}
    row.update(max_voltage_difference_mv=float(abs(delta).max()), rms_voltage_difference_mv=float(np.sqrt(np.mean(delta**2))), spike_times_exact=spikes == baseline_spikes)
    row['integration_speedup'] = rows[0]['integration_seconds'] / row['integration_seconds'] if rows else 1.
    rows.append(row)
    Path('outputs/benchmark.json').write_text(json.dumps({'hostname': platform.node(), 'available_cpus': cpus, 'runs': rows, 'note': 'One run per count; timings include initialization. Numerical differences are reported, not assumed absent.'}, indent=2) + '\n')
    print(json.dumps(row), flush=True)
with Path('outputs/benchmark.csv').open('w', newline='') as stream:
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
    writer.writeheader(); writer.writerows(rows)
PY
REMOTE
results="$repo/outputs/$(basename -- "$remote_dir")"
mkdir -p "$results"
ssh "${ssh_options[@]}" "$host" "tar czf - -C '$remote_dir' outputs" | tar xzf - -C "$results"
echo "Results: $results/outputs/benchmark.csv"
echo "Remote files retained at $host:$remote_dir"
