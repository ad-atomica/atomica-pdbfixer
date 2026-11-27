#!/usr/bin/env bash
# Reproduce the deposited-structure comparison in docs/2026-11-26-rebuilt-chirality-restraint.md.
#
# Usage: devtools/verify_rebuilt_chirality.sh BENCH_DIR [quick]
#
# BENCH_DIR receives a virtual environment, the deposited inputs, copies of the release and
# parent revisions, and the measured CSVs.  Rerunning reuses the environment, inputs, and
# revision copies.  The fix revision is this checkout.  "quick" measures only ubiquitin with
# seed 1 on the parent and fix revisions.

set -euo pipefail

if [ $# -lt 1 ] || [ $# -gt 2 ] || { [ $# -eq 2 ] && [ "$2" != quick ]; }; then
    echo "usage: $0 BENCH_DIR [quick]" >&2
    exit 2
fi
mode=${2:-full}
mkdir -p "$1"
bench=$(cd "$1" && pwd)
repo=$(cd "$(dirname "$0")/.." && pwd)
cd "$repo"
python="$bench/venv/bin/python"
runner=devtools/verify_rebuilt_chirality.py

# Both commands are idempotent, so rerunning them also repairs an interrupted setup.

mise exec -- uv venv --python 3.12.16 --managed-python --allow-existing "$bench/venv"
mise exec -- uv pip install --python "$bench/venv" 'numpy==2.5.5' 'openmm==8.6.1'
"$python" "$runner" fetch "$bench/inputs"

# Extract a commit's pdbfixer package into BENCH_DIR/REVISION unless an earlier run did.  The
# package is moved into place only once complete, so an interrupted run is not reused.

extract() {
    local revision=$1 commit=$2
    if [ ! -d "$bench/$revision" ]; then
        rm -rf "$bench/$revision.partial"
        mkdir "$bench/$revision.partial"
        git archive "$commit" pdbfixer | tar -x -C "$bench/$revision.partial"
        mv "$bench/$revision.partial" "$bench/$revision"
    fi
}
extract release v1.12
extract parent bdf98bc

# Measure one revision, importing pdbfixer from its extracted copy or, for fix, this checkout.

measure() {
    local revision=$1 source="$bench/$1"
    shift
    if [ "$revision" = fix ]; then
        source=$repo
    fi
    OPENMM_CPU_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONPATH="$source" \
        "$python" "$runner" measure --revision "$revision" "$bench/inputs" "$@"
}

if [ "$mode" = quick ]; then
    measure parent 1 --pdb 1UBQ
    measure fix 1 --pdb 1UBQ | tail -n +2
    exit
fi

# The runs are independent single-threaded processes with fixed seeds, so run them in parallel.
# They are waited for in launch order, so a failed run stops the script once every run launched
# before it has finished.  set -m gives each run its own process group, so the trap then stops
# the unfinished runs, with their Python processes.  It does the same after Ctrl-C or SIGTERM.

revisions=(release parent fix)
seeds=(1 17 12345)
set -m
pids=()
trap 'for pid in ${pids[@]+"${pids[@]}"}; do kill -- "-$pid" 2> /dev/null || true; done' EXIT
for revision in "${revisions[@]}"; do
    for seed in "${seeds[@]}"; do
        measure "$revision" "$seed" > "$bench/$revision-$seed.csv" &
        pids+=($!)
    done
    measure "$revision" 1 --pdb 1UBQ --platform CPU > "$bench/$revision-cpu.csv" &
    pids+=($!)
done
for pid in "${pids[@]}"; do
    wait "$pid"
done

# Combine the runs in the committed order.

{
    head -n 1 "$bench/${revisions[0]}-${seeds[0]}.csv"
    for revision in "${revisions[@]}"; do
        for seed in "${seeds[@]}"; do tail -n +2 "$bench/$revision-$seed.csv"; done
    done
    for revision in "${revisions[@]}"; do tail -n +2 "$bench/$revision-cpu.csv"; done
} > "$bench/rebuilt-chirality.csv"

"$python" "$runner" summarize "$bench/rebuilt-chirality.csv"
"$python" "$runner" compare docs/data/rebuilt-chirality.csv "$bench/rebuilt-chirality.csv"
