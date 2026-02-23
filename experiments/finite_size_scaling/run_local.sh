#!/bin/bash
# Run finite-size scaling simulations locally.
#
# Usage:
#   ./run_local.sh                          # run all sizes sequentially
#   ./run_local.sh 250                      # run N=250 only
#   ./run_local.sh 250 500                  # run N=250 and N=500
#   ./run_local.sh --parallel 4             # run all with 4 concurrent jobs (GNU parallel)
#   ./run_local.sh --parallel 4 250 500     # specific sizes with parallelism
#   ./run_local.sh --baseline 250           # baseline (no cavity) for N=250
#
# Requires: Python with config.py importable from this directory.

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

PARALLEL_JOBS=0
BASELINE=false
SIZES=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --parallel)
            PARALLEL_JOBS="$2"
            shift 2
            ;;
        --baseline)
            BASELINE=true
            shift
            ;;
        *)
            SIZES+=("$1")
            shift
            ;;
    esac
done

# If no sizes specified, read all from config.py
if [[ ${#SIZES[@]} -eq 0 ]]; then
    SIZES=($(python3 -c "from config import SYSTEM_SIZES; print(' '.join(str(s) for s in SYSTEM_SIZES))"))
fi

echo "=== Finite-Size Scaling Campaign (local) ==="
echo "Sizes: ${SIZES[*]}"
echo "Baseline: $BASELINE"
echo "Parallel jobs: $PARALLEL_JOBS (0 = sequential)"
echo ""

BASELINE_FLAG=""
if $BASELINE; then
    BASELINE_FLAG="--no-cavity"
fi

run_size() {
    local N=$1
    local NREP
    NREP=$(python3 -c "from config import num_replicas; print(num_replicas($N))")
    echo "[N=$N] Running $NREP replicas ..."

    for ((r=0; r<NREP; r++)); do
        echo "  [N=$N, replica=$r/$NREP]"
        python3 run_single.py --N "$N" --replica "$r" $BASELINE_FLAG
    done
}

run_size_parallel() {
    local N=$1
    local NREP
    NREP=$(python3 -c "from config import num_replicas; print(num_replicas($N))")
    echo "[N=$N] Submitting $NREP replicas to GNU parallel ($PARALLEL_JOBS workers) ..."

    seq 0 $((NREP - 1)) | parallel -j "$PARALLEL_JOBS" \
        "python3 run_single.py --N $N --replica {} $BASELINE_FLAG"
}

for N in "${SIZES[@]}"; do
    if [[ $PARALLEL_JOBS -gt 0 ]]; then
        if command -v parallel &>/dev/null; then
            run_size_parallel "$N"
        else
            echo "WARNING: GNU parallel not found; falling back to sequential"
            run_size "$N"
        fi
    else
        run_size "$N"
    fi
done

echo ""
echo "=== Campaign complete ==="
