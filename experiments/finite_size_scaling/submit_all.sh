#!/bin/bash
# Submit all finite-size scaling jobs to SLURM (NYU Torch).
#
# Usage:
#   ./submit_all.sh                   # submit cavity runs for all sizes
#   ./submit_all.sh --baseline        # submit baseline (no-cavity) for all sizes
#   ./submit_all.sh --sizes 250 500   # specific sizes only
#   ./submit_all.sh --dry-run         # print commands without submitting

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

DRY_RUN=false
BASELINE=false
SIZES=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        --baseline)
            BASELINE=true
            shift
            ;;
        --sizes)
            shift
            while [[ $# -gt 0 && ! "$1" =~ ^-- ]]; do
                SIZES+=("$1")
                shift
            done
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

if [[ ${#SIZES[@]} -eq 0 ]]; then
    SIZES=($(python3 -c "from config import SYSTEM_SIZES; print(' '.join(str(s) for s in SYSTEM_SIZES))"))
fi

mkdir -p logs

echo "=== Finite-Size Scaling: SLURM Submission ==="
echo "Baseline: $BASELINE"
echo "Sizes: ${SIZES[*]}"
echo ""

# Wall-time estimates: scale with N (larger systems take longer per step)
get_walltime() {
    local N=$1
    if   [[ $N -le 500 ]];    then echo "12:00:00"
    elif [[ $N -le 2000 ]];   then echo "24:00:00"
    elif [[ $N -le 8000 ]];   then echo "48:00:00"
    elif [[ $N -le 32000 ]];  then echo "72:00:00"
    else                           echo "96:00:00"
    fi
}

get_mem() {
    local N=$1
    if   [[ $N -le 2000 ]];   then echo "32GB"
    elif [[ $N -le 16000 ]];  then echo "64GB"
    else                           echo "128GB"
    fi
}

for N in "${SIZES[@]}"; do
    NREP=$(python3 -c "from config import num_replicas; print(num_replicas($N))")
    ARRAY_SPEC="0-$((NREP - 1))"
    WALLTIME=$(get_walltime "$N")
    MEM=$(get_mem "$N")

    if $BASELINE; then
        JOBNAME="fss_base_N${N}"
        EXTRA="--no-cavity"
    else
        JOBNAME="fss_N${N}"
        EXTRA=""
    fi

    CMD="sbatch --job-name=${JOBNAME} --array=${ARRAY_SPEC} --time=${WALLTIME} --mem=${MEM} submit_slurm.sh ${N} ${EXTRA}"

    echo "$CMD"
    if ! $DRY_RUN; then
        eval "$CMD"
    fi
done

echo ""
if $DRY_RUN; then
    echo "(dry run -- no jobs submitted)"
else
    echo "All jobs submitted. Monitor with: squeue --me"
fi
