#!/usr/bin/env bash
set -euo pipefail
REPO="$(cd "$(dirname "$0")/.." && pwd)"
# shellcheck source=/dev/null
source "$REPO/scripts/activate_openmm.sh"
cd "$REPO"
exec jupyter lab --ip=0.0.0.0 --port="${JUPYTER_PORT:-8888}" --no-browser --allow-root "$@"
