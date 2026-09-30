#!/bin/bash
set -eo pipefail
[[ -n "${SLURM_JOB_ID:-}" ]] || { echo 'Use a Biowulf compute allocation.' >&2; exit 2; }
source "$(dirname "$0")/module-init.sh"
module load GCTA/1.94.3
if [[ "${1:-}" == --check ]]; then command -v gcta64; exit; fi
exec gcta64 "$@"
