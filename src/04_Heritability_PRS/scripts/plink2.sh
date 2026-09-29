#!/bin/bash
set -eo pipefail
[[ -n "${SLURM_JOB_ID:-}" ]] || { echo 'Use a Biowulf compute allocation.' >&2; exit 2; }
source "$(dirname "$0")/module-init.sh"
module load plink/7.1-alpha
if [[ "${1:-}" == --check ]]; then command -v plink2; exit; fi
exec plink2 "$@"
