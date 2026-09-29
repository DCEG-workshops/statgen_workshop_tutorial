#!/bin/bash
set -eo pipefail
[[ -n "${SLURM_JOB_ID:-}" ]] || { echo 'Use a Biowulf compute allocation.' >&2; exit 2; }
source "$(dirname "$0")/module-init.sh"
module load ldsc/1.0.1-20200724
if [[ "${1:-}" == --check ]]; then command -v ldsc.py; command -v munge_sumstats.py; exit; fi
case "${1:-}" in
  ldsc.py|munge_sumstats.py) script="$1"; shift; exec "$script" "$@" ;;
  *) echo 'Expected ldsc.py or munge_sumstats.py.' >&2; exit 2 ;;
esac
