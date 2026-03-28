#!/bin/bash
set -euo pipefail
PBS_SCRIPT="/home/581/da1339/AFIM/src/AFIM/scripts/downloading/CAWCR_spec_by_year.pbs"
for YEAR in $(seq 1993 1993); do
    echo "Submitting year ${YEAR} at $(date)"
    qsub -v YEAR="${YEAR}" -W block=true "${PBS_SCRIPT}"
    echo "Completed year ${YEAR} at $(date)"
done
