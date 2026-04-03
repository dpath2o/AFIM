#!/bin/bash
set -euo pipefail

PBS_SCRIPT="/home/581/da1339/AFIM/src/AFIM/scripts/waves/CAWCR_spec_reG.pbs"

usage() {
    echo "Usage: $0 YEAR MONTH [--ow_nc] [--ow_wgt]"
    exit 1
}

[[ $# -ge 2 ]] || usage

YEAR="$1"
MONTH="$(printf "%02d" "$2")"
shift 2

OW_NC="false"
OW_WGT="false"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --ow_nc)
            OW_NC="true"
            ;;
        --ow_wgt)
            OW_WGT="true"
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
    shift
done

echo "Submitting year ${YEAR} month ${MONTH} at $(date)"
echo "  overwrite nc      = ${OW_NC}"
echo "  overwrite weights = ${OW_WGT}"

qsub -N "cawcr_spec_reG_${YEAR}-${MONTH}" \
     -v YEAR="${YEAR}",MONTH="${MONTH}",OW_NC="${OW_NC}",OW_WGT="${OW_WGT}" \
     "${PBS_SCRIPT}"

echo "Submitted year ${YEAR} month ${MONTH} at $(date)"
