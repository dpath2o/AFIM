#!/bin/bash
SIM_NAME=${1}
YEAR=${2}
ISPD_THRESH=${3}
ISPD_TYPE=${4}

for m in {01..12}; do
    DT0="${YEAR}-${m}-01"
    DTN=$(date -d "${DT0} +1 month -1 day" +%F)
    if [ -z "$ISPD_TYPE" ]; then
        # No ISPD_TYPE — omit from -v and from CLI
        qsub -N fi_${SIM_NAME}_${YEAR}_${m} \
            -v SIM_NAME=${SIM_NAME},MONTH=${m},YEAR=${YEAR},START_DATE=${DT0},END_DATE=${DTN},ISPD_THRESH=${ISPD_THRESH} \
            process_daily_fast_ice_execute.pbs
    else
        # Include ISPD_TYPE
        qsub -N fi_${SIM_NAME}_${YEAR}_${m} \
            -v SIM_NAME=${SIM_NAME},MONTH=${m},YEAR=${YEAR},START_DATE=${DT0},END_DATE=${DTN},ISPD_THRESH=${ISPD_THRESH},ISPD_TYPE=${ISPD_TYPE} \
            process_daily_fast_ice_execute.pbs
    fi
done

