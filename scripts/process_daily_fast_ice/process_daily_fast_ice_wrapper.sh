#!/bin/bash
for m in {01..12}; do
    start="${2}-${m}-01"
    end=$(date -d "$start +1 month -1 day" +%F)
    qsub -N fi_${1}_${2}_${m} -v SIM_NAME=${1},MONTH=${m},YEAR=${2},START_DATE=${start},END_DATE=${end} process_daily_fast_ice_execute.pbs
done
