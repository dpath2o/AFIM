#!/bin/bash
SIM_NAME=${1}
ISPD_THRESH=${2}
ISPD_TYPE=${3:-None}
for yr in {1993..1999}; do
    ./process_daily_fast_ice_wrapper.sh ${SIM_NAME} $yr ${ISPD_THRESH} ${ISP_TYPE}
done
