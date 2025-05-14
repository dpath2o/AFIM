#!/bin/bash
for yr in {1993..1999}; do
    ./process_daily_fast_ice_wrapper.sh ${1} $yr
done
