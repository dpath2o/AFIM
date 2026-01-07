#!/bin/bash
SIM_ROOT="/g/data/gv90/da1339/afim_output"
for sim_dir in "$SIM_ROOT"/*/; do
    sim_name=$(basename "$sim_dir")
    if [[ "$sim_name" != "AOM2-ERA5" && "$sim_name" != "pack_ice.zarr" && -d "$sim_dir" ]]; then
        ./SI_metrics_pbs_wrapper.sh -s "$sim_name" -t 5.0e-4 -i SI_BT -S 1993-01-01 -E 1999-12-31
    fi
done
