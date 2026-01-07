#!/bin/bash
SIM_ROOT="/g/data/gv90/da1339/afim_output"
for sim_dir in "$SIM_ROOT"/*/; do
    sim_name=$(basename "$sim_dir")
    if [[ "$sim_name" != "AOM2-ERA5" && "$sim_name" != "pack_ice.zarr" && -d "$sim_dir" ]]; then
        ./classify_PI_pbs_wrapper.sh -s "$sim_name" -t 5.0e-4 -i BT -S 1993-01-01 -E 1999-12-31 -d -r -y -k -z
    fi
done
