#!/bin/bash

D_graph="/g/data/jk72/da1339/GRAPHICAL"
model_run_names=("CICE6-NILbath" "CICE6-bath" "CICE6-FI_tensile" "CICE6-GI_10m" "CICE6-GI_neg0p2m" "CICE6-GI_NIL_basal")
geographical_regions=("circumpolar" "mawson" "olav" "sabrina")
graphical_types=("zero_vel_days" "frazil_grow_days")
numerical_months=("01" "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12")
abbreviated_months=("Jan" "Feb" "Mar" "Apr" "May" "Jun" "Jul" "Aug" "Sep" "Oct" "Nov" "Dec")

for region in "${geographical_regions[@]}"; do
    for graphical_type in "${graphical_types[@]}"; do
        for ((i = 0; i < ${#numerical_months[@]}; i++)); do
            num_mon="${numerical_months[$i]}"
            alpha_mon="${abbreviated_months[$i]}"
            
            F_inputs=()
            for run_name in "${model_run_names[@]}"; do
                F_new="${D_graph}/${run_name}/${region}/${graphical_type}/monthly/month${num_mon}_${alpha_mon}_${run_name}_${graphical_type}_SH_jra55do_aom2_0p25.png"
                if [ -f "$F_new" ]; then
                    F_inputs+=("$F_new")
                else
                    echo "Input file not found: ${F_new}"
                fi
            done
            F_out="${D_graph}/${graphical_type}_${region}_mon${num_mon}.png"
            echo "attempting to montage the following files: ${F_inputs[@]}"
            montage "${F_inputs[@]}" -quality 100 -depth 16 -tile 2x -geometry 700x500! "$F_out"
            if [ -f "$F_out" ]; then
                echo "Created montage for ${graphical_type}, ${region}, Month ${num_mon}"
            else
                echo "something went wrong! did not create montage for ${graphical_type}, ${region}, Month ${num_mon}"
            fi
        done
    done
done
