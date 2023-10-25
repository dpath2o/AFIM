#!/bin/bash

D_graph="/g/data/jk72/da1339/GRAPHICAL"
model_run_names=("CICE6-NILbath" "CICE6-bath" "CICE6-FI_tensile" "CICE6-GI_10m" "CICE6-GI_neg0p2m" "CICE6-GI_NIL_basal")
geographical_regions=("circumpolar" "mawson" "olav" "sabrina")
graphical_types=("aice_m" "atmspd_m" "congel_m" "frazil_m" "hi_m" "iage_m" "icespd_m" "sice_m" "snow_m" "sst_m" "strairx_m" "strcorx_m" "strength_m" "strintx_m" "strocnx_m" "taubx_m")
numerical_months=("01" "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12")
abbreviated_months=("Jan" "Feb" "Mar" "Apr" "May" "Jun" "Jul" "Aug" "Sep" "Oct" "Nov" "Dec")

for run_name in "${model_run_names[@]}"; do
    for region in "${geographical_regions[@]}"; do
        for graphical_type in "${graphical_types[@]}"; do
            for month in "${numerical_months[@]}"; do
                D_ani="${D_graph}/${run_name}/${region}/${graphical_type}/monthly/*.png"
                F_out="${D_graph}/animations/CICE6/${run_name}_${graphical_type}_${region}.mp4"
                echo "Attempting to create animation for: ${D_ani}"
                /g/data/hh5/public/apps/miniconda3/envs/analysis3-23.04/bin/ffmpeg -framerate 10 -pattern_type glob -i "$D_ani" -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" -c:v libx264 -y -r 30 -pix_fmt yuv420p "$F_out"
                if [ -f "$F_out" ]; then
                    echo "Created animation for ${graphical_type}, ${region}, Month ${month}"
                else
                    echo "Something went wrong! Did not create animation for ${graphical_type}, ${region}, Month ${month}"
                fi
            done
        done
    done
done



