#!/bin/bash
#PBS -P jk72
#PBS -l walltime=48:00:00
#PBS -q normal
#PBS -l mem=40GB
#PBS -l ncpus=8
#PBS -l storage=gdata/cj50+gdata/jk72+gdata/ik11+gdata/hh5+gdata/rt52+gdata/gb6+gdata/qv56
#PBS -M daniel.atwater@utas.edu.au
U_base="https://noaadata.apps.nsidc.org/NOAA/G02202_V4"
D_base="/g/data/jk72/da1339/SeaIce/nsdic/G02202_V4"
hemis=("south" "north")
freqs=("monthly" "aggregate" "daily")
YR0=1979
YRN=2022
for freq in "${freqs[@]}"; do
    for hemi in "${hemis[@]}"; do
        if [ "$hemi" = "south" ]; then
            hem="sh"
        else
            hem="nh"
        fi
        for yr in $(seq "$YR0" "$YRN"); do
            if [ "$freq" = "monthly" ]; then
                U_down="${U_base}/${hemi}/${freq}"
                D_down="${D_base}/${hemi}/${freq}"
                cd "$D_down"
                for mo in {01..12}; do
                    if [ $yr -lt 1988 ]; then
                        if [ $yr -eq 1987 and $mo -gt 7 ]; then
                            F_name="seaice_conc_monthly_${hem}_${yr}${mo}_f08_v04r00.nc"    
                        else
                            F_name="seaice_conc_monthly_${hem}_${yr}${mo}_n07_v04r00.nc"
                        fi
                    elif [ $yr -lt 1992 ]; then
                        F_name="seaice_conc_monthly_${hem}_${yr}${mo}_f08_v04r00.nc"
                    elif [ $yr -lt 1996 ]; then
                        if [ $yr -eq 1995 and $mo -gt 9 ]; then
                            F_name="seaice_conc_monthly_${hem}_${yr}${mo}_f13_v04r00.nc"
                        else
                            F_name="seaice_conc_monthly_${hem}_${yr}${mo}_f11_v04r00.nc"
                        fi
                    elif [ $yr -lt 2008 ]; then
                        F_name="seaice_conc_monthly_${hem}_${yr}${mo}_f13_v04r00.nc"
                    else
                        F_name="seaice_conc_monthly_${hem}_${yr}${mo}_f17_v04r00.nc"
                    fi
                    if [[ -s "${F_name}" && $(stat -c%s "${F_name}") -gt 100000 ]]; then
                        echo "${F_name} exists. Skipping download."
                    else
                        echo "Downloading ${U_down}/${F_name}"
                        wget -nc "${U_down}/${F_name}"
                    fi
                done
            elif [ "$freq" = "daily" ]; then
                U_down="${U_base}/${hemi}/${freq}/${yr}"
                D_down="${D_base}/${hemi}/${freq}/${yr}"
                if [ ! -d "$D_down" ]; then
                    mkdir -pv "$D_down"
                fi
                cd "$D_down"
                wget -q "${U_down}/"
                F_index="index.html"
                F_names=($(grep -oP '(?<=<a href=")[^"]+' "$F_index"))
                # Loop through the filenames and download the files
                for F_name in "${F_names[@]}"; do
                    if [[ -s "${F_name}" && $(stat -c%s "${F_name}") -gt 100000 ]]; then
                        echo "${F_name} exists. Skipping download."
                    else
                        echo "Downloading ${F_name}"
                        wget -nc "${U_down}/${F_name}"
                    fi
                done
                rm "${F_index}"
            elif [ "$freq" = "aggregate" ]; then
                U_down="${U_base}/${hemi}/${freq}"
                D_down="${D_base}/${hemi}/${freq}"
                F_name="seaice_conc_daily_${hem}_${yr}_v04r00.nc"
                cd "$D_down"
                if [[ -s "${F_name}" && $(stat -c%s "${F_name}") -gt 100000 ]]; then
                    echo "${F_name} exists. Skipping download."
                else
                    echo "Downloading ${F_name}"
                    wget -nc "${U_down}/${F_name}"
                fi
            fi
        done
    done
done
