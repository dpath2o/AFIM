#!/bin/bash

U_base="https://noaadata.apps.nsidc.org/NOAA/G02202_V4"
D_base="/g/data/jk72/da1339/SeaIce/nsidc/G02202_V4"
hemis=("south" "north")
freqs=("daily" "monthly" "aggregate")
YR0=1979
YRN=2024
for freq in "${freqs[@]}"; do
    for hemi in "${hemis[@]}"; do
        if [ "${hemi}" = "south" ]; then
            hem="sh"
        else
            hem="nh"
        fi
        for yr in $(seq "$YR0" "$YRN"); do
            echo $yr
            if [ "${freq}" = "monthly" ]; then
                U_down="${U_base}/${hemi}/${freq}"
                D_down="${D_base}/${hemi}/${freq}"
                cd "$D_down"
                for mo in {01..12}; do
                    echo $mo
                    if [ $yr -lt 1988 ]; then
                        if [ $yr -eq 1987 ] && [ $mo -gt 7 ]; then
                            F_name="seaice_conc_monthly_${hem}_${yr}${mo}_f08_v04r00.nc"    
                        else
                            F_name="seaice_conc_monthly_${hem}_${yr}${mo}_n07_v04r00.nc"
                        fi
                    elif [ $yr -lt 1992 ]; then
                        F_name="seaice_conc_monthly_${hem}_${yr}${mo}_f08_v04r00.nc"
                    elif [ $yr -lt 1996 ]; then
                        if [ $yr -eq 1995 ] && [ $mo -gt 9 ]; then
                            F_name="seaice_conc_monthly_${hem}_${yr}${mo}_f13_v04r00.nc"
                        else
                            F_name="seaice_conc_monthly_${hem}_${yr}${mo}_f11_v04r00.nc"
                        fi
                    elif [ $yr -lt 2008 ]; then
                        F_name="seaice_conc_monthly_${hem}_${yr}${mo}_f13_v04r00.nc"
                    else
                        F_name="seaice_conc_monthly_${hem}_${yr}${mo}_f17_v04r00.nc"
                    fi
                    if [[ -s "${F_name}" && $(stat -c%s "${F_name}") -gt 10000 ]]; then
                        echo "${F_name} exists. Skipping download."
                    else
                        echo "Downloading ${U_down}/${F_name}"
                        wget -nc "${U_down}/${F_name}"
                    fi
                done
            elif [ "${freq}" = "daily" ]; then
                U_down="${U_base}/${hemi}/${freq}/${yr}"
                D_down="${D_base}/${hemi}/${freq}"
                if [ ! -d "$D_down" ]; then
                    mkdir -pv "$D_down"
                fi
                cd "${D_down}"
                for mo in {01..12}; do
                    case $mo in
                        01|03|05|07|08|10|12)
                            ndy=31
                            ;;
                        04|06|09|11)
                            ndy=30
                            ;;
                        02)
                            # Leap year calculation
                            if [ $(( yr % 4 )) -eq 0 ] && [ $(( yr % 100 )) -ne 0 ] || [ $(( yr % 400 )) -eq 0 ]; then
                                ndy=29
                            else
                                ndy=28
                            fi
                            ;;
                    esac
                    for dy in $(seq -w 01 $ndy); do
                        if [ $yr -lt 1988 ]; then
                            if [ $yr -eq 1987 ]; then
                                if [ $mo -eq 7 ] && [ $dy -gt 9 ]; then
                                    F_name="seaice_conc_daily_${hem}_${yr}${mo}${dy}_f08_v04r00.nc"
                                elif [ $mo -gt 7 ]; then
                                    F_name="seaice_conc_daily_${hem}_${yr}${mo}${dy}_f08_v04r00.nc"
                                else
                                    F_name="seaice_conc_daily_${hem}_${yr}${mo}${dy}_n07_v04r00.nc"
                                fi
                            else
                                F_name="seaice_conc_daily_${hem}_${yr}${mo}${dy}_n07_v04r00.nc"
                            fi
                        fi
                        if [ $yr -lt 1992 ] && [ $yr -gt 1987 ]; then
                            if [ $yr -eq 1991 ]; then
                                if [ $mo -eq 12 ] && [ $dy -gt 2 ]; then
                                    F_name="seaice_conc_daily_${hem}_${yr}${mo}${dy}_f11_v04r00.nc"
                                else
                                    F_name="seaice_conc_daily_${hem}_${yr}${mo}${dy}_f08_v04r00.nc"
                                fi
                            else
                                F_name="seaice_conc_daily_${hem}_${yr}${mo}${dy}_f08_v04r00.nc"
                            fi
                        fi
                        if [ $yr -lt 1996 ] && [ $yr -gt 1991 ]; then
                            if [ $yr -eq 1995 ] && [ $mo -gt 9 ]; then
                                F_name="seaice_conc_daily_${hem}_${yr}${mo}${dy}_f13_v04r00.nc"
                            else
                                F_name="seaice_conc_daily_${hem}_${yr}${mo}${dy}_f11_v04r00.nc"
                            fi
                        fi
                        if [ $yr -lt 2008 ] && [ $yr -gt 1995 ]; then
                            F_name="seaice_conc_daily_${hem}_${yr}${mo}${dy}_f13_v04r00.nc"
                        fi
                        if [ $yr -gt 2007 ]; then
                            F_name="seaice_conc_daily_${hem}_${yr}${mo}${dy}_f17_v04r00.nc"
                        fi
                        if [[ -f "${F_name}" && $(stat -c%s "${F_name}") -gt 10000 ]]; then
                            echo "${F_name} exists. Skipping download."
                        else
                            echo "Downloading ${U_down}/${F_name}"
                            wget -nc "${U_down}/${F_name}"
                        fi
                    done
                done
            elif [ "${freq}" = "aggregate" ]; then
                U_down="${U_base}/${hemi}/${freq}"
                D_down="${D_base}/${hemi}/${freq}"
                F_name="seaice_conc_daily_${hem}_${yr}_v04r00.nc"
                cd "${D_down}"
                if [[ -s "${F_name}" && $(stat -c%s "${F_name}") -gt 10000 ]]; then
                    echo "${F_name} exists. Skipping download."
                else
                    echo "Downloading ${F_name}"
                    wget -nc "${U_down}/${F_name}"
                fi
            fi
        done
    done
done
