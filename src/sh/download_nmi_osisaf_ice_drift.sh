#!/bin/bash
yr0=2009
yrN=2022
hems=("sh" "nh")
Dloc="/g/data/jk72/da1339/SeaIce/nmi/osisaf/ice_drift"
Durl="https://thredds.met.no/thredds/fileServer/osisaf/met.no/ice/drift_lr/merged"
for hem in "${hems[@]}"; do
  for yr in $(seq $yr0 $yrN); do
    for mo in $(seq -w 01 12); do
      if [ "${yr}" -lt "2013" ] || { [ "${yr}" -eq "2013" ] && [ "${mo}" -lt "03" ]; } && [ "${hem}" = "sh" ]; then
        continue
      fi
      for dy0 in $(seq -w 01 31); do
        dyN=$(printf "%02d" $((10#$dy0 + 2)))  # Increment day by 2 and zero-pad
        if [ "${dyN}" -ge "32" ] && [ "${mo}" = "01" -o "${mo}" = "03" -o "${mo}" = "05" -o "${mo}" = "07" -o "${mo}" = "08" -o "${mo}" = "10" -o "${mo}" = "12" ]; then continue
        elif [ "${dyN}" -ge "31" ] && [ "${mo}" = "04" -o "${mo}" = "06" -o "${mo}" = "09" -o "${mo}" = "11" ]; then continue
        elif [ "${dyN}" -ge "30" ] && [ "${mo}" = "02" ]; then
          if [ $((yr % 4)) -ne 0 ] || { [ $((yr % 100)) -eq 0 ] && [ $((yr % 400)) -ne 0 ]; }; then
            continue
          fi
        fi
        Furl="${Durl}/${yr}/${mo}/ice_drift_${hem}_polstere-625_multi-oi_${yr}${mo}${dy0}1200-${yr}${mo}${dyN}1200.nc"
        Floc="${Dloc}/${hem}/$(basename ${Furl})"
        # Check if the file already exists
        if [ -e "${Floc}" ]; then
          echo "File ${Floc} already exists. Skipping..."
        else
          mkdir -p "$(dirname "${Floc}")"
          echo "Downloading ${Floc}..."
          wget --content-disposition "${Furl}" -P "$(dirname ${Floc})"
        fi
      done
    done
  done
done