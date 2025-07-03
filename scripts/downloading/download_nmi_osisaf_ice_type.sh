#!/bin/bash
#PBS -N DL_OSISAF_itype
#PBS -P gv90
#PBS -q copyq
#PBS -l walltime=10:00:00
#PBS -l ncpus=1
#PBS -l mem=4GB
#PBS -l storage=gdata/gv90
#PBS -o DL_OSISAF_itype.out
#PBS -e DL_OSISAF_itype.err
yr0=2005
yrN=2022
hems=("sh" "nh")
Dloc="/g/data/gv90/da1339/SeaIce/OSI_SAF/ice_type"
Durl="https://thredds.met.no/thredds/fileServer/osisaf/met.no/ice/type/"
for hem in "${hems[@]}"; do
  for yr in $(seq $yr0 $yrN); do
    for mo in $(seq -w 01 12); do
      for dy0 in $(seq -w 01 31); do
        if [ "${dy0}" -ge "31" ] && [ "${mo}" = "04" -o "${mo}" = "06" -o "${mo}" = "09" -o "${mo}" = "11" ]; then continue
        elif [ "${dy0}" -ge "30" ] && [ "${mo}" = "02" ]; then
          if [ $((yr % 4)) -ne 0 ] || { [ $((yr % 100)) -eq 0 ] && [ $((yr % 400)) -ne 0 ]; }; then
            continue
          fi
        fi
        Furl="${Durl}/${yr}/${mo}/ice_type_${hem}_polstere-100_multi_${yr}${mo}${dy0}1200.nc"
        Floc="${Dloc}/${hem}/$(basename ${Furl})"
        # Check if the file already exists
        if [ -e "${Floc}" ]; then
          echo "File ${Floc} already exists. Skipping..."
        else
          mkdir -p "$(dirname "${Floc}")"
          echo "Downloading ${Floc}..."
          echo "Attempting: ${Furl}"
          if wget --quiet --content-disposition "${Furl}" -P "$(dirname "${Floc}")"; then
            echo "[OK] Downloaded: $(basename "${Furl}")"
          else
            echo "[FAIL] Could not download: ${Furl}" >&2
          fi
        fi
      done
    done
  done
done
echo "Download complete!"