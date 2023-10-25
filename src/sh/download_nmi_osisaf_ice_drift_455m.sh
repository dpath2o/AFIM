yr0=1991
yrN=2020
hems=("sh" "nh")
Dloc="/g/data/jk72/da1339/SeaIce/nmi/osisaf/ice_drift_455m"
Durl="https://thredds.met.no/thredds/fileServer/osisaf/met.no/reprocessed/ice/drift_455m_files/merged"
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
        Furl="${Durl}/${yr}/${mo}/ice_drift_${hem}_ease2-750_cdr-v1p0_24h-${yr}${mo}${dy0}1200.nc"
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
echo "Download complete!"