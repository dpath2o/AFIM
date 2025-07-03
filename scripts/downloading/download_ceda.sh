#!/bin/bash
#PBS -N DL_ESACCI
#PBS -P gv90
#PBS -q copyq
#PBS -l walltime=10:00:00
#PBS -l ncpus=1
#PBS -l mem=4GB
#PBS -l storage=gdata/gv90
#PBS -o DL_ESACCI.out
#PBS -e DL_ESACCI.err
umask 002
levs=("L2P" "L3C")
sats=("envisat" "cryosat2")
hems=("sh" "nh")
Dloc="/g/data/gv90/da1339/SeaIce/ESA_CCI"
Durl="https://dap.ceda.ac.uk/neodc/esacci/sea_ice/data/sea_ice_thickness"
for lev in "${levs[@]}"; do
  for sat in "${sats[@]}"; do
    if [ "$sat" == "cryosat2" ]; then
      yr0=2010; yrN=2017
    else
      yr0=2002; yrN=2012
    fi
    for hem in "${hems[@]}"; do
      for yr in $(seq $yr0 $yrN); do
        mo0=01; moN=12
        [ "$sat" == "cryosat2" ] && [ "$yr" == "2010" ] && mo0=11
        [ "$sat" == "cryosat2" ] && [ "$yr" == "2017" ] && moN=04
        [ "$sat" == "envisat" ] && [ "$yr" == "2002" ] && mo0=10
        [ "$sat" == "envisat" ] && [ "$yr" == "2012" ] && moN=03
        for mo in $(seq -w $mo0 $moN); do
          for dy in $(seq -w 01 31); do
            case "$mo" in
              04|06|09|11) [ "$dy" == "31" ] && continue ;;
              02)
                if [ "$dy" -ge 29 ]; then
                  if ! date -d "$yr-$mo-$dy" &>/dev/null; then
                    continue
                  fi
                fi
                ;;
            esac
            case "$lev-$sat-$hem" in
              L2P-cryosat2-*) Furl="ESACCI-SEAICE-L2P-SITHICK-SIRAL_${sat^^}-${hem^^}-${yr}${mo}${dy}-fv2.0.nc" ;;
              L2P-envisat-*)  Furl="ESACCI-SEAICE-L2P-SITHICK-RA2_${sat^^}-${hem^^}-${yr}${mo}${dy}-fv2.0.nc" ;;
              L3C-cryosat2-nh) Furl="ESACCI-SEAICE-L3C-SITHICK-SIRAL_${sat^^}-${hem^^}25KMEASE2-${yr}${mo}-fv2.0.nc" ;;
              L3C-cryosat2-sh) Furl="ESACCI-SEAICE-L3C-SITHICK-SIRAL_${sat^^}-${hem^^}50KMEASE2-${yr}${mo}-fv2.0.nc" ;;
              L3C-envisat-nh)  Furl="ESACCI-SEAICE-L3C-SITHICK-RA2_${sat^^}-${hem^^}25KMEASE2-${yr}${mo}-fv2.0.nc" ;;
              L3C-envisat-sh)  Furl="ESACCI-SEAICE-L3C-SITHICK-RA2_${sat^^}-${hem^^}50KMEASE2-${yr}${mo}-fv2.0.nc" ;;
              *) echo "Invalid combination: $lev $sat $hem" >&2; continue ;;
            esac
            Purl="${Durl}/${lev}/${sat}/v2.0/${hem^^}/${yr}/${mo}/${Furl}"
            Floc="${Dloc}/${lev}/${sat}/${hem}/${Furl}"
            if [ -f "$Floc" ]; then
              echo "[SKIP] Already exists: $Floc"
              continue
            fi
            mkdir -p "$(dirname "$Floc")"
            logfile="${Dloc}/download_log.txt"
            echo "[TRY] $Furl" | tee -a "$logfile"
            wget --quiet --timeout=30 --tries=3 -O "$Floc" "$Purl"
            exit_code=$?
            if [ $exit_code -eq 0 ]; then
              echo "[OK] Downloaded $Furl" | tee -a "$logfile"
            else
              echo "[FAIL] $Furl (exit code $exit_code)" | tee -a "$logfile"
              rm -f "$Floc"
            fi
            sleep 1  # Prevent server hammering
          done
        done
      done
    done
  done
done
echo "All downloads attempted. Script finished."

