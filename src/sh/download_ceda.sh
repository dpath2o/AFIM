#!/bin/bash
levs=("L2P" "L3C")
sats=("cryosat2" "envisat")
hems=("sh" "nh")
Dloc="/g/data/jk72/da1339/SeaIce/ceda"
Durl="https://dap.ceda.ac.uk/neodc/esacci/sea_ice/data/sea_ice_thickness"
for lev in "${levs[@]}"; do
  for sat in "${sats[@]}"; do
    if [ "${sat}" == "cryosat2" ]; then
      yr0=2010
      yrN=2017
    else
      yr0=2002
      yrN=2012
    fi
    for hem in "${hems[@]}"; do
      for yr in $(seq $yr0 $yrN); do
        mo0=01
        moN=12
        if [ "${yr}" == "2010" ] && [ "${sat}" == "cryosat2" ]; then
          mo0=11
        elif [ "${yr}" == "2017" ] && [ "${sat}" == "cryosat2" ]; then
          moN=04
        elif [ "${yr}" == "2002" ] && [ "${sat}" == "envisat" ]; then
          mo0=10
        elif [ "${yr}" == "2012" ] && [ "${sat}" == "envisat" ]; then
          moN=03
        fi
        for mo in $(seq -w $mo0 $moN); do
          for dy in $(seq -w 01 31); do
            if [ "${dy}" -ge "31" ] && [ "${mo}" == "04" -o "${mo}" == "06" -o "${mo}" == "09" -o "${mo}" == "11" ]; then continue
            elif [ "${dy}" -ge "29" ] && [ "${mo}" == "02" ]; then
              if [ $((yr % 4)) -ne 0 ] || { [ $((yr % 100)) -eq 0 ] && [ $((yr % 400)) -ne 0 ]; }; then
                continue
              fi
            fi
            if [ "$lev" == "L2P" ] && [ "$sat" == "cryosat2" ]; then
              Furl="ESACCI-SEAICE-${lev}-SITHICK-SIRAL_${sat^^}-${hem^^}-${yr}${mo}${dy}-fv2.0.nc" #?download=1"
            elif [ "$lev" == "L2P" ] && [ "$sat" == "envisat" ]; then
              Furl="ESACCI-SEAICE-${lev}-SITHICK-RA2_${sat^^}-${hem^^}-${yr}${mo}${dy}-fv2.0.nc" #?download=1"
            elif [ "$lev" == "L3C" ] && [ "$sat" == "cryosat2" ] && [ "$hem" == "nh" ]; then
              Furl="ESACCI-SEAICE-${lev}-SITHICK-SIRAL_${sat^^}-${hem^^}25KMEASE2-${yr}${mo}-fv2.0.nc" #?download=1"
            elif [ "$lev" == "L3C" ] && [ "$sat" == "cryosat2" ] && [ "$hem" == "sh" ]; then
              Furl="ESACCI-SEAICE-${lev}-SITHICK-SIRAL_${sat^^}-${hem^^}50KMEASE2-${yr}${mo}-fv2.0.nc" #?download=1"
            elif [ "$lev" == "L3C" ] && [ "$sat" == "envisat" ] && [ "$hem" == "nh" ]; then
              Furl="ESACCI-SEAICE-${lev}-SITHICK-RA2_${sat^^}-${hem^^}25KMEASE2-${yr}${mo}-fv2.0.nc" #?download=1"
            elif [ "$lev" == "L3C" ] && [ "$sat" == "envisat" ] && [ "$hem" == "sh" ]; then
              Furl="ESACCI-SEAICE-${lev}-SITHICK-RA2_${sat^^}-${hem^^}50KMEASE2-${yr}${mo}-fv2.0.nc" #?download=1"
            fi
            Purl="${Durl}/${lev}/${sat}/v2.0/${hem^^}/${yr}/${mo}/${Furl}"
            Floc="${Dloc}/${lev}/${sat}/${hem}/$(basename "${Furl}")"
            if [ ! -f "$Floc" ]; then
              if wget --spider "$Purl" 2>/dev/null; then
                if ! wget --spider --quiet --method=HEAD "$Purl"; then
                  echo "File $Furl does not exist on the server."
                elif [ ! -f "$Floc" ]; then
                  mkdir -p "$(dirname "$Floc")"
                  echo "Downloading ${Furl}..."
                  wget "$Purl" -O "$Floc"
                fi
              else
                echo "URL $Purl is not valid."
              fi
            else
              echo "File $Floc exists ... skipping download."
            fi
          done
        done
      done
    done
  done
done
echo "Download complete!"
