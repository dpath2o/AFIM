#!/bin/bash
#PBS -N download-CMEMS-0p083-seaice
#PBS -P gv90
#PBS -l walltime=10:00:00
#PBS -q copyq
#PBS -l mem=10GB
#PBS -l ncpus=1
#PBS -l storage=gdata/xp65+gdata/gv90
#PBS -l wd
#PBS -m abe
#PBS -M daniel.atwater@utas.edu.au
module purge
module use /g/data/xp65/public/modules
module load conda/analysis3-25.05
cd ~/.local/bin
for year in {1993..2023}; do
    F_="${year}0101_${year}1231_CMEMS_0p083_org.nc"
    D_="/g/data/gv90/da1339/SeaIce/CMEMS/0p083/daily/"
    if [ -f "$D_/$F_" ]; then
        echo "$F_ already exists. Skipping..."
        continue
    fi
    dt0="${year}-01-01"
    dtN="${year}-12-31"
    ./copernicusmarine subset\
                       --dataset-id cmems_mod_glo_phy_my_0.083deg_P1D-m\
                       --dataset-version 202311\
                       --variable siconc\
                       --variable sithick\
                       --variable usi\
                       --variable vsi\
                       -f "$F_" \
                       -o "$D_" \
                       -t "$dt0" \
                       -T "$dtN" \
                       --minimum-longitude -180\
                       --maximum-longitude 179.9166717529297\
                       --minimum-latitude -80\
                       --maximum-latitude 90\
                       --minimum-depth 0.49402499198913574\
                       --maximum-depth 0.49402499198913574\
                       --coordinates-selection-method strict-inside\
                       --disable-progress-bar\
                       --log-level ERROR
    echo "Processed year: $year"
done
