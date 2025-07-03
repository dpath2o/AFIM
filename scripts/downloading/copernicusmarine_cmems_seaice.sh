#!/bin/bash
#PBS -N download-CMEMS-seaice
#PBS -P gv90
#PBS -l walltime=10:00:00
#PBS -q copyq
#PBS -l mem=10GB
#PBS -l ncpus=1
#PBS -l storage=gdata/xp65+gdata/gv90
#PBS -l wd
#PBS -m abe
#PBS -M daniel.atwater@utas.edu.au
moudle purge
module use /g/data/xp65/public/modules
module load conda/analysis3-25.05
cd ~/.local/bin
res="0p083" #0p25
datasetid="cmems_mod_glo_phy_my_0.083deg_P1D-m" #"cmems_mod_glo_phy-all_my_0.25deg_P1D-m"
D_="/g/data/gv90/da1339/SeaIce/CMEMS/${res}/daily/"
for year in {2023..2025}; do
    F_="${year}0101_${year}1231_CMEMS_org.nc"
    if [ -f "$D_/$F_" ]; then
        echo "$F_ already exists. Skipping..."
        continue
    fi
    dt0="${year}-01-01"
    dtN="${year}-12-31"
    ./copernicusmarine subset \
                       --dataset-id $datasetid\
                       -f "$F_" \
                       -o "$D_" \
                       -t "$dt0" \
                       -T "$dtN" \
                       -v siconc_oras -v sithick_oras \
                       -v siconc_cglo -v sithick_cglo \
                       -v siconc_glor -v sithick_glor
    echo "Processed year: $year"
done
