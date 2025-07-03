#!/bin/bash
#PBS -N download-CMEMS-ORAS
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
for year in {2023..2025}; do
    F1 ="${year}0101_${year}0630_ORAS_org.nc"
    F3 ="${year}0701_${year}1231_ORAS_org.nc"
    D_="/g/data/gv90/da1339/afim_input/ORAS/daily/org"
    # Check if the output file already exists
    if [ -f "$D_/$F1" ]; then
        echo "$F1 already exists. Skipping..."
        continue
    fi
    dt0="${year}-01-01"
    dtN="${year}-06-30"
    ./copernicusmarine subset \
                       --dataset-id "cmems_mod_glo_phy-all_my_0.25deg_P1D-m" \
                       -f "$F1" \
                       -o "$D_" \
                       -t "$dt0" \
                       -T "$dtN" \
                       -v uo_oras -v vo_oras -v thetao_oras -v so_oras
    if [ -f "$D_/$F2" ]; then
        echo "$F2 already exists. Skipping..."
        continue
    fi
    dt0="${year}-07-01"
    dtN="${year}-12-31"
    ./copernicusmarine subset \
                       --dataset-id "cmems_mod_glo_phy-all_my_0.25deg_P1D-m" \
                       -f "$F2" \
                       -o "$D_" \
                       -t "$dt0" \
                       -T "$dtN" \
                       -v uo_oras -v vo_oras -v thetao_oras -v so_oras
    echo "Processed year: $year"
done
