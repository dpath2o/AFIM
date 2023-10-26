#!/bin/bash
#PBS -N regrid-era5
#PBS -P jk72
#PBS -l walltime=48:00:00
#PBS -q hugemembw
#PBS -l mem=1020GB
#PBS -l ncpus=28
#PBS -l storage=gdata/cj50+gdata/jk72+gdata/ik11+gdata/hh5+gdata/rt52+gdata/gb6
#PBS -M daniel.atwater@utas.edu.au
# Modules setup
module purge
module use /g/data/hh5/public/modules
#module load conda/analysis3-unstable
module load nco
# Setup logging
LOGFILE="/home/581/da1339/logs/regrid_era5_for_cice6.log"
exec 3>&1 1>>${LOGFILE} 2>&1
# Variables
G_origin="aom2"
G_res="0p25"
regrid_method="nearest_s2d"
F_weights="/g/data/jk72/da1339/grids/weights/map_ERA5_${G_origin}_${G_res}_${regrid_method}.nc"
G_CICE="/g/data/ik11/inputs/access-om2/input_20200530/cice_025deg/grid.nc"
D_era5="/g/data/rt52/era5/single-levels/reanalysis"
D_frcg="/scratch/jk72/da1339/cice-dirs/input/AFIM/forcing/0p25/ERA5/24XDAILY"
D_tmp="/scratch/jk72/da1339/afim_input/ERA5/0p25"
var_names_in=("2t" "msdwlwrf" "msdwswrf" "mtpr" "10u" "10v" "sp")
var_names_out=("airtmp" "dlwsfc" "glbrad" "ttlpcp" "wndewd" "wndnwd" "spchmd")
cd $D_tmp
# Loop through years: 2005-2018
for yr_str in {2005..2018}; do
    echo "$(date) : DEBUG : Starting regridding for year ${yr_str}"
    # loop through months: 01-12
    for month in {01..12}; do
        # Note: Ensure that the ERA5 files are named in the format yyyyMMdd-* in the directory. 
        # Regridding individual variables
        for index in "${!var_names_in[@]}"; do
            var_in="${var_names_in[$index]}"
            var_out="${var_names_out[$index]}"
            if [ "$var_in" == "sp" ]; then
                echo "$(date) : INFO : Starting regridding ${var_in} for year ${yr_str} month ${month}"
                ncremap -t 28 -i ${D_era5}/${var_in}/${yr_str}/${yr_str}${month}01-*.nc -m $F_weights -o sp_reG_${month}.nc
                echo "$(date) : INFO : Starting regridding 2d for year ${yr_str} month ${month}"
                ncremap -t 28 -i ${D_era5}/2d/${yr_str}/${yr_str}${month}01-*.nc -m $F_weights -o d2m_reG_${month}.nc
                # Merge sp_regridded and d2m_regridded for specific humidity calculation
                echo "$(date) : INFO : Merging sp_regridded and d2m_regridded for year ${yr_str} month ${month}"
                ncks -t 28 -A sp_reG_${month}.nc d2m_reG_${month}.nc -o merged_${month}.nc
                # Compute specific humidity
                # Define constants
                echo "$(date) : INFO : Setting constants for specific humidity calculation for year ${yr_str} month ${month}"
                ncap2 -O -s "Rdry=287.0597" -s "Rvap=461.5250" -s "a1=611.21" -s "a3=17.502" -s "a4=32.19" -s "T0=273.16" merged_${month}.nc intermediate1_${month}.nc
                # Calculate E
                echo "$(date) : INFO : computing E for specific humidity calculation for year ${yr_str} month ${month}"
                ncap2 -O -s "E=a1*exp(a3*(d2m-T0)/(d2m-a4))" intermediate1_${month}.nc intermediate2_${month}.nc
                # Compute qsat
                echo "$(date) : INFO : computing specific humidity for year ${yr_str} month ${month}"
                ncap2 -O -s "spchmd=(Rdry/Rvap)*E/(sp-((1-Rdry/Rvap)*E))" intermediate2_${month}.nc ${var_out}_${month}.nc
                ncatted -a units,spchmd,o,c,"kg/kg" ${var_out}_${month}.nc
            else
                echo "$(date) : INFO : Starting regridding ${var_in} for year ${yr_str} month ${month}"
                ncremap -t 28 -i ${D_era5}/${var_in}/${yr_str}/${yr_str}${month}01-*.nc -m $F_weights -o ${var_out}_reG_${month}.nc
                ncrename -v ${var_in},${var_out} ${var_out}_reG_${month}.nc
            fi
        done
    done
    # Merge into a single yearly file
    echo "$(date) : INFO : concating monthly files into a yearly file for ${yr_str}"
    ncrcat airtmp_reG_*.nc dlwsfc_reG_*.nc glbrad_reG_*.nc ttlpcp_reG_*.nc wndewd_reG_*.nc wndnwd_reG_*.nc spchmd_*.nc -O yearly_output_${yr_str}.nc
    # Clean-up monthly files
    rm airtmp_reG_*.nc dlwsfc_reG_*.nc glbrad_reG_*.nc ttlpcp_reG_*.nc wndewd_reG_*.nc wndnwd_reG_*.nc spchmd_*.nc sp_reG_*.nc d2m_reG_*.nc merged_*.nc int*_*.nc
    echo "$(date) : INFO : removed temporary monthly files"
    # Additional metadata
    ncatted -a author,global,o,c,'dpath2o, daniel.atwater@utas.edu.au' yearly_output_${yr_str}.nc
    # Final output path
    F_out="era5_for_cice6_${yr_str}.nc"
    P_out="${D_frcg}/${F_out}"
    # Move to the final path
    mv yearly_output_${yr_str}.nc $P_out
    echo "$(date) : INFO : Finished writing to ${P_out}"
done
