#!/bin/bash
#PBS -P jk72
#PBS -l walltime=8:00:00
#PBS -l ncpus=384
#PBS -l mem=1TB
#PBS -l storage=gdata/cj50+gdata/jk72+gdata/ik11+gdata/hh5+gdata/rt52+gdata/gb6
#PBS -M daniel.atwater@utas.edu.au
module use /g/data/hh5/public/modules
module load conda/analysis3
source ~/.bashrc
conda activate afim
P_ERA5_grid=/g/data/jk72/da1339/grids/ERA5/ERA5_grid.nc
P_CICE_grid=/g/data/jk72/da1339/grids/CICE/g0p1_cice_tgrid_only_latlon.nc
P_weights=/g/data/jk72/da1339/grids/weights/map_ERA5_to_access-om2_cice_bilnr_pole_none.20221017.nc
ESMF_RegridWeightGen -m bilinear -p none -s ${P_ERA5_grid} -d ${P_CICE_grid} -w $P_weights