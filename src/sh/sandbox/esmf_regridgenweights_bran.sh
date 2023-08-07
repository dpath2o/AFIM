#!/bin/bash
#PBS -P jk72
#PBS -l walltime=48:00:00
#PBS -q hugemem
#PBS -l mem=756GB
#PBS -l ncpus=48
#PBS -l storage=gdata/cj50+gdata/jk72+gdata/ik11+gdata/hh5+gdata/rt52+gdata/gb6
#PBS -M daniel.atwater@utas.edu.au
module purge
module load esmf
#module load openmpi
ulimit -s unlimited
D_grids="/g/data/jk72/da1339/grids"
P_G_BRAN=$D_grids"/BRAN/ocean_grid_for_regridding.nc"
P_G_CICE=$D_grids"/CICE/acom2_0p1_cice5_grid_for_regridding.nc"
P_W_BRAN=$D_grids"/weights/map_BRAN_tgrid_acom2_cice_tgrid_0p1_bilinear.nc"
ESMF_RegridWeightGen -p none -m bilinear -s $P_G_BRAN -d $P_G_CICE -w $P_W_BRAN