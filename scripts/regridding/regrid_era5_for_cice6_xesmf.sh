#!/bin/bash
#PBS -N regrid-ERA5-xesmf
#PBS -P jk72
#PBS -l walltime=10:00:00
#PBS -q hugemembw
#PBS -l mem=1020GB
#PBS -l ncpus=28
#PBS -l storage=gdata/cj50+gdata/jk72+gdata/ik11+gdata/rt52+gdata/gb6+gdata/gv90+gdata/xp65
#PBS -M daniel.atwater@utas.edu.au
# Modules setup
module purge
module use /g/data/xp65/public/modules
module load conda/analysis3-25.05
module load nco
python3 ~/AFIM/src/python/sandbox/regrid_era5_for_cice6.py
