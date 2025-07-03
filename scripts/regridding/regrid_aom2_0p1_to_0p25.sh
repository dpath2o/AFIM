#!/bin/bash
#PBS -N regrid-AOM2-ncl
#PBS -P jk72
#PBS -l walltime=48:00:00
#PBS -q normalbw
#PBS -l mem=256GB
#PBS -l ncpus=28
#PBS -l storage=gdata/cj50+gdata/jk72+gdata/ik11+gdata/hh5+gdata/rt52+gdata/gb6
#PBS -M daniel.atwater@utas.edu.au
# Modules setup
module purge
module use /g/data/hh5/public/modules
module load conda/analysis3
module load nco
module load ncl
ncl ~/AFIM/src/ncl/regrid_aom2_0p1_to_0p25.ncl
#python3 ~/AFIM/src/python/sandbox/regrid_aom2_0p1_to_0p25.py
