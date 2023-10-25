#!/bin/bash
#PBS -P jk72
#PBS -l walltime=48:00:00
#PBS -q megamembw
#PBS -l mem=3000GB
#PBS -l ncpus=32
#PBS -l storage=gdata/cj50+gdata/jk72+gdata/ik11+gdata/hh5+gdata/rt52+gdata/gb6
#PBS -M daniel.atwater@utas.edu.au
module purge
module use /g/data/hh5/public/modules
module load conda/analysis3-unstable
python3 ~/AFIM/src/python/regrid_era5_for_cice6.py
