#!/bin/bash
#PBS -N afim-analysis-zero_vel_days-alt
#PBS -P jk72
#PBS -l walltime=24:00:00
#PBS -q hugemembw
#PBS -l mem=256GB
#PBS -l ncpus=7
#PBS -l storage=gdata/cj50+gdata/jk72+gdata/ik11+gdata/hh5+gdata/rt52+gdata/gb6
#PBS -M daniel.atwater@utas.edu.au
module purge
module use /g/data/hh5/public/modules
module load conda/analysis3
python3 ~/AFIM/src/python/plot_cice_zero_velocity_days_alt.py