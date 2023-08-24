#!/bin/bash
#PBS -N afim-analysis
#PBS -P jk72
#PBS -l walltime=48:00:00
#PBS -q normal
#PBS -l mem=128GB
#PBS -l ncpus=8
#PBS -l storage=gdata/cj50+gdata/jk72+gdata/ik11+gdata/hh5+gdata/rt52+gdata/gb6
#PBS -M daniel.atwater@utas.edu.au
module purge
module use /g/data/hh5/public/modules
module load conda/analysis3
python3 ~/src/python/plot_cice_analysis.py