#!/bin/bash
#PBS -P jk72
#PBS -l walltime=48:00:00
#PBS -q hugemembw
#PBS -l mem=1020GB
#PBS -l ncpus=28
#PBS -l storage=gdata/cj50+gdata/jk72+gdata/ik11+gdata/hh5+gdata/rt52+gdata/gb6
#PBS -M daniel.atwater@utas.edu.au
module purge
module use /g/data/hh5/public/modules
module load conda/analysis3-22.01
python3 ~/src/python/acom2_coallate_for_cice_forcing.py
