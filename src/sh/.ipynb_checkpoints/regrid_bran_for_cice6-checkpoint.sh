#!/bin/bash
#PBS -P jk72
#PBS -l walltime=48:00:00
#PBS -q megamem
#PBS -l mem=2TB
#PBS -l ncpus=48
#PBS -l storage=gdata/cj50+gdata/jk72+gdata/ik11+gdata/hh5+gdata/rt52+gdata/gb6
#PBS -M daniel.atwater@utas.edu.au
source ~/.bashrc
purge_and_load_conda
python ./regrid_bran_for_cice6.py
