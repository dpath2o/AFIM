#!/bin/bash
#PBS -P jk72
#PBS -l walltime=6:00:00
#PBS -q hugemem
#PBS -l mem=512GB
#PBS -l ncpus=8
#PBS -l storage=gdata/cj50+gdata/jk72+gdata/ik11+gdata/hh5+gdata/rt52+gdata/gb6+gdata/qv56
#PBS -M daniel.atwater@utas.edu.au
module purge
module use /g/data/hh5/public/modules
module load conda/analysis3-22.01
python3 ~/src/python/merge_regridded_jra55.py