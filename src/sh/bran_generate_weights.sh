#!/bin/bash
#PBS -P jk72
#PBS -l walltime=8:00:00
#PBS -l ncpus=384
#PBS -l mem=1TB
#PBS -l storage=gdata/cj50+gdata/jk72+gdata/ik11+gdata/hh5+gdata/rt52+gdata/gb6
#PBS -M daniel.atwater@utas.edu.au
module use /g/data/hh5/public/modules
module load conda/analysis3
module load nco
source ~/.bashrc
conda activate afim
python src/python/generate_weights_bran.py

