#!/bin/bash
#PBS -v PROJECT
module use /g/data/hh5/public/modules
module load conda/analysis3
conda activate afim
python ./regrid_standalone_forcing_cice6.py
