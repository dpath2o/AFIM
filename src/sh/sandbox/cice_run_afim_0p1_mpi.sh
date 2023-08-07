#!/bin/bash
#PBS -P jk72
#PBS -l walltime=6:00:00
#PBS -q megamem
#PBS -l mem=3TB
#PBS -l ncpus=96
#PBS -l storage=gdata/cj50+gdata/jk72+gdata/ik11+gdata/hh5+gdata/rt52+gdata/gb6
#PBS -M daniel.atwater@utas.edu.au
module purge
module load intel-compiler/2021.6.0 netcdf/4.9.0 openmpi/4.1.4 pnetcdf
ulimit -s unlimited
cd ${HOME}/src/CICE/tx1_vanilla
F_log=/home/581/da1339/logs/cice_tx1_vanilla.log
./cice.run >&! $F_log
