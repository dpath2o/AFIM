#!/bin/bash
#PBS -P jk72
#PBS -l walltime=12:00:00
#PBS -q megamembw
#PBS -l mem=3TB
#PBS -l ncpus=64
#PBS -l storage=gdata/cj50+gdata/jk72+gdata/ik11+gdata/hh5+gdata/rt52+gdata/gb6
#PBS -M daniel.atwater@utas.edu.au
module purge
module load pbs intel-compiler openmpi netcdf pnetcdf
ulimit -s unlimited
cd ${HOME}/src/CICE/afim_0p1_mpi/
pwd 
./cice.run
