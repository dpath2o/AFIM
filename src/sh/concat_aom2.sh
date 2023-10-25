#!/bin/bash
#PBS -P jk72
#PBS -l walltime=36:00:00
#PBS -q hugemem
#PBS -l mem=1024GB
#PBS -l ncpus=16
#PBS -l storage=gdata/cj50+gdata/jk72+gdata/ik11+gdata/hh5+gdata/rt52+gdata/gb6+gdata/qv56
#PBS -M daniel.atwater@utas.edu.au
module load nco
ncrcat -O /scratch/jk72/da1339/cice-dirs/input/AFIM/forcing/0p25/ocean/daily/aom2_2*.nc /scratch/jk72/da1339/cice-dirs/input/AFIM/forcing/0p25/ocean/daily/aom2_as_hycom_ocn_frcg_cice6_0p25_2005_2018.nc