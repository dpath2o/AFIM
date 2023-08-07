#!/bin/bash
#PBS -P jk72
#PBS -l walltime=48:00:00
#PBS -q hugemem
#PBS -l mem=1024GB
#PBS -l ncpus=16
#PBS -l storage=gdata/cj50+gdata/jk72+gdata/ik11+gdata/hh5+gdata/rt52+gdata/gb6+gdata/qv56
#PBS -M daniel.atwater@utas.edu.au
module load nco
cd /scratch/jk72/da1339/cice-dirs/input/AFIM/forcing/0p25/JRA55/regridded_raw/
ncrcat -v huss,prra,prsn,rlds,rsds,tas,uas,vas,psl *2005.nc ../jra55do_v1p5_merged_2005.nc
ncrename -v huss,spchmd -v prra,ttlpcp -v rlds,dlwsfc -v tsds,glbrad -v tas,airtmp -v uas,wnd