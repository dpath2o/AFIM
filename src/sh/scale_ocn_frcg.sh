#!/bin/bash
#PBS -P jk72
#PBS -l walltime=12:00:00
#PBS -q hugemembw
#PBS -l mem=1020GB
#PBS -l ncpus=28
#PBS -l storage=gdata/cj50+gdata/jk72+gdata/ik11+gdata/hh5+gdata/rt52+gdata/gb6
#PBS -M daniel.atwater@utas.edu.au
module load nco 
ulimit -s unlimited
cd /scratch/jk72/da1339/cice-dirs/input/AFIM/forcing/0p1/daily
ncflint -C -O -v qdp -w const_fact,0.025 acom2_ocn_frcg_cice6_0p1_2005.nc acom2_ocn_frcg_cice6_0p1_2005.nc acom2_ocn_frcg_cice6_0p1_2005_scaled_qdp.nc
ncflint -C -O -v qdp -w const_fact,0.025 bran_ocn_frcg_cice6_0p1_2005_revised.nc bran_ocn_frcg_cice6_0p1_2005_revised.nc bran_ocn_frcg_cice6_0p1_2005_scaled_qdp.nc
