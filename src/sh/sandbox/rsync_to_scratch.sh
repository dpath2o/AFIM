#!/bin/bash
#PBS -P jk72
#PBS -l walltime=48:00:00
#PBS -q hugemembw
#PBS -l mem=1020GB
#PBS -l ncpus=28
#PBS -l storage=gdata/jk72
#PBS -M daniel.atwater@utas.edu.au
rsync -az /scratch/jk72/da1339/cice-dirs/runs/afim_025_jra55do_aom2/ /g/data/jk72/da1339/afim_output/20230708_jra55_aom2_ndte_3710_mld_ice_dynamic/
rsync -az /scratch/jk72/da1339/cice-dirs/runs/afim_025_jra55do_bran_ncar/ /g/data/jk72/da1339/afim_output/20230708_jra55_bran_ndte_3710_mld_ice_dynamic/
