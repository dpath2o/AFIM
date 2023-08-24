#!/bin/bash
#PBS -P jk72
#PBS -l walltime=48:00:00
#PBS -q normal
#PBS -l mem=32GB
#PBS -l ncpus=4
#PBS -l storage=gdata/jk72
#PBS -M daniel.atwater@utas.edu.au
rsync -az /scratch/jk72/da1339/cice-dirs/runs/afim_025_jra55do_aom2/ /g/data/jk72/da1339/afim_output/20230814_jra55_aom2_3yr_GIneg0p2m/

