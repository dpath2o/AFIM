#!/bin/bash
#PBS -l ncpus=10
#PBS -l walltime=10:00:00
#PBS -l mem=100mb
#PBS -v PROJECT

module load cdo
ncremap -m ~/grids/0p1/map_ERA5_to_access-om2-3700x3600_nco.20220909.nc -I /g/data/rt52/era5/single-levels/reanalysis/msr/2005 -O /scratch/jk72/da1339
