#!/bin/bash
#PBS -P jk72
#PBS -l walltime=48:00:00
#PBS -q normalbw
#PBS -l mem=8GB
#PBS -l ncpus=2
#PBS -l storage=gdata/jk72
#PBS -M daniel.atwater@utas.edu.au
rsync -az /g/data/jk72/da1339/afim_input/ /scratch/jk72/da1339/afim_input/

