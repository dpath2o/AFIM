#!/bin/bash
#PBS -N era5_parser
#PBS -P jk72
#PBS -l walltime=48:00:00
#PBS -q normalbw
#PBS -l mem=32GB
#PBS -l ncpus=28
#PBS -l storage=gdata/cj50+gdata/jk72+gdata/ik11+gdata/hh5+gdata/rt52+gdata/gb6
#PBS -M daniel.atwater@utas.edu.au
module purge
module use /g/data/hh5/public/modules
module load nco
module load openmpi
module load conda/analysis3
ulimit -s unlimited

var_Fnames=(2t 2d 100u 100v 10u 10v z cin sst cape sp msl blh tcc lcc mcc hcc crr lsrr csfr lssfr iews inss ishf msror mssror mser msmr mlspf mlspr mcpr msr 10fg)

for i in "${var_Fnames[@]}"; do
    for j in {1959..2022}; do
        Dname=$(printf "/g/data/rt52/era5/single-levels/reanalysis/%s/%04d/%s_*" $i $j $i);
        mpirun -np 28 python ./AFIM/src/python/era5_parser.py $i $Dname
    done;
done