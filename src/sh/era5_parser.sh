#!/bin/bash
#PBS -P jk72
#PBS -l walltime=12:00:00
#PBS -q normal
#PBS -l mem=32GB
#PBS -l ncpus=2
#PBS -l storage=gdata/cj50+gdata/jk72+gdata/ik11+gdata/hh5+gdata/rt52+gdata/gb6
#PBS -M daniel.atwater@utas.edu.au
module load nco 
ulimit -s unlimited

var_Fnames=(2t 2d 100u 100v 10u 10v z cin sst cape sp msl blh tcc lcc mcc hcc crr lsrr csfr lssfr iews inss ishf msror mssror mser msmr mlspf mlspr mcpr msr 10fg)

cd /g/data/jk72/da1339/era5_sfc/

for i in "${var_Fnames[@]}"; do 
    mkdir $i;
    for j in {1959..2022}; do
	Dname=$(printf "/g/data/rt52/era5/single-levels/reanalysis/%s/%04d/%s_*" $i $j $i);
	for k in $Dname; do
	    Fname=$(echo "${k##*/}");
	    if [[ ! -f "$i/$Fname" ]]; then
		echo $i/$Fname;
		ncks -O -d longitude,1120,1399,1 -d latitude,158,360,1 $k $i/$Fname;
	    fi
	done;
    done; 
done
