#!/bin/bash
#PBS -l ncpus=10
#PBS -l walltime=24:00:00
#PBS -l mem=100GB
#PBS -v PROJECT
module load nco
Dnative=/scratch/jk72/da1339/ERA5/native
Dregrid=/scratch/jk72/da1339/ERA5/regrid
tvars=('2d' '2t' 'mror' 'msdrswrf' 'msdwlwrf' 'msdwswrf' 'msl' 'msr' 'msror' 'mtdwswrf' 'mtnlwrf' 'mtnswrf' 'mtpr' 'sp' 'z')
uvars=('10u' '10v' 'metss' 'mntss')
years=`seq 2005 2021`
Mt=${HOME}/grids/0p1/map_ERA5_to_access-om2-3700x3600_nco.20220909.nc
Mu=${HOME}/grids/0p1/map_ERA5_to_access-om2-3700x3600_nco_ugrid.20220909.nc
# t_vars
for year in ${years[@]}; do
	for var in ${tvars[@]}; do
		for Fnc in $Dnative/${var}/${year}/*.nc; do 
		    Fin=`basename ${Fnc}`
		    mkdir -p $Dregrid/${year}/${var}
		    if [[ ! -f $Dregrid/${year}/${var}/${Fin} ]]; then
			echo "re-gridding: $Dregrid/${year}/${var}/${Fin}"
			ncremap -m ${Mt} $Dnative/${var}/${year}/${Fin} $Dregrid/${year}/${var}/${Fin}
		    fi
		done
	done
done
# u_vars
for year in ${years[@]}; do
        for var in ${uvars[@]}; do
              	for Fnc in $Dnative/${var}/${year}/*.nc; do
                    Fin=`basename ${Fnc}`
		    mkdir -p $Dregrid/${year}/${var}
                    if [[ ! -f $Dregrid/${year}/${var}/${Fin} ]]; then
			echo "re-gridding: $Dregrid/${year}/${var}/${Fin}"
                        ncremap -m ${Mu} $Dnative/${var}/${year}/${Fin} $Dregrid/${year}/${var}/${Fin}
                    fi
		done
       	done
done
