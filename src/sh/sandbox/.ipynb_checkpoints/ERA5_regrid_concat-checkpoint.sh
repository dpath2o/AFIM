#!/bin/bash
#PBS -P jk72
#PBS -l walltime=36:00:00
#PBS -q megamem
#PBS -l mem=2TB
#PBS -l ncpus=48
#PBS -l storage=gdata/cj50+gdata/jk72+gdata/ik11+gdata/hh5+gdata/rt52+gdata/gb6
#PBS -M daniel.atwater@utas.edu.au
module use /g/data/hh5/public/modules
module load conda/analysis3
module load nco
source ~/.bashrc
#ncap_opts="-s 'airtmp=float(airtmp)' -s 'dlwsfc=float(dlwsfc)' -s 'glbrad=float(glbrad)' -s 'spchmd=float(spchmd)' -s 'ttlpcp=float(ttlpcp)' -s 'wndewd=float(wndewd)' -s 'wndnwd=float(wndnwd)'"
ncrcat_opts="--cnk_dmn time,1"
D_data_mo="/g/data/jk72/da1339/data/ERA5/monthly"
D_data_yr="/g/data/jk72/da1339/cice-dirs/input/AFIM/forcing/0p1/8XDAILY/"
year=2006
#cd $D_data_mo
#for i in `ls -h ./double`; do
#    ncap2 $ncap_opts $D_data_tmp/double/${i} $D_data_tmp/single/${i}
#done
#for i in {1..12}; do
#    F_i=$(printf "JRA55_03hr_forcing_%s%02d" $year $i)
#    F_o=$(printf "JRA55_03hr_forcing_%s_%02d.nc" $year $i)
#    ncrcat $ncrcat_opts $D_data_tmp/single/${F_i} $D_data_yr/${F_o}
#done
ncrcat $ncrcat_opts $D_data_mo/JRA55_03hr_forcing_${year}_* $D_data_yr/JRA55_03hr_forcing_${year}.nc
