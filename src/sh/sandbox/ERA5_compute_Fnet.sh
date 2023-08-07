#!/bin/bash
#PBS -P jk72
#PBS -l walltime=2:00:00
#PBS -q megamembw
#PBS -l mem=3TB
#PBS -l ncpus=64
#PBS -l storage=gdata/cj50+gdata/jk72+gdata/ik11+gdata/hh5+gdata/rt52+gdata/gb6
#PBS -M daniel.atwater@utas.edu.au
module purge
module load nco
module load esmf
D_ERA5="/g/data/rt52/era5/single-levels/reanalysis"
ncrcat_opts="-4 -L 1 --cnk_csh=1000000000 --cnk_plc=g3d --cnk_dmn=time,168"
D_data="/g/data/jk72/da1339/"
D_out="/g/data/jk72/da1339/data/ERA5/Fnet"
year=2005
vars=(msdwlwrf msdwswrf msshf mslhf)
# first add time record dimension
for v in ${vars[@]}; do
    for F in `ls ${D_ERA5}/${v}/${year}`; do
        mkdir -pv ${D_out}/data/ERA5/${v}/${year}
        if [ ! -f "$D_out/data/ERA5/$v/$year/$F" ]; then
            ncks --mk_rec_dmn time ${D_ERA5}/${v}/${year}/${F} ${D_out}/${v}/${year}/${F}
        fi
    done
    if [ ! -f "$D_out/$v/$v_$year.nc" ]; then
        ncrcat $ncrcat_opts ${D_out}/${v}/${year}/*.nc ${D_out}/${v}/${v}_${year}.nc
    fi
done
if [ ! -f "${D_out}/Fnet_${year}.nc" ]; then
    ncks -h -A ${D_out}/msdwlwrf/msdwlwrf_${year}.nc ${D_out}/msdwswrf/msdwswrf_${year}.nc ${D_out}/msshf/msshf_${year}.nc ${D_out}/mslhf/mslhf_${year}.nc ${D_out}/Fnet_${year}.nc
fi 
#ESMF_regrid -m bilinear -s ${D_out}/Fnet_${year}.nc -d ${D_data}/grids/CICE/g0p1_cice_tgrid.nc
#ncrcat $ncrcat_opts $D_data_mo/JRA55_03hr_forcing_${year}_* $D_data_yr/JRA55_03hr_forcing_${year}.nc
