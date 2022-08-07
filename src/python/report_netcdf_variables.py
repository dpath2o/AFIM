#!/usr/env python

# LIBRARIES
import os
import argparse
import subprocess
import xarray as xr

# INPUT ARGUMENTS
parser = argparse.ArgumentParser(prog='regrid_via_cdo',
                                 description='Use CDO to regrid netcdf')
parser.add_argument("-i",
                    help="i [FULL/PATH/TO/ORIGINAL/NetCDF.nc]",
                    required=True,
                    type=str,
                    dest='netcdf_in')
args  = parser.parse_args()
Fin   = args.netcdf_in

# OPEN NETCDF: parse variable name and time
ds = xr.open_dataset(Fin, engine='netcdf4')

# show what the dataset looks like before re-gridding
print('DATASET via xarray print: ', ds)

for i in ds.keys():
    print('')
    print('variable name       : ', i)
    print('variable attributes : ', ds[i].attrs)
    #print('variable min        : ', ds[i].min())
    #print('variable max        : ', ds[i].max())
    #print('variable mean       : ', ds[i].mean())
