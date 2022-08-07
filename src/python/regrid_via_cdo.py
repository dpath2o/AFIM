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
parser.add_argument("-o",
                    help="o [FULL/PATH/TO/REGRIDDED/NetCDF.nc]",
                    required=True,
                    type=str,
                    dest='netcdf_out')
parser.add_argument("-g",
                    help="g [FULL/PATH/TO/GRID/FILE/TO/REGRID/FROM/NetCDF.nc]",
                    required=True,
                    type=str,
                    dest='grid_file')
parser.add_argument("-v",
                    help="v [NAME OF THE VARIABLE IN THE NETCDF TO REGRID]",
                    required=False,
                    type=str,
                    dest='variable_name')
parser.add_argument("-t",
                    help="t [YYYY-MM-DD]",
                    required=False,
                    type=str,
                    dest='slice_time')
args  = parser.parse_args()
Fin   = args.netcdf_in
Fout  = args.netcdf_out
Fgrd  = args.grid_file
Vname = args.variable_name
tslce = args.slice_time

# OPEN NETCDF: parse variable name and time
ds = xr.open_dataset(Fin, engine='netcdf4')
if tslce:
    if Vname:
        ds = ds[Vname].sel(time=tslce, method='nearest')
    else:
        ds = ds.sel(time=tslce, method='nearest')
else:
    if Vname:
        ds = ds[Vname]
    else:
        ds = ds

# show what the dataset looks like before re-gridding
print('DATASET BEFORE REGRIDDING: ', ds)

# TEMPORARY NETCDF FILE: create a temporary NetCDF for CDO to read from
Ftmp = os.path.join('/', 'Users', 'dpath2o', 'tmp', 'tmp_pre-regrid.nc')
if os.path.exists(Ftmp):
    os.remove(Ftmp)
ds.to_netcdf(Ftmp)

# IF OUTPUT FILE ALREADY EXISTS CLOBBER IT
if os.path.exists(Fout):
    os.remove(Fout)

# REGRID
SYS_CALL = '/Users/dpath2o/opt/anaconda3/envs/afim/bin/cdo remapbic,{Fgrd:s} {Fin:s} {Fout:s}'.format(Fgrd=Fgrd, Fin=Ftmp, Fout=Fout)
print('System call: ', SYS_CALL)
SYS_RET = subprocess.call(SYS_CALL, shell=True)
print('System call return: ', SYS_RET)

# show what the file looks like after re-gridding
ds = xr.open_dataset(Fout, engine='netcdf4')
print('DATASET AFTER REGRIDDING: ', ds)
