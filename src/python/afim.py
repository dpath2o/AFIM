'''
This is a python module for functions and objects (classes) that are beneficial
to the Antarctica Fast Ice Modelling that I'm doing for my PhD.

author: dpath2o, daniel.atwater@utas.edu.au, May 2022 
'''

import os
import subprocess
import xarray as xr

############################################################################
# a simple function to regrid a NetCDF file using CDO
def regrid_netcdf_with_cdo(Fin, Fout, Fgrd, Vname='', tslice='',
                           cdo_regrid_type='remapbic', cdo_options=''):
    '''
    This function takes a gridded NetCDF and attempts to regrid it using 2D
    NetCDF grid file. The function does not return any values. The function
    requires the following inputs:

    'Fin'  :: the full path name of the NetCDF file to be regridded
    'Fout' :: the full path name of the *new* (re-gridded) NetCDF file
    'Fgrd' :: the full path name of the 2D NetCDF grid file that will be used
              to preform the re-gridding

    Additional/Optional inputs:
    'Vname'  :: the short name of the particular variable in 'Fin' which will be
                re-gridded. All other variables in the file will be ignored.
    'tslice' :: if the user is only concerned with a particular timeframe of
                the variable then this is place to put this in string format
                YYYY-MM-DD[:YYYY-MM-DD]
    'cdo_regrid_type' :: type of regridding interpolation to perform. See CDO
                         manual for more detailed help.
    'cdo_options'     :: specific options passed to CDO. See CDO manual for
                         more detailed help.
    '''

    cdo  = os.path.join('/','Users','dpath2o','opt','anaconda3','envs','afim','bin','cdo')
    Ftmp = os.path.join('/', 'Users', 'dpath2o', 'tmp', 'tmp_pre-regrid.nc')

    # OPEN NETCDF: parse variable name and time
    ds = xr.open_dataset(Fin, engine='netcdf4')
    if tslice:
        if Vname:
            ds = ds[Vname].sel(time=tslice, method='nearest')
        else:
            ds = ds.sel(time=tslice, method='nearest')
    else:
        if Vname:
            ds = ds[Vname]

    # show what the dataset looks like before re-gridding
    print('DATASET BEFORE REGRIDDING: ', ds)

    # TEMPORARY NETCDF FILE: create a temporary NetCDF for CDO to read from
    if os.path.exists(Ftmp):
        os.remove(Ftmp)
    ds.to_netcdf(Ftmp)

    # IF OUTPUT FILE ALREADY EXISTS CLOBBER IT
    if os.path.exists(Fout):
        os.remove(Fout)

    # REGRID
    if cdo_options:
        SYS_CALL = '{cdo:s} {opts:s} {type:s},{Fgrd:s} {Fin:s} {Fout:s}'.format(cdo=cdo,
                                                                                type=cdo_regrid_type
                                                                                opts=cod_options,
                                                                                Fgrd=Fgrd,
                                                                                Fin=Ftmp,
                                                                                Fout=Fout)
    else:
        SYS_CALL = '{cdo:s} {type:s},{Fgrd:s} {Fin:s} {Fout:s}'.format(cdo=cdo,
                                                                       type=cdo_regrid_type
                                                                       Fgrd=Fgrd,
                                                                       Fin=Ftmp,
                                                                       Fout=Fout)
    print('System call: ', SYS_CALL)
    SYS_RET = subprocess.call(SYS_CALL, shell=True)
    print('System call return: ', SYS_RET)

    # show what the file looks like after re-gridding
    ds = xr.open_dataset(Fout, engine='netcdf4')
    print('DATASET AFTER REGRIDDING: ', ds)

############################################################################
# a simple function to report on NetCDF variables
def report_netcdf_variables(Fin, mmm=False):
    '''
    Given a full path name to a NetCDF file then report on the variables' names
    and attributes. Optional switch for the function is to enable reporting on
    mean, minimum and maximum of each variable. Note, that if the NetCDF
    variable is large in time and space then this mean, min, max switch can
    take a long time to complete. The function does not return any values.
    '''

    # OPEN NETCDF: parse variable name and time
    ds = xr.open_dataset(Fin, engine='netcdf4')

    # show what the dataset looks like before re-gridding
    print('DATASET via xarray print: ', ds)

    for i in ds.keys():
        print('')
        print('variable name       : ', i)
        print('variable attributes : ', ds[i].attrs)
        if mmm:
            print('variable min        : ', ds[i].min())
            print('variable max        : ', ds[i].max())
            print('variable mean       : ', ds[i].mean())
