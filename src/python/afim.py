'''
This is a python module for functions and objects (classes) that are beneficial
to the Antarctica Fast Ice Modelling that I'm doing for my PhD.

author: dpath2o, daniel.atwater@utas.edu.au, May 2022 
'''

import os
import subprocess
import json
import numpy       as np
import xarray      as xr
import pandas      as pd
import metpy.calc  as mpc
import cartopy.crs as ccrs
from cartopy       import feature as cfeature

############################################################################
# globals
cdo = os.path.join('/','Users','dpath2o','opt','anaconda3','envs','afim','bin','cdo')

############################################################################
# DIRECTORIES
class cice_prep:
    '''
    '''
    def __init__(self,FJSON='',**kwargs):
        '''
        '''
        # read in the JSON file
        self.FJSON = FJSON
        with open(FJSON) as f:
            PARAMS = json.load(f)
        # miscellaneous
        self.CICE_ver     = PARAMS['CICE_ver']
        self.start_date   = PARAMS['start_date']
        self.n_months     = PARAMS['n_months']
        self.n_years      = PARAMS['n_years']
        self.cdo_reG_opts = PARAMS['cdo_reG_opts']
        # these are switches for conditional computations
        self.switches = PARAMS['switches']
        # grid
        self.G_res = PARAMS['G_res']
        # directories
        self.D_01    = PARAMS['D01']
        self.D_02    = PARAMS['D02']
        self.D_03    = PARAMS['D03']
        self.D_BRAN  = PARAMS['D_BRAN']
        self.D_JRA55 = PARAMS['D_JRA55']
        self.D_MERRA = PARAMS['D_MERRA']
        self.D_CAWCR = PARAMS['D_CAWCR']
        self.D_ERA5  = PARAMS['D_ERA5']
        self.D_SOSE  = PARAMS['D_SOSE']
        self.D_IC    = PARAMS['D_IC']
        self.D_frcg  = os.path.join(PARAMS['D_frcg'],self.G_res)
        self.D_grid  = PARAMS['D_grid']
        self.D_graph = PARAMS['D_graph']
        self.D_modin = PARAMS['D_modin']
        self.D_reG   = os.path.join(self.D_modin,'grids',self.G_res,'regrid')
        # files
        self.F_gx1ic   = PARAMS['F_gx1ic']
        self.F_gx1bath = PARAMS['F_gx1bath']
        self.F_Gt      = os.path.join(self.D_modin,'grids',self.G_res,'g{:s}_cice_tgrid.nc'.format(self.G_res))
        self.F_Gu      = os.path.join(self.D_modin,'grids',self.G_res,'g{:s}_cice_ugrid.nc'.format(self.G_res))
        # directory names that are associated with variables
        self.BRAN_dir_var_list  = PARAMS['BRAN_dir_var_list']
        self.JRA55_dir_var_list = PARAMS['JRA55_dir_var_list']
        self.MERRA_dir_var_list = PARAMS['MERRA_dir_var_list']
        self.ERA5_dir_var_list  = PARAMS['ERA5_dir_var_list']
        # variables
        self.ERA5_frcg_tvars     = PARAMS['ERA5_frcg_tvars']
        self.ERA5_frcg_uvars     = PARAMS['ERA5_frcg_uvars']
        self.BRAN_frcg_tvars     = PARAMS['BRAN_frcg_tvars']
        self.BRAN_frcg_uvars     = PARAMS['BRAN_frcg_uvars']
        self.ERA5_alt_var_names  = PARAMS['ERA5_alt_var_names']
        self.ERA5toCICE_varnames = PARAMS['ERA5toCICE_varnames']
        self.ERA5_filename_form  = PARAMS['ERA5_filename_form']
        self.CICE6_atm_var_names = PARAMS['CICE6_atm_var_names']
        self.ERA5_derv_rho_name  = PARAMS['ERA5_derv_rho_name']
        self.ERA5_derv_qsat_name = PARAMS['ERA5_derv_qsat_name']
        self.ERA5_yrmo_Ftmp_form = PARAMS['ERA5_yrmo_Ftmp_form']
        self.ERA5_F_force_form   = PARAMS['ERA5_F_force_form']

    ############################################################################################
    def regrid_wrapper(self,ds_name='ERA5'):
        '''
        '''
        # Directories and files
        if ds_name=='ERA5':
            D_ERA5 = self.D_ERA5
            D_reG  = os.path.join(self.D_reG,ds_name)
            D_frcg = os.path.join(self.D_frcg,'hourly')
        F_Gt  = self.F_Gt
        F_Gu  = self.F_Gu
        G_res = self.G_res
        # variable lists
        list_of_var_names = self.ERA5_frcg_tvars + self.ERA5_frcg_uvars
        # time arrays
        stdts = pd.date_range(self.start_date, freq='MS', periods=self.n_months*self.n_years)
        spdts = pd.date_range(self.start_date, freq='M', periods=self.n_months*self.n_years)
        for j in np.arange(0,len(stdts),1):
            yr_str   = stdts[j].strftime('%Y')
            yrmo_str = stdts[j].strftime('%Y%m')
            stt_str  = stdts[j].strftime('%Y%m%d')
            stp_str  = spdts[j].strftime('%Y%m%d')
            if self.switches["compute_sfc_qsat"]:
                proc  = 0
                N_sp  = 0
                N_d2m = 0
            for var_name in list_of_var_names:
                if var_name=='sp': N_sp += 1
                if var_name=='d2m': N_d2m += 1
                Fn_reG = self.ERA5_filename_form.format(field=var_name,stt=stt_str,stp=stp_str)
                Pn_reG = os.path.join(D_ERA5,var_name,yr_str,Fn_reG)
                Do_reG = os.path.join(D_reG,yr_str,var_name)
                if not(os.path.exists(Do_reG)): os.makedirs(Do_reG)
                Fo_reG = self.ERA5_yrmo_Ftmp_form.format(G_res=G_res,field=var_name,yrmo=yrmo_str)
                Po_reG = os.path.join(Do_reG,Fo_reG)
                if self.switches['regridding']:
                    Prename = ''
                    if var_name in self.ERA5toCICE_varnames:
                        Frename = self.ERA5_yrmo_Ftmp_form.format(G_res=G_res,field=self.ERA5toCICE_varnames[var_name],yrmo=yrmo_str)
                        Prename = os.path.join(Do_reG,Frename)
                    if not(os.path.exists(Prename)):
                        if not(os.path.exists(Po_reG)):
                            print("\nAttempting to regrid\nFROM: {:s}\nTO: {:s}".format(Pn_reG,Po_reG))
                            if var_name in self.ERA5_frcg_tvars:
                                regrid_nc(Pn_reG, Po_reG, F_Gt, cdo_options=self.cdo_reG_opts)
                            elif var_name in self.ERA5_frcg_uvars:
                                regrid_nc(Pn_reG, Po_reG, F_Gu, cdo_options=self.cdo_reG_opts)
                    if (var_name in self.ERA5toCICE_varnames) and not(os.path.exists(Prename)):
                        print("Renaming internal variable \nFROM: {:s}\nTO: {:s}".format(var_name,self.ERA5toCICE_varnames[var_name]))
                        rename_nc_var(Po_reG,Prename,self.ERA5_alt_var_names[var_name],self.ERA5toCICE_varnames[var_name])
                if N_sp==1 and N_d2m==1 and proc==0:
                    print('\nAttempting to compute SPECIFIC HUMIDITY (or "qsat")')
                    self.compute_derivations(ds_name=ds_name,compute='qsat',yr_str=yr_str,yrmo_str=yrmo_str)
                    proc=1
            if spdts[j].is_year_end:
                # merge months into one year for each field
                for i in self.ERA5toCICE_var_names:
                    Fin_merge  = os.path.join(D_reG,yr_str,i,'*.nc')
                    Fout_merge = self.ERA5_yr_Ftmp_form.format(G_res=G_res,field=i,yr=yr_str)
                    Pout_merge = os.path.join(D_reG,yr_str,Fout_merge)
                    merge_nc(Fin_merge,Pout_merge,self.switches['delete_merge_input'],merge_along_time_axis=True)
                # merge each field into one NetCDF for the year as per CICE forcing file requirements
                Fin_merge  = os.path.join(D_reG,yr_str,'*.nc')
                Fout_merge = self.ERA5_F_force_form.format(G_res=G_res,yr=yr_str)
                Pout_merge = os.path.join(D_reG,Fout_merge)
                merge_nc(Fin_merge,Pout_merge,self.switches['delete_merge_input'])

    #########################################################################################
    def compute_derivations(self,ds_name='ERA5',compute='qsat',yr_str='',yrmo_str=''):
        '''
        '''
        #stdts = pd.date_range(self.start_date, freq='MS', periods=self.n_months*self.n_years)
        #spdts = pd.date_range(self.start_date, freq='M', periods=self.n_months*self.n_years)
        D_reG = os.path.join(self.D_reG,ds_name)
        G_res = self.G_res
        if compute=='airrho':
            var_names = ['d2m','sp','t2m']
            vout      = 'rho'
            vnew      = self.ERA5_derv_rho_name
        elif compute=='qsat':
            var_names = ['d2m','sp']
            vout      = 'qsat'
            vnew      = self.ERA5toCICE_varnames['qsat']
        Dnew = os.path.join(D_reG,yr_str,vnew)
        Fnew = self.ERA5_yrmo_Ftmp_form.format(G_res=G_res,field=vnew,yrmo=yrmo_str)
        Pnew = os.path.join(Dnew,Fnew)
        if os.path.exists(Pnew):
            print('Derived file already exists -- skipping derivation: {:s}'.format(Pnew))
            return
        #for j in np.arange(0,len(stdts),1):
        #    yrmo_str = stdts[j].strftime('%Y%m')
        D1 = os.path.join(D_reG,yr_str,var_names[0])
        F1 = self.ERA5_yrmo_Ftmp_form.format(G_res=G_res,field=var_names[0],yrmo=yrmo_str)
        print("LOADING: {:s}".format(os.path.join(D1,F1)))
        ds1 = xr.open_dataset(os.path.join(D1,F1))
        D2   = os.path.join(D_reG,yr_str,var_names[1])
        F2   = self.ERA5_yrmo_Ftmp_form.format(G_res=G_res,field=var_names[1],yrmo=yrmo_str)
        print("LOADING: {:s}".format(os.path.join(D2,F2)))
        ds2  = xr.open_dataset(os.path.join(D2,F2))
        if compute=='airrho':
            D3 = os.path.join(D_reG,yr_str,var_names[2])
            F3 = self.ERA5_yrmo_Ftmp_form.format(G_res=G_res,field=var_names[2],yrmo=yrmo_str)
            print("LOADING: {:s}".format(os.path.join(D3,F3)))
            ds3 = xr.open_dataset(os.path.join(D3,F3))
            print("COMPUTING AIR DENSITY at the surface")
            ds = compute_sfc_airrho(ds1.d2m, ds2.sp, ds3.t2m)
        elif compute=='qsat':
            print("COMPUTING SPECIFIC HUMIDITY at the surface")
            ds = compute_sfc_qsat(ds1.d2m, ds2.sp)
        Dout = os.path.join(D_reG,yr_str,vout)
        Fout = self.ERA5_yrmo_Ftmp_form.format(G_res=G_res,field=vout,yrmo=yrmo_str)
        Pout = os.path.join(Dout,Fout)
        print("SAVING to NetCDF: {:s}".format(Pout))
        ds.to_netcdf(Pout)
        print("RENAMING variable {:s} to {:s} and saving to NetCDF: {:s}".format('__xarray_dataarray_variable__',vnew,Pnew))
        rename_nc_var(Pout,Pnew,'__xarray_dataarray_variable__',vnew)

############################################################################
# a simple function to call subprocess commands
def call_subprocess(sys_str_in):
    '''
    Given a string that represents a valid command line input attempt to run
    that command via subprocess module call. The string and subprocess call
    return is reported to STDOUT
    '''
    print('System Call: ',sys_str_in)
    sys_str_out = subprocess.call(sys_str_in,shell=True)
    print('System call return: ',sys_str_out)

############################################################################
# a simple function to regrid a NetCDF file using CDO
def regrid_nc(Fin, Fout, Fgrd, Vname='', tslice='',
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
    print('DATASET BEFORE REGRIDDING: ', ds.sizes)

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
                                                                                type=cdo_regrid_type,
                                                                                opts=cdo_options,
                                                                                Fgrd=Fgrd,
                                                                                Fin=Fin,
                                                                                Fout=Fout)
    else:
        SYS_CALL = '{cdo:s} {type:s},{Fgrd:s} {Fin:s} {Fout:s}'.format(cdo=cdo,
                                                                       type=cdo_regrid_type,
                                                                       Fgrd=Fgrd,
                                                                       Fin=Fin,
                                                                       Fout=Fout)
    call_subprocess(SYS_CALL)
 # show what the file looks like after re-gridding
    ds = xr.open_dataset(Fout, engine='netcdf4')
    print('DATASET AFTER REGRIDDING: ', ds.sizes)

############################################################################
# a simple function to merge NetCDF along the time axis using CDO
def merge_nc(Fin,Fout,del_in,merge_along_time_axis=False):
    '''
    Given a wildcard path string, merge those files into a single NetCDF file
    along the time axis using CDO. Optionally, choose to delete/remove input
    files after being merged. 
    '''
    if merge_along_time_axis:
        SYS_CALL = '{cdo:s} mergetime {Fin:s} {Fout:s}'.format(cdo=cdo,Fin=Fin,Fout=Fout)
    else:
        SYS_CALL = '{cdo:s} merge {Fin:s} {Fout:s}'.format(cdo=cdo,Fin=Fin,Fout=Fout)
    call_subprocess(SYS_CALL)
    if del_in:
        SYS_CALL = 'rm -f {Fin:s}'.format(Fin=Fin)
        call_subprocess(SYS_CALL)

############################################################################
# a simple function to rename NetCDF variable names using CDO
def rename_nc_var(Fin,Fout,var_in,var_out):
    '''
    '''
    SYS_CALL = '{cdo:s} chname,{vin:s},{vout:s} {Fin:s} {Fout:s}'.format(cdo=cdo,
                                                                         vin=var_in,
                                                                         vout=var_out,
                                                                         Fin=Fin,
                                                                         Fout=Fout)
    call_subprocess(SYS_CALL)
    # remove/delete input file
    SYS_CALL = 'rm -f {Fin:s}'.format(Fin=Fin)
    call_subprocess(SYS_CALL)

############################################################################
# a simple function to report on NetCDF variables
def report_nc_var_names(Fin, mmm=False):
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

############################################################################
# a simple function to compute the diffuse short/long wave EM at surface
def compute_diffuse_sfc_em(em_sfc, em_toa):
    '''
    '''
    k_t = em_sfc / em_toa
    return 0.952 - 1.041 * np.exp(np.exp((2.3 - 4.702*k_t)))

############################################################################
# a simple function transfer/compute wind speed at 2-metres from reported
# wind speed at 10-metres
def compute_u2_from_u10(u10):
    '''
    '''
    return (u10 * 4.87) / np.log((67.8 * 10) - 5.42)

############################################################################
# compute atmospheric air density at the surface based on temperature,
# dewpoint and surface pressure using metpy calculation
def compute_sfc_airrho(t2m, d2m, sp):
    '''
    '''
    RH  = mpc.relative_humidity_from_dewpoint(t2m,d2m)
    r   = mpc.mixing_ratio_from_relative_humidity(sp, t2m, RH)
    return mpc.density(sp, t2m, r)

############################################################################
# compute specific humidity at 2-metres based on dewpoint and surface
# pressure
def compute_sfc_qsat(d2m, sp):
    '''
    '''
    Rdry = 287.0597
    Rvap = 461.5250
    a1   = 611.21
    a3   = 17.502
    a4   = 32.19
    T0   = 273.16
    E    = a1 * np.exp(a3 * (d2m-T0) / (d2m-a4) )
    return (Rdry/Rvap) * E / (sp - ( (1-Rdry/Rvap) * E) )

#############################################################################
# compute ocean temperature from potential temperature
def ocean_temp_from_theta(p0, p):
    '''
    p0 :: reference pressure [db]
    p  :: pressure at level [db]

    '''
    dp = p0 - p
