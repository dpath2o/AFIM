'''
This is a python module for functions and objects (classes) that are beneficial
to the Antarctica Fast Ice Modelling that I'm doing for my PhD.

author: dpath2o, daniel.atwater@utas.edu.au, May 2022 
'''

import os
import pdb
import subprocess
import json
import pygmt
import numpy       as np
import xarray      as xr
import pandas      as pd
import metpy.calc  as mpc

############################################################################
# globals
cdo      = os.path.join('/','Users','dpath2o','opt','anaconda3','envs','afim','bin','cdo')
ncks     = os.path.join('/','opt','homebrew','bin','ncks')
ncwa     = os.path.join('/','opt','homebrew','bin','ncwa')
ncrename = os.path.join('/','opt','homebrew','bin','ncrename')
ncdump   = os.path.join('/','opt','homebrew','bin','ncdump')

############################################################################
# DIRECTORIES
class cice_prep:
    '''
    This class is written to prepare datasets for CICE6 input, either forcing files or initial conditions. It
    works from a supplied JSON file for defining the parameters required.

    External software dependencies:
    - CDO
    *** With this dependency, the user must edit this module's global variable "cdo" (just a few lines above these words),
        and provide the path location to the binary cdo executables.

    Internal software dependencies:
    *** Author recommends using package manager "Anaconda" to download, install and manage the following packages. Version
        dependencies with the following packages are not being rigorously accounted for by the author of this module and
        hence errors and crashes of internal functions and classes of this module may occur do to differences in versions
        in which the author wrote the package and versions that the current user is using.
    - xarray
    - metpy
    - pygmt
    - numpy
    - pandas
    - json
    
    '''
    def __init__(self,FJSON='',**kwargs):
        '''
        Below is a description of the required/expected fields of the JSON file and the contents of those fields.
        Initialising this class in essence gives the user access to those fields and sub-modules that rely
        heavily on the definitions.

        Given a JSON file define the following:

        CICE version        :: 'CICE_ver'
                            :: version number of CICE. Currently on CICE6 is supported
                            :: default is 6 
        Start date          :: 'start_date'
                            :: the date that is used to determine when to start looking for data
                            :: default is 2010-01-01
        Number of Months    :: 'n_months'
                            :: the number of dates from the start date with which to consider biulding
                            :: a period within a year. For instance if one was just considering June, July
                            :: and August, then they would have a start date of 2010-06-01 and set this
                            :: value to 3
                            :: defualt is 12
        Number of Years     :: much like the number of months, but this defines the duration
                            :: default is 2
        CDO Re-grid Options :: these are the optional parameters that are passed to CDO during regridding.
                            :: Consult CDO documentation for more informaiton.
                            :: default is "-f nc -b F64"
        CDO Re-grid Type    :: this the type of regridding to perform
                            :: Consult CDO documentation for more informaiton.
                            :: default is "remapbic"
        '''
        # read in the JSON file
        self.FJSON = FJSON
        with open(FJSON) as f:
            PARAMS = json.load(f)
        # miscellaneous
        self.CICE_ver         = PARAMS['CICE_ver']
        self.start_date       = PARAMS['start_date']
        self.n_months         = PARAMS['n_months']
        self.n_years          = PARAMS['n_years']
        self.cdo_reG_opts     = PARAMS['cdo_reG_opts']
        self.cdo_reG_type     = PARAMS['cdo_reG_type']
        self.qdp_depth        = PARAMS['qdp_depth']
        self.BRAN_sfc_lev_idx = PARAMS['BRAN_sfc_lev_idx'] 
        # these are switches for conditional computations
        self.switches = PARAMS['switches']
        # grid
        self.G_res = PARAMS['G_res']
        self.G_scale_fact = PARAMS['G_scale_fact']
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
        self.D_CICE  = PARAMS['D_CICE']
        self.D_graph = PARAMS['D_graph']
        self.D_modin = PARAMS['D_modin']
        self.D_IC    = os.path.join(PARAMS['D_CICE'],'input','CICE_data','ic',self.G_res)
        self.D_frcg  = os.path.join(PARAMS['D_CICE'],'input','CICE_data','forcing',self.G_res)
        self.D_grid  = os.path.join(PARAMS['D_CICE'],'input','CICE_data','grid',self.G_res)
        self.D_reG   = os.path.join(self.D_modin,'regridded',self.G_res)
        # files
        self.F_gx1ic    = os.path.join(self.D_grid,'gx1','iced_gx1_v6.2005-01-01.nc')
        self.F_gx1bath  = os.path.join(self.D_grid,'gx1','global_gx1.bathy.nc')
        self.F_aom2bath = PARAMS['F_aom2bath']
        self.F_Gt       = os.path.join(self.D_modin,'grids',self.G_res,'g{:s}_cice_tgrid.nc'.format(self.G_res))
        self.F_Gu       = os.path.join(self.D_modin,'grids',self.G_res,'g{:s}_cice_ugrid.nc'.format(self.G_res))
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
        self.BRAN_3D_var_names   = PARAMS['BRAN_3D_var_names']
        self.ERA5_rename_vars  = PARAMS['ERA5_rename_vars']
        self.BRAN_rename_vars  = PARAMS['BRAN_rename_vars']
        self.ERA5toCICE_varnames = PARAMS['ERA5toCICE_varnames']
        self.BRANtoCICE_varnames = PARAMS['BRANtoCICE_varnames']
        self.ERA5_filename_form  = PARAMS['ERA5_filename_form']
        self.BRAN_filename_form  = PARAMS['BRAN_filename_form']
        self.ERA5_F_force_form   = PARAMS['ERA5_F_force_form']
        self.BRAN_F_force_form   = PARAMS['BRAN_F_force_form']
        self.Ftmp_yrmo_form      = PARAMS['Ftmp_yrmo_form']
        self.ERA5_yr_Ftmp_form   = PARAMS['ERA5_yr_Ftmp_form']
        self.BRAN_yr_Ftmp_form   = PARAMS['BRAN_yr_Ftmp_form']
        self.CICE6_atm_var_names = PARAMS['CICE6_atm_var_names']
        self.ERA5_derv_rho_name  = PARAMS['ERA5_derv_rho_name']
        self.ERA5_derv_qsat_name = PARAMS['ERA5_derv_qsat_name']        

    ############################################################################################
    def regrid_wrapper(self,ds_name='ERA5'):
        '''
        After defining the JSON fields attempt to regrid a dataset using CDO and format that regridding
        so that it is directly usable by CICE6. An example for regridding ERA5 and BRAN data for a test
        of fast ice on the Mawson Coast is simple as:
        import afim
        cice_prep = afim.cice_prep("./fastice_mawson_coast.json")
        cice_prep.regrid_wrapper(ds_name="ERA5")
        cice_prep.regrid_warpper(ds_name="BRAN")
        
        '''
        # Directories and files
        if ds_name=='ERA5':
            D_ERA5 = self.D_ERA5
            D_reG  = os.path.join(self.D_reG,ds_name)
            D_frcg = os.path.join(self.D_frcg,'hourly')
            list_of_var_names = self.ERA5_frcg_tvars + self.ERA5_frcg_uvars
        elif ds_name=='BRAN':
            D_BRAN = self.D_BRAN
            D_reG  = os.path.join(self.D_reG,ds_name)
            D_frcg = os.path.join(self.D_frcg,'daily')
            list_of_var_names = self.BRAN_frcg_tvars + self.BRAN_frcg_uvars
        F_Gt  = self.F_Gt
        F_Gu  = self.F_Gu
        G_res = self.G_res
        G     = xr.open_dataset(self.F_Gt).hte
        # time arrays
        stdts = pd.date_range(self.start_date, freq='MS', periods=self.n_months*self.n_years)
        spdts = pd.date_range(self.start_date, freq='M', periods=self.n_months*self.n_years)
        for j in np.arange(0,len(stdts),1):
            # create temporal related strings to reduce clutter in code
            yr_str   = stdts[j].strftime('%Y')
            mo_str   = stdts[j].strftime('%m')
            yrmo_str = stdts[j].strftime('%Y%m')
            stt_str  = stdts[j].strftime('%Y%m%d')
            stp_str  = spdts[j].strftime('%Y%m%d')
            # create loop switches/counters so that when a value of 1 is reached with any of
            # the following variable in a conditional statement then processing of a particular
            # derivation can proceeed
            N_sp   = 0
            N_2d   = 0
            N_2t   = 0
            N_eta  = 0
            N_temp = 0
            N_salt = 0
            N_mld  = 0
            # create the name of the full-year-all-variables file first because if it exists I'm
            # making the assumption that I don't want it to be re-made
            Pn_mrg_var = os.path.join(D_reG,yr_str,'*.nc')
            if ds_name=='ERA5':
                Fo_mrg_var = self.ERA5_F_force_form.format(G_res=G_res,yr=yr_str)
            elif ds_name=='BRAN':
                Fo_mrg_var = self.BRAN_F_force_form.format(G_res=G_res,yr=yr_str)
            Po_mrg_var = os.path.join(D_reG,Fo_mrg_var)
            # skip this year if the output file aready exists
            if os.path.exists(Po_mrg_var):
                print('Merged full-year-all-variables already exists -- skip *ALL* regridding and merging.\n{:s}'.format(Po_mrg_var))
                continue
            for var_name in list_of_var_names:
                print("\n\nWorking on {vin:s} for the year-month {date:s}".format(vin=var_name,date=yrmo_str))
                # as variables are regridded turn-on their switches for processing derived variables
                if var_name=='sp'    : N_sp   += 1
                if var_name=='2d'    : N_2d   += 1
                if var_name=='2t'    : N_2t   += 1
                if var_name=='eta_t' : N_eta  += 1
                if var_name=='temp'  : N_temp += 1
                if var_name=='salt'  : N_salt += 1
                if var_name=='mld'   : N_mld  += 1
                 # create the temporary path name for regrdding the file before it has it's variable name change
                Do_reG = os.path.join(D_reG,yr_str,var_name)
                if not(os.path.exists(Do_reG)): os.makedirs(Do_reG)
                Fo_reG = self.Ftmp_yrmo_form.format(G_res=G_res,field=var_name,yrmo=yrmo_str)
                Po_reG = os.path.join(Do_reG,Fo_reG)
                # create the input path name depending on the dataset
                # this is data that is from the reanalysis source
                if ds_name=='ERA5':
                    Fn_reG  = self.ERA5_filename_form.format(field=var_name,stt=stt_str,stp=stp_str)
                    Pn_reG  = os.path.join(D_ERA5,var_name,yr_str,Fn_reG)
                    if not(os.path.exists(Po_reG)):
                        if var_name in self.ERA5_frcg_tvars:
                            regrid_ERA5(Pn_reG, Po_reG, F_Gt, Vname=var_name, var_renames=self.ERA5_rename_vars,
                                        cdo_options=self.cdo_reG_opts, cdo_regrid_type=self.cdo_reG_type)
                        else:
                            regrid_ERA5(Pn_reG, Po_reG, F_Gu, Vname=var_name, var_renames=self.ERA5_rename_vars,
                                        cdo_options=self.cdo_reG_opts, cdo_regrid_type=self.cdo_reG_type)
                    if N_sp==1 and N_2d==1:
                        Do = os.path.join(D_reG,yr_str,'qsat')
                        if not(os.path.exists(Do)): os.makedirs(Do)
                        Fo = self.Ftmp_yrmo_form.format(G_res=G_res,field='qsat',yrmo=yrmo_str)
                        Po = os.path.join(Do,Fo)
                        if not(os.path.exists(Po)):
                            print('\nAttempting to compute SPECIFIC HUMIDITY (or "qsat")')
                            D1 = os.path.join(D_reG,yr_str,'2d')
                            F1 = self.Ftmp_yrmo_form.format(G_res=G_res,field='2d',yrmo=yrmo_str)
                            print("LOADING: {:s}".format(os.path.join(D1,F1)))
                            d2m = xr.open_dataset(os.path.join(D1,F1))
                            D2 = os.path.join(D_reG,yr_str,'sp')
                            F2 = self.Ftmp_yrmo_form.format(G_res=G_res,field='sp',yrmo=yrmo_str)
                            print("LOADING: {:s}".format(os.path.join(D2,F2)))
                            sp = xr.open_dataset(os.path.join(D2,F2))
                            qsat = compute_sfc_qsat(d2m.d2m,sp.sp)
                            print("SAVING to NetCDF: {:s}".format(Po))
                            qsat.to_netcdf(Po)
                            SYS_CALL = '{ncrename:s} -O -v __xarray_dataarray_variable__,spchmd {Fout:s} {Fout:s}'.format(ncrename=ncrename,Fout=Po)
                            call_subprocess(SYS_CALL)
                elif ds_name=='BRAN':
                    Fn_reG = self.BRAN_filename_form.format(field=var_name,yr=yr_str,mo=mo_str)
                    Pn_reG = os.path.join(D_BRAN,var_name,Fn_reG)
                    if not(os.path.exists(Po_reG)):
                        if var_name in self.BRAN_frcg_tvars:
                            regrid_BRAN(Pn_reG, Po_reG, F_Gt, var_type=var_name, var_3D_names=self.BRAN_3D_var_names.keys(),
                                        cdo_options=self.cdo_reG_opts, cdo_regrid_type=self.cdo_reG_type, var_renames=self.BRAN_rename_vars)
                        else:
                            regrid_BRAN(Pn_reG, Po_reG, F_Gu, var_type=var_name, var_3D_names=self.BRAN_3D_var_names.keys(),
                                        cdo_options=self.cdo_reG_opts, cdo_regrid_type=self.cdo_reG_type, var_renames=self.BRAN_rename_vars)
                    if N_eta==1:
                        Do = os.path.join(D_reG,yr_str,'dhdx')
                        if not(os.path.exists(Do)): os.makedirs(Do)
                        Fo = self.Ftmp_yrmo_form.format(G_res=G_res,field='dhdx',yrmo=yrmo_str)
                        Po = os.path.join(Do,Fo)
                        if not(os.path.exists(Po)):
                            print('\nAttempting to compute OCEAN SURFACE SLOPE in the X-DIRECTION (or "dhdx")')
                            D1 = os.path.join(D_reG,yr_str,'eta_t')
                            F1 = self.Ftmp_yrmo_form.format(G_res=G_res,field='eta_t',yrmo=yrmo_str)
                            print("LOADING: {:s}".format(os.path.join(D1,F1)))
                            eta = xr.open_dataset(os.path.join(D1,F1))
                            dhdx = compute_ocn_sfc_slope(eta.eta_t,G,direction='x',grid_scale_factor=self.G_scale_fact)
                            print("SAVING to NetCDF: {:s}".format(Po))
                            dhdx.to_netcdf(Po)
                            SYS_CALL = '{ncrename:s} -O -v __xarray_dataarray_variable__,dhdx -v lat,LAT -v lon,LON {Fout:s} {Fout:s}'.format(ncrename=ncrename,Fout=Po)
                            call_subprocess(SYS_CALL)
                        Do = os.path.join(D_reG,yr_str,'dhdy')
                        if not(os.path.exists(Do)): os.makedirs(Do)
                        Fo = self.Ftmp_yrmo_form.format(G_res=G_res,field='dhdy',yrmo=yrmo_str)
                        Po = os.path.join(Do,Fo)
                        if not(os.path.exists(Po)):
                            print('\nAttempting to compute OCEAN SURFACE SLOPE in the X-DIRECTION (or "dhdy")')
                            dhdy = compute_ocn_sfc_slope(eta.eta_t,G,direction='y',grid_scale_factor=self.G_scale_fact)
                            print("SAVING to NetCDF: {:s}".format(Po))
                            dhdy.to_netcdf(Po)
                            SYS_CALL = '{ncrename:s} -O -v __xarray_dataarray_variable__,dhdy -v lat,LAT -v lon,LON {Fout:s} {Fout:s}'.format(ncrename=ncrename,Fout=Po)
                            call_subprocess(SYS_CALL)
                    if N_temp==1 and N_salt==1 and N_mld==1:
                        Do = os.path.join(D_reG,yr_str,'qdp')
                        if not(os.path.exists(Do)): os.makedirs(Do)
                        Fo = self.Ftmp_yrmo_form.format(G_res=G_res,field='qdp',yrmo=yrmo_str)
                        Po = os.path.join(Do,Fo)
                        if not(os.path.exists(Po)):
                            print('\nAttempting to compute OCEAN HEAT FLUX (or "qdp")')
                            Fn = self.BRAN_filename_form.format(field='temp',yr=yr_str,mo=mo_str)
                            Pn = os.path.join(D_BRAN,'temp',Fn)
                            print("LOADING: {:s}".format(Pn))
                            T = xr.open_dataset(Pn)
                            Fn = self.BRAN_filename_form.format(field='salt',yr=yr_str,mo=mo_str)
                            Pn = os.path.join(D_BRAN,'salt',Fn)
                            print("LOADING: {:s}".format(Pn))
                            S = xr.open_dataset(Pn)
                            Fn = self.BRAN_filename_form.format(field='mld',yr=yr_str,mo=mo_str)
                            Pn = os.path.join(D_BRAN,'mld',Fn)
                            print("LOADING: {:s}".format(Pn))
                            MLD = xr.open_dataset(Pn)
                            print("COMPUTING OCEAN HEAT_FLUX")
                            QDP = compute_ocn_heat_flux(T,S,MLD)
                            print("SAVING to NetCDF: {:s}".format(Po))
                            QDP.to_netcdf(Po)
                            SYS_CALL = '{ncks:s} -C -O -x -v st_ocean {Fout:s} {Fout:s}'.format(ncks=ncks,Fout=Po)
                            call_subprocess(SYS_CALL)
                            SYS_CALL = '{ncrename:s} -O -v __xarray_dataarray_variable__,qdp -v lat,LAT -v lon,LON {Fout:s} {Fout:s}'.format(ncrename=ncrename,Fout=Po)
                            call_subprocess(SYS_CALL)
                            Pt = os.path.join(Do,'tmp.nc')
                            SYS_CALL = '{cdo:s} {opts:s} {type:s},{Fgrd:s} {Fin:s} {Fout:s}'.format(cdo=cdo,
                                                                                                    type=self.cdo_reG_type,
                                                                                                    opts=self.cdo_reG_opts,
                                                                                                    Fgrd=F_Gt,
                                                                                                    Fin=Po,
                                                                                                    Fout=Pt)
                            call_subprocess(SYS_CALL)
                            os.rename(Pt,Po)
            # at the end of each year go back into each variable and merge that variable across the months
            # then once merged across the months, merge all the variables into one large yearly file as is
            # required by CICE6 
            if spdts[j].is_year_end:
                if ds_name=='ERA5':
                    tmp_list = self.ERA5toCICE_varnames
                elif ds_name=='BRAN':
                    tmp_list = self.BRANtoCICE_varnames
                for var_name in tmp_list:
                    Pn_mrg_mo = os.path.join(D_reG,yr_str,var_name,'*.nc')
                    if ds_name=='ERA5':
                        Fo_mrg_mo = self.ERA5_yr_Ftmp_form.format(G_res=G_res,field=self.ERA5toCICE_varnames[var_name],yr=yr_str)
                    elif ds_name=='BRAN':
                        Fo_mrg_mo = self.BRAN_yr_Ftmp_form.format(G_res=G_res,field=self.BRANtoCICE_varnames[var_name],yr=yr_str)
                    Po_mrg_mo = os.path.join(D_reG,yr_str,Fo_mrg_mo)
                    if not(os.path.exists(Po_mrg_mo)):
                        print('CREATING full-year single variable -- APPROXIMATE TIME 45MINS.\n{:s}'.format(Po_mrg_mo))
                        merge_nc(Pn_mrg_mo,Po_mrg_mo,self.switches['delete_merge_input'],merge_along_time_axis=True)
                    else:
                        print('Merged full-year single variable already exists -- skip merging.\n{:s}'.format(Po_mrg_mo))
                print('CREATING full-year-all-variables -- APPROXIMATE TIME 4.5HRS.\n{:s}'.format(Po_mrg_var))
                merge_nc(Pn_mrg_var,Po_mrg_var,self.switches['delete_merge_input'])


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
#
def info_nc(Fin,all_info=False):
    '''
    '''
    if all_info:
        SYS_CALL = '{cdo:s} infon {Fin:s}'.format(cdo=cdo,Fin=Fin)
    else:
        SYS_CALL = '{cdo:s} sinfon {Fin:s}'.format(cdo=cdo,Fin=Fin)
    call_subprocess(SYS_CALL)

############################################################################
# a simple function to regrid a NetCDF file using CDO
def regrid_ERA5(Fin, Fout, Fgrd, Vname='', tslice='', var_renames={},
              cdo_regrid_type='remapbic', cdo_options='', del_in=False):
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
    Do   = os.path.split(Fout)[0]
    Fo   = os.path.split(Fout)[1]
    Po   = os.path.join(Do,Fout)
    Ftmp = os.path.join(Do,"tmp.nc")
    # REGRID
    if cdo_options:
        SYS_CALL = '{cdo:s} {opts:s} {type:s},{Fgrd:s} {Fin:s} {Fout:s}'.format(cdo=cdo,
                                                                                type=cdo_regrid_type,
                                                                                opts=cdo_options,
                                                                                Fgrd=Fgrd,
                                                                                Fin=Fin,
                                                                                Fout=Ftmp)
    else:
        SYS_CALL = '{cdo:s} {type:s},{Fgrd:s} {Fin:s} {Fout:s}'.format(cdo=cdo,
                                                                       type=cdo_regrid_type,
                                                                       Fgrd=Fgrd,
                                                                       Fin=Fin,
                                                                       Fout=Ftmp)
    call_subprocess(SYS_CALL)
    os.rename(Ftmp,Fout)
    SYS_CALL = '{ncks:s} -C -O -x -v Time_bnds {Fout:s} {Fout:s}'.format(ncks=ncks,Fin=Fout,Fout=Fout)
    #call_subprocess(SYS_CALL)
    SYS_CALL = '{ncrename:s} -O -v lon,LON -v lat,LAT {Fout:s} {Fout:s}'.format(ncrename=ncrename,Fout=Fout)
    #call_subprocess(SYS_CALL)
    if Vname in var_renames.keys():
        SYS_CALL = '{ncrename:s} -O -v {vin:s},{vout:s} {Fout:s} {Fout:s}'.format(ncrename=ncrename,
                                                                                  vin=Vname,
                                                                                  vout=var_renames[Vname],
                                                                                  Fout=Fout)
        call_subprocess(SYS_CALL)

############################################################################
# a simple function to regrid a NetCDF file using CDO
def regrid_BRAN(Fin, Fout, Fgrd, var_type='', var_3D_names='', var_renames={},
                Vname='', tslice='', cdo_regrid_type='remapbic', cdo_options='',
                del_in=False):
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
    if var_type in var_3D_names:
        SYS_CALL = '{cdo:s} sellevidx,1 {Fin:s} {Fout:s}'.format(cdo=cdo,Fin=Fin,Fout=Fout)
        call_subprocess(SYS_CALL)
        SYS_CALL = '{ncwa:s} -O -a st_ocean {Fout:s} {Fout:s}'.format(ncwa=ncwa,Fout=Fout)
        call_subprocess(SYS_CALL)
        SYS_CALL = '{ncks:s} -C -O -x -v st_ocean {Fout:s} {Fout:s}'.format(ncks=ncks,Fout=Fout)
        call_subprocess(SYS_CALL)
        Fin = Fout
    Do   = os.path.split(Fout)[0]
    Fo   = os.path.split(Fout)[1]
    Po   = os.path.join(Do,Fout)
    Ftmp = os.path.join(Do,"tmp.nc")
    # REGRID
    if cdo_options:
        SYS_CALL = '{cdo:s} {opts:s} {type:s},{Fgrd:s} {Fin:s} {Fout:s}'.format(cdo=cdo,
                                                                                type=cdo_regrid_type,
                                                                                opts=cdo_options,
                                                                                Fgrd=Fgrd,
                                                                                Fin=Fin,
                                                                                Fout=Ftmp)
    else:
        SYS_CALL = '{cdo:s} {type:s},{Fgrd:s} {Fin:s} {Fout:s}'.format(cdo=cdo,
                                                                       type=cdo_regrid_type,
                                                                       Fgrd=Fgrd,
                                                                       Fin=Fin,
                                                                       Fout=Ftmp)
    call_subprocess(SYS_CALL)
    # rename filename
    os.rename(Ftmp,Fout)
    SYS_CALL = '{ncks:s} -C -O -x -v Time_bnds {Fout:s} {Fout:s}'.format(ncks=ncks,Fin=Fout,Fout=Fout)
    call_subprocess(SYS_CALL)
    SYS_CALL = '{ncrename:s} -O -v lon,LON -v lat,LAT {Fout:s} {Fout:s}'.format(ncrename=ncrename,Fout=Fout)
    call_subprocess(SYS_CALL)
    if var_type in var_renames.keys():
        SYS_CALL = '{ncrename:s} -O -v {vin:s},{vout:s} {Fout:s} {Fout:s}'.format(ncrename=ncrename,
                                                                                  vin=var_type,
                                                                                  vout=var_renames[var_type],
                                                                                  Fout=Fout)
        call_subprocess(SYS_CALL)

############################################################################
# a simple function slice NetCDF along a *level* index
def slice_nc(Fin,Fout,level_index=0,ds_name=''):
    '''
    Given a NetCDF with three spatial dimensions, slice along the level at
    a given index.
    '''
    SYS_CALL = '{cdo:s} sellevidx,{idx:d} {Fin:s} {Fout:s}'.format(cdo=cdo,Fin=Fin,Fout=Fout,idx=level_index)
    call_subprocess(SYS_CALL)
    if ds_name=='BRAN':
        SYS_CALL = '{ncwa:s} -O -a st_ocean {Fout:s} {Fout:s}'.format(ncwa=ncwa,Fout=Fout)
        call_subprocess(SYS_CALL)
        SYS_CALL = '{ncks:s} -C -O -x -v st_ocean {Fout:s} {Fout:s}'.format(ncks=ncks,Fout=Fout)
        call_subprocess(SYS_CALL)

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
#
def compute_ocn_heat_flux(T,S,mld,lat_name='xt_ocean'):
    '''
    Given ocean temperature and salinity, compute the deep ocean heat flux. The
    approach is first to accurately compute heat capacity and then use the relationship
    Q = cp + dT ; where T is the temperature across the layer and cp is the
    heat capacity calculated at a depth. The heat capacity comes from the reference:
    Fofonff, P. and Millard, R.C. Jr
    Unesco 1983. Algorithms for computation of fundamental properties of seawater,
    1983. Unesco Tech. Pap. in Mar. Sci., No. 44, 53 pp.
    '''
    sst    = T.sel(st_ocean=0      ,method='nearest').temp
    T_mld  = T.sel(st_ocean=mld.mld,method='nearest').temp
    S_mld  = S.sel(st_ocean=mld.mld,method='nearest').salt
    # depth to pressure, convert lats to radians
    latrad = np.sin(np.abs(T[lat_name])*(np.pi/180)) 
    C1     = 5.92e-3 + latrad**(2*5.25e-3)
    P_mld  = ( (1-C1) - np.sqrt( ( (1-C1)**2 ) - (8.84e-6*mld.mld) ) ) / 4.42e-6
    # to convert db to Bar as used in Unesco routines
    P_mld  = P_mld/10
    # HEAT CAPACITY
    c0     =  4.2174e3
    c1     = -3.720283
    c2     =   .1412855
    c3     = -2.654387e-3
    c4     = 2.093236e-5
    a0     = -7.64357
    a1     =   .1072763
    a2     = -1.38385e-3
    b0     =   .1770383
    b1     = -4.07718e-3
    b2     =  5.148e-5
    Cpst0  = c0 + c1*T_mld + c2*T_mld**2 + c3*T_mld**3 + c4*T_mld**4 +\
             (a0 + a1*T_mld + a2*T_mld**2)*S_mld +\
             (b0 + b1*T_mld + b2*T_mld**2)*S_mld*np.sqrt(S_mld)
    a0     = -4.9592e-1
    a1     =  1.45747e-2
    a2     = -3.13885e-4
    a3     =  2.0357e-6
    a4     =  1.7168e-8
    b0     =  2.4931e-4
    b1     = -1.08645e-5
    b2     =  2.87533e-7
    b3     = -4.0027e-9
    b4     =  2.2956e-11
    c0     = -5.422e-8
    c1     =  2.6380e-9
    c2     = -6.5637e-11
    c3     =  6.136e-13
    dCp0t0 = (a0 + a1*T_mld + a2*T_mld**2 + a3*T_mld**3 + a4*T_mld**4)*P_mld +\
                (b0 + b1*T_mld + b2*T_mld**2 + b3*T_mld**3 + b4*T_mld**4)*P_mld**2 +\
                (c0 + c1*T_mld + c2*T_mld**2 + c3*T_mld**3)*P_mld**3
    d0     =  4.9247e-3
    d1     = -1.28315e-4
    d2     =  9.802e-7
    d3     =  2.5941e-8
    d4     = -2.9179e-10
    e0     = -1.2331e-4
    e1     = -1.517e-6
    e2     =  3.122e-8
    f0     = -2.9558e-6
    f1     =  1.17054e-7
    f2     = -2.3905e-9
    f3     =  1.8448e-11
    g0     =  9.971e-8
    h0     =  5.540e-10
    h1     = -1.7682e-11
    h2     =  3.513e-13
    j1     = -1.4300e-12
    S3_2   = S_mld*np.sqrt(S_mld)
    dCpstp = ((d0 + d1*T_mld + d2*T_mld**2 + d3*T_mld**3 + d4*T_mld**4)*S_mld + (e0 + e1*T_mld + e2*T_mld**2)*S3_2)*P_mld +\
             ((f0 + f1*T_mld + f2*T_mld**2 + f3*T_mld**3)*S_mld + g0*S3_2)*P_mld**2 +\
             ((h0 + h1*T_mld + h2*T_mld**2)*S_mld + j1*T_mld*S3_2)*P_mld**3
    cp     = Cpst0 + dCp0t0 + dCpstp
    # HEAT FLUX
    return cp * (sst-T_mld)

############################################################################
#
def compute_ocn_sfc_slope(eta,grid,direction='x',grid_scale_factor=100):
    '''
    '''
    if direction=='x':
        dhdx                    = eta / (grid/grid_scale_factor)
        dhdx.attrs['units']     = 'meters'
        dhdx.attrs['long_name'] = 'sea surface slope in x-direction'
        return dhdx
    if direction=='y':
        dhdy                    = eta / (grid/grid_scale_factor)
        dhdy.attrs['units']     = 'meters'
        dhdy.attrs['long_name'] = 'sea surface slope in y-direction'
        return dhdy

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

#############################################################################
# quick global plot
def map_global_gridded_field(dat,region='g',projection='W180/15c',tit_str='',D_plt='',):
    '''
    '''
    fig = pygmt.Figure()
    fig.basemap(region=region,projection=projection) # Cyl_stere/30/-20/12c
    pygmt.makecpt(cmap="batlow", series=[dat.min().values, dat.max().values])
    fig.grdimage(grid=dat, region="g", frame="a",transparency=25)
    fig.coast(region="g", land="black",frame=["af", f'WSne+t"{tit_str}"'])
    fig.colorbar(frame='af')
    fig.savefig(os.path.join(D_plt,'map_global_gridded_field_{:s}.png'.format(tit_str)))
