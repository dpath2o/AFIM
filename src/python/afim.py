'''
This is a python module for functions and objects (classes) that are beneficial
to the Antarctica Fast Ice Modelling that I'm doing for my PhD.

author: dpath2o, daniel.atwater@utas.edu.au, May 2022 
'''

import os
import glob
import pdb
import subprocess
import json
import pygmt
import numpy             as np
import xarray            as xr
import xesmf             as xe
import pandas            as pd
import metpy.calc        as mpc
import matplotlib.pyplot as plt
import cartopy.crs       as ccrs
from dask.distributed    import Client

############################################################################
# globals
ncrcat = os.path.join('/','apps','nco','5.0.5','bin','ncrcat')

########################################################################################################################################################
############################################################# NON-CLASS FUNCTIONS ######################################################################
########################################################################################################################################################
def compute_sfc_qsat(d2m, sp):
    '''
    compute specific humidity at 2-metres based on dewpoint and surface pressure
    '''
    Rdry = 287.0597
    Rvap = 461.5250
    a1   = 611.21
    a3   = 17.502
    a4   = 32.19
    T0   = 273.16
    E    = a1 * np.exp(a3 * (d2m-T0) / (d2m-a4) )
    return (Rdry/Rvap) * E / (sp - ( (1-Rdry/Rvap) * E) )

############################################################################
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
    Cpst0  = c0 + c1*T_mld + c2*T_mld**2 + c3*T_mld**3 + c4*T_mld**4 + (a0 + a1*T_mld + a2*T_mld**2)*S_mld + (b0 + b1*T_mld + b2*T_mld**2)*S_mld*np.sqrt(S_mld)
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
    dCp0t0 = (a0 + a1*T_mld + a2*T_mld**2 + a3*T_mld**3 + a4*T_mld**4)*P_mld + (b0 + b1*T_mld + b2*T_mld**2 + b3*T_mld**3 + b4*T_mld**4)*P_mld**2 + (c0 + c1*T_mld + c2*T_mld**2 + c3*T_mld**3)*P_mld**3
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
    dCpstp = ((d0 + d1*T_mld + d2*T_mld**2 + d3*T_mld**3 + d4*T_mld**4)*S_mld + (e0 + e1*T_mld + e2*T_mld**2)*S3_2)*P_mld + ((f0 + f1*T_mld + f2*T_mld**2 + f3*T_mld**3)*S_mld + g0*S3_2)*P_mld**2 + ((h0 + h1*T_mld + h2*T_mld**2)*S_mld + j1*T_mld*S3_2)*P_mld**3
    cp     = Cpst0 + dCp0t0 + dCpstp
    # HEAT FLUX
    return cp + (sst-T_mld)

############################################################################
def compute_diffuse_sfc_em(em_sfc, em_toa):
    '''
    a simple function to compute the diffuse short/long wave EM at surface
    '''
    k_t = em_sfc / em_toa
    return 0.952 - 1.041 * np.exp(np.exp((2.3 - 4.702*k_t)))

############################################################################
def compute_u2_from_u10(u10):
    '''
    a simple function transfer/compute wind speed at 2-metres from reported wind speed at 10-metres
    '''
    return (u10 * 4.87) / np.log((67.8 * 10) - 5.42)

############################################################################
def compute_sfc_airrho(t2m, d2m, sp):
    '''
    compute atmospheric air density at the surface based on temperature, dewpoint and surface pressure using metpy calculation
    '''
    RH  = mpc.relative_humidity_from_dewpoint(t2m,d2m)
    r   = mpc.mixing_ratio_from_relative_humidity(sp, t2m, RH)
    return mpc.density(sp, t2m, r)

############################################################################
def compute_sfc_qsat(d2m, sp):
    '''
    compute specific humidity at 2-metres based on dewpoint and surface pressure
    '''
    Rdry = 287.0597
    Rvap = 461.5250
    a1   = 611.21
    a3   = 17.502
    a4   = 32.19
    T0   = 273.16
    E    = a1 * np.exp(a3 * (d2m-T0) / (d2m-a4) )
    return (Rdry/Rvap) * E / (sp - ( (1-Rdry/Rvap) * E) )

############################################################################
def ocean_temp_from_theta(p0, p):
    '''
    compute ocean temperature from potential temperature
    p0 :: reference pressure [db]
    p  :: pressure at level [db]

    '''
    dp = p0 - p

########################################################################################################################################################
######################################################################### CLASSES ######################################################################
########################################################################################################################################################

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
        CICE grid resolution:: this is string that defines the grid resolution of CICE and is strictly 
                            :: used in naming of files and directories
                            :: default is "0p1"
        Interpolation method:: "regrid_interp_type"
        for re-gridding     :: this is a string that is passed to ESMF_RegridWeightGen for defining the
                            :: type of re-gridding method used for interpolation. See 
                            :: ESMF_RegridWeightGen for more details on methods of interpolation. This 
                            :: is used in the method cice_prep.esmf_generate_weights(), however, that
                            :: method does have an over-ride option 'interp_type' that when provided
                            :: will over-ride the setting in this JSON file.
                            :: default is "bilinear"
        Directories:
           The following inputs are useful directories
           base data        :: "D_data"
                            :: defines the base directory where data will be output
           BRAN data        :: "D_BRAN", location of native BRAN climatology data
           ERA5 data        :: "D_ERA5", location of native ERA5 climatology data
           Model Input      :: "D_modin", probably no longer necessary as it was used in previous 
                            :: version of this module
           Graphical        :: "D_graph", this is figures/plots/animations are stored
           ACCESS-OM2 output:: "D_access-om2_out", this defines which ACCESS-OM2 cycle data will be 
                            :: in the concatenation routines and then subsequently forcing the CICE
                            :: stand-alone model
        Files:
           The followin are specific files that are used by the module
           Bathymetry       :: "F_amo2bath", this file is no presently used in any of the methods,
                            :: but is highly relevant and might want to be used for plotting
           CICE grid        :: "F_G_CICE", this contains the spherical grid (lat/lon in decimal degrees) 
                            :: and is an AFIM-compatible-version of the grid that stand-alone CICE uses
                            :: when running. It is used here in the  generation of re-gridding weights and 
                            :: re-gridding of climatology datasets
           BRAN grid        :: "F_G_BRAN", this contains the spherical grid (lat/lon in decimal degrees) 
                            :: of the BRAN climatology and is used in the generation of re-gridding weights
                            :: and re-gridding of climatology datasets
           ERA5 grid        :: "F_G_ERA5", this contains the spherical grid (lat/lon in decimal degrees) 
                            :: of the ERA5 climatology and is used in the generation of re-gridding weights
                            :: and re-gridding of climatology datasets
           BRAN weights     :: "F_BRAN_weights", once created via the AFIM python module
                            :: cice_prep.esmf_generate_weights this is the name of the file used in the 
                            :: method cice_prep.regrid_climatology() for BRAN climatology
           ERA5 weights     :: "F_ERA5_weights", once created via the AFIM python module
                            :: cice_prep.esmf_generate_weights this is the name of the file used in the 
                            :: method cice_prep.regrid_climatology() for ERA5 climatology
           
        '''
        # read in the JSON file
        self.FJSON = FJSON
        with open(FJSON) as f:
            PARAMS = json.load(f)
        # miscellaneous
        self.CICE_ver           = PARAMS['CICE_ver']
        self.start_date         = PARAMS['start_date']
        self.n_months           = PARAMS['n_months']
        self.n_years            = PARAMS['n_years']
        self.regrid_interp_type = PARAMS['regrid_interp_type']
        # grid
        self.G_res        = PARAMS['G_res']
        # directories
        self.D_01         = PARAMS['D01']
        self.D_data       = os.path.join(PARAMS['D01'],'data')
        self.D_grids      = os.path.join(PARAMS['D01'],'grids')
        self.D_weights    = os.path.join(self.D_grids,'weights')
        self.D_BRAN       = PARAMS['D_BRAN']
        self.D_ERA5       = PARAMS['D_ERA5']
        self.D_CICE       = PARAMS['D_CICE']
        self.D_graph      = PARAMS['D_graph']
        self.D_modin      = PARAMS['D_modin']
        self.D_ac_om2_out = PARAMS['D_access-om2_out']
        self.D_CICE_IC    = os.path.join(PARAMS['D_CICE'],'input','AFIM','ic',self.G_res)
        self.D_CICE_force = os.path.join(PARAMS['D_CICE'],'input','AFIM','forcing',self.G_res)
        self.D_CICE_grid  = os.path.join(PARAMS['D_CICE'],'input','AFIM','grid',self.G_res)
        # files
        self.F_gx1ic    = os.path.join(self.D_CICE_grid,'gx1','iced_gx1_v6.2005-01-01.nc')
        self.F_gx1bath  = os.path.join(self.D_CICE_grid,'gx1','global_gx1.bathy.nc')
        self.F_aom2bath = PARAMS['F_aom2bath']
        self.F_G_CICE   = PARAMS['F_G_CICE']
        self.F_G_BRAN   = PARAMS['F_G_BRAN']
        self.F_G_ERA5   = PARAMS['F_G_ERA5']
        self.F_ERA5_weights = PARAMS['F_ERA5_weights']
        self.F_BRAN_weights = PARAMS['F_BRAN_weights']
        self.F_ERA5_reG_form = PARAMS['F_ERA5_reG_form']
        self.F_BRAN_reG_form = PARAMS['F_BRAN_reG_form']

    ############################################################################################
    def esmf_generate_weights(self, ds_name='ERA5', weight_file_prefix='access-om2_cice',
                              interp_type='', regrid_options='-p none'):
        '''
        '''
        if not interp_type:
            interp_type = self.regrid_interp_type
        if ds_name=='ERA5':
            P_wgt = os.path.join( self.D_weights, 'map_{from_name:s}_{to_name:s}_{grid_res:s}_{type:s}.nc'.format(from_name=ds_name,
                                                                                                                  to_name=weight_file_prefix,
                                                                                                                  grid_res=self.G_res,
                                                                                                                  type=interp_type))
            str_esmf_gen_weights = 'ESMF_RegridWeightGen -m {type:s} {opts:s} -s {src:s} -d {dst:s} -w {wgt:s}'.format(type=interp_type,
                                                                                                              opts=regrid_options,
                                                                                                              src=self.F_G_ERA5,
                                                                                                              dst=self.F_G_CICE,
                                                                                                              wgt=P_wgt)
        elif ds_name=='BRAN':
            P_wgt = os.path.join( self.D_weights, 'map_{from_name:s}_{to_name:s}_{grid_res:s}_{type:s}.nc'.format(from_name=ds_name,
                                                                                                                  to_name=weight_file_prefix,
                                                                                                                  grid_res=self.G_res,
                                                                                                                  type=interp_type))
            str_esmf_gen_weights = 'ESMF_RegridWeightGen -m {type:s} {opts:s} -s {src:s} -d {dst:s} -w {wgt:s}'.format(type=interp_type,
                                                                                                              opts=regrid_options,
                                                                                                              src=self.F_G_BRAN,
                                                                                                              dst=self.F_G_CICE,
                                                                                                              wgt=P_wgt)
        print("generating weights with the command:\n",str_esmf_gen_weights)
        os.system(str_esmf_gen_weights)

    ############################################################################################
    def regrid_climatology(self,ds_name='ERA5'):
        '''
        After defining the JSON fields attempt to regrid a dataset using CDO and format that regridding
        so that it is directly usable by CICE6. An example for regridding ERA5 and BRAN data for a test
        of fast ice on the Mawson Coast is simple as:
        import afim
        cice_prep = afim.cice_prep("./fastice_mawson_coast.json")
        cice_prep.regrid_climatology(ds_name="ERA5")
        cice_prep.regrid_climatology(ds_name="BRAN")
        '''
        # load in the grid file to grid to
        G_CICE                       = xr.open_dataset(self.F_G_CICE,chunks={'ny':100,'nx':100})
        G_CICE                       = G_CICE.drop(('clon_t','clat_t','clon_u','clat_u','ulat','ulon','htn','hte','angle','angleT','tarea','uarea'))
        ln_tmp                       = G_CICE.tlon.values[0,:]
        lt_tmp                       = G_CICE.tlat.values[:,0]
        G_CICE                       = G_CICE.drop(('tlon','tlat'))
        G_CICE['lon']                = ln_tmp
        G_CICE['lat']                = lt_tmp
        G_CICE['lon']                = G_CICE.lon * 180/np.pi
        G_CICE['lat']                = G_CICE.lat * 180/np.pi
        G_CICE['lon'].attrs['units'] = 'degrees_east'
        G_CICE['lat'].attrs['units'] = 'degrees_north'
        # time arrays
        stdts = pd.date_range(self.start_date, freq='MS', periods=self.n_months*self.n_years)
        spdts = pd.date_range(self.start_date, freq='M', periods=self.n_months*self.n_years)
        # Directories and files
        if ds_name=='ERA5':
            for i in range(self.n_years):
                yr_str   = stdts[i].strftime('%Y')
                stt_str  = stdts[i].strftime('%Y%m%d')
                stp_str  = spdts[i].strftime('%Y%m%d')
                F_ERA5 = '{var:s}_era5_oper_sfc_{stt:s}-{stp:s}.nc'.format(var='2t',stt=stt_str,stp=stp_str)
                P_ERA5 = os.path.join(self.D_ERA5,
                                      self.ERA5_dir_var_list[0],
                                      yr_str,
                                      F_ERA5)
                #G_ERA5      = xr.open_dataset(P_ERA5).rename({'longitude': 'lon', 'latitude': 'lat'})
                G_ERA5      = xr.open_dataset(F_G_ERA5).rename({'longitude': 'lon', 'latitude': 'lat'})
                rg          = xe.Regridder(G_ERA5,G_CICE,method="bilinear",periodic=True,filename=F_ERA5_weights,reuse_weights=True)
                t2m_na      = xr.open_mfdataset(os.path.join(self.D_ERA5,'2t',yr_str,'*.nc'), parallel=True, chunks={"time": 1})
                t2m         = rg(t2m_na)
                msdwlwrf_na = xr.open_mfdataset(os.path.join(self.D_ERA5,'msdwlwrf',yr_str,'*.nc'), parallel=True, chunks={"time": 1})
                msdwlwrf    = rg(msdwlwrf_na)
                msdwswrf_na = xr.open_mfdataset(os.path.join(self.D_ERA5,'msdwswrf',yr_str,'*.nc'), parallel=True, chunks={"time": 1})
                msdwswrf    = rg(msdwswrf_na)
                mtpr_na     = xr.open_mfdataset(os.path.join(self.D_ERA5,'mtpr',yr_str,'*.nc'), parallel=True, chunks={"time": 1})
                mtpr        = rg(mtpr_na)
                u10_na      = xr.open_mfdataset(os.path.join(self.D_ERA5,'u10',yr_str,'*.nc'), parallel=True, chunks={"time": 1})
                u10         = rg(u10_na)
                v10_na      = xr.open_mfdataset(os.path.join(self.D_ERA5,'v10',yr_str,'*.nc'), parallel=True, chunks={"time": 1})
                v10         = rg(v10_na)
                d2m_na      = xr.open_mfdataset(os.path.join(self.D_ERA5,'2d',yr_str,'*.nc'), parallel=True, chunks={"time": 1})
                d2m         = rg(d2m_na)
                sp_na       = xr.open_mfdataset(os.path.join(self.D_ERA5,'sp',yr_str,'*.nc'), parallel=True, chunks={"time": 1})
                sp          = rg(d2m_na)
                qsat        = compute_sfc_qsat(d2m.d2m, sp.sp)
                ATM         = xr.Dataset({ "airtmp" : t2m.t2m,
                                           "dlwsfc" : msdwlwrf.msdwlwrf,
                                           "glbrad" : msdwswrf.msdwswrf,
                                           "spchmd" : qsat,
                                           "ttlpcp" : mtpr.mtpr,
                                           "wndewd" : u10.u10,
                                           "wndnwd" : v10.v10 },
                                         coords = { "LON"  : (["ny","nx"],LON_rg),
                                                    "LAT"  : (["ny","nx"],LAT_rg),
                                                    "time" : time.values })
                ATM         = ATM.rename_dims({"ny":"nj","nx":"ni"})
                P_ATM       = os.path.join(self.D_frcg,self.F_ERA5_reG_form.format(yr=yr_str))
                ATM.to_netcdf(P_ATM)
        elif ds_name=='BRAN':
            for i in range(self.n_years):
                G_BRAN  = xr.open_dataset(F_G_BRAN).rename({'xt_ocean': 'lon', 'yt_ocean': 'lat'}).chunk({'lon':100,'lat':100})
                G_BRAN  = G_BRAN.drop(('st_edges_ocean','tmask','umask','kmu','kmt','hu','ht','xu_ocean','yu_ocean'))
                rg      = xe.Regridder(G_BRAN,G_CICE,method="bilinear",periodic=True,filename=F_bran_t_wgt,reuse_weights=True)
                temp_na = xr.open_mfdataset(os.path.join(D_bran,'ocean_temp_{yr:s}*.nc'.format(yr=yr_str)), parallel=True, chunks={"Time":1,'st_ocean':1}).rename({'xt_ocean': 'lon', 'yt_ocean': 'lat'})
                temp_na = temp_na.drop(('st_edges_ocean','nv','Time_bounds','average_DT','average_T2','average_T1'))
                temp    = rg(temp_na)
                salt_na = xr.open_mfdataset(os.path.join(D_bran,'ocean_salt_{yr:s}*.nc'.format(yr=yr_str)), parallel=True, chunks={"Time":1,'st_ocean':1}).rename({'xt_ocean': 'lon', 'yt_ocean': 'lat'})
                salt_na = salt_na.drop(('Time_bounds','average_DT','average_T2','average_T1','nv'))
                salt    = rg(salt_na)
                mld_na  = xr.open_mfdataset(os.path.join(D_bran,'ocean_mld_{yr:s}*.nc'.format(yr=yr_str)), parallel=True, chunks={'Time':1}).rename({'xt_ocean': 'lon', 'yt_ocean': 'lat'})
                mld_na  = mld_na.drop(('Time_bounds','average_T1','average_T2','average_DT','nv'))
                mld     = rg(mld_na)
                uocn_na = xr.open_mfdataset(os.path.join(D_bran,'ocean_u_{yr:s}*.nc'.format(yr=yr_str)), parallel=True, chunks={"Time":1,'st_ocean':1}).rename({'xu_ocean': 'lon', 'yu_ocean': 'lat'})
                uocn_na = uocn_na.drop(('st_edges_ocean','Time_bounds','average_T1','average_T2','average_DT','nv'))
                uocn    = rg(uocn_na)
                vocn_na = xr.open_mfdataset(os.path.join(D_bran,'ocean_v_{yr:s}*.nc'.format(yr=yr_str)), parallel=True, chunks={"Time":1,'st_ocean':1}).rename({'xu_ocean': 'lon', 'yu_ocean': 'lat'})
                vocn_na = vocn_na.drop(('st_edges_ocean','Time_bounds','average_T1','average_T2','average_DT','nv'))
                vocn    = rg(vocn_na)
                eta_na  = xr.open_mfdataset(os.path.join(D_bran,'ocean_eta_t_{yr:s}*.nc'.format(yr=yr_str)), parallel=True, chunks={"Time":1}).rename({'xt_ocean': 'lon', 'yt_ocean': 'lat'})
                eta_na  = eta_na.drop(('Time_bounds','average_T1','average_T2','average_DT','nv'))
                eta     = rg(eta_na)
                sst     = temp.sel(st_ocean=0,method='nearest').values
                sss     = salt.sel(st_ocean=0,method='nearest').values
                u       = uocn.sel(st_ocean=0,method="nearest").values
                v       = vocn.sel(st_ocean=0,method="nearest").values
                dhdx    = compfute_ocn_sfc_slope(eta, G_CICE, direction='x', grid_scale_factor=100)
                dhdy    = compfute_ocn_sfc_slope(eta, G_CICE, direction='y', grid_scale_factor=100)
                qdp     = np.zeros( ( len(temp.Time), len(temp.lat), len(temp.lon) ) )
                for j in range(len(temp.Time)):
                    qdp[j,:,:] = compute_ocn_heat_flux( temp.isel(Time=j), salt.isel(Time=j), mld.isel(Time=j), lat_name='lat')
                ln    = sst.lon
                lt    = sst.lat
                LN,LT = np.meshgrid(ln,lt)
                time  = sst.Time
                OCN   = xr.Dataset({ "T"    : sst,
                                     "S"    : sss,
                                     "hblt" : mld.mld.values,
                                     "u"    : u,
                                     "v"    : v,
                                     "dhdx" : dhdx,
                                     "dhdy" : dhdy,
                                     "qdp"  : qdp},
                                   coords = { "LON"  : (["ny","nx"],LN),
                                              "LAT"  : (["ny","nx"],LT),
                                              "time" : time.values })
                OCN   = OCN.rename_dims({"ny":"nj","nx":"ni"})
                P_ATM = os.path.join(self.D_frcg,self.F_BRAN_reG_form.format(yr=yr_str))
                OCN.to_netcdf(P_ocn)

    ############################################################################################
    def concat_access_om2(self,start_year='2005',stop_year='2018'):
        '''
        '''
        import re
        import fnmatch
        F_2d_prefixes_to_concat = ['ocean-2d-surface_temp-1-daily-mean-ym_',
                                   'ocean-2d-surface_salt-1-daily-mean-ym_',
                                   'ocean-2d-mld-1-daily-mean-ym_']
        F_3d_prefixes_to_concat = ['ocean-3d-u-1-daily-mean-ym_',
                                   'ocean-3d-v-1-daily-mean-ym_']
        search_pattern = '*{yr:s}_01.nc'.format(yr=start_year)
        matches        = []
        for root, dirnames, filenames in os.walk(self.D_ac_om2_out):
            for basename in filenames:
                if fnmatch.fnmatch(basename,search_pattern):
                    matches.append(os.path.join(root, basename))
        D_out0_parsed = matches[0].split(os.sep)
        D_out0        = D_out0_parsed[8]
        D_out0_n      = int(re.split('(\d+)',D_out0)[1])
        D_outN_n      = ( (int(stop_year) - int(start_year))*4 ) + D_out0_n
        yr_out        = int(start_year)
        for i in np.arange(D_out0_n,D_outN_n,4):
            D_cats = [os.path.join(self.D_ac_om2_out,'output{:d}'.format(i),'ocean'),
                      os.path.join(self.D_ac_om2_out,'output{:d}'.format(i+1),'ocean'),
                      os.path.join(self.D_ac_om2_out,'output{:d}'.format(i+2),'ocean'),
                      os.path.join(self.D_ac_om2_out,'output{:d}'.format(i+3),'ocean')]
            for j in F_2d_prefixes_to_concat:
                P_wilds = [os.path.join(D_cats[0],'{:s}*.nc'.format(j)),
                           os.path.join(D_cats[1],'{:s}*.nc'.format(j)),
                           os.path.join(D_cats[2],'{:s}*.nc'.format(j)),
                           os.path.join(D_cats[3],'{:s}*.nc'.format(j))]
                P_out = os.path.join(self.D_data,'ac-om2-{var:s}{yr:s}.nc'.format(var=j,yr=str(yr_out)))
                str_concat = 'ncrcat {:s} {:s} {:s} {:s} {:s}'.format(P_wilds[0],P_wilds[1],P_wilds[2],P_wilds[3],P_out)
                if os.exists(P_out):
                    print('File exists: {:s}\n SKIPPING concatentation')
                else:
                    print("concatenating: ",str_concat)
                    os.system(str_concat)
            for j in F_3d_prefixes_to_concat:
                P_wilds = [os.path.join(D_cats[0],'{:s}*.nc'.format(j)),
                           os.path.join(D_cats[1],'{:s}*.nc'.format(j)),
                           os.path.join(D_cats[2],'{:s}*.nc'.format(j)),
                           os.path.join(D_cats[3],'{:s}*.nc'.format(j))]
                P_out = os.path.join(self.D_data,'ac-om2-{var:s}{yr:s}.nc'.format(var=j,yr=str(yr_out)))
                str_concat = 'ncrcat -d st_ocean,0,1 {:s} {:s} {:s} {:s} {:s}'.format(P_wilds[0],P_wilds[1],P_wilds[2],P_wilds[3],P_out)
                if os.exists(P_out):
                    print('File exists: {:s}\n SKIPPING concatentation')
                else:
                    print("concatenating: ",str_concat)
                    os.system(str_concat)
            yr_out+=1
