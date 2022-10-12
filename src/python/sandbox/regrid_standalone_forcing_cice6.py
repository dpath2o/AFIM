#!/usr/bin/env python
# coding: utf-8
import os
import xesmf             as xe
import numpy             as np
import xarray            as xr
import cosima_cookbook   as cc
import matplotlib.pyplot as plt
import cartopy.crs       as ccrs
from dask.distributed import Client

V_era5       = '2t' 
V_bran       = 'temp'
yr_str       = '2005'
mo_str       = '01'
md_stt       = '01'
md_stp       = '31'
dt_stt       = '{yr:s}{mo:s}{md:s}'.format(yr=yr_str,mo=mo_str,md=md_stt)
dt_stp       = '{yr:s}{mo:s}{md:s}'.format(yr=yr_str,mo=mo_str,md=md_stp)
D_bran       = '/g/data/gb6/BRAN/BRAN2020/daily/'
D_era5       = '/g/data/rt52/era5/single-levels/reanalysis'
D_rgrd       = '/g/data/jk72/da1339/regrid/'
F_grd        = '/g/data/ik11/inputs/access-om2/input_20200530/cice_01deg/grid.nc'
F_Gbran      = '/g/data/jk72/da1339/grids/BRAN/ocean_grid.nc'
F_era5_t_wgt = '/g/data/jk72/da1339/grids/0p1/bilnr_wgts_xesmf_era5_0p25_1440x720_3600x7200_0p1_t.nc'
F_era5_u_wgt = '/g/data/jk72/da1339/grids/0p1/bilnr_wgts_xesmf_era5_0p25_1440x720_3600x7200_0p1_t.nc'
F_bran_t_wgt = '/g/data/jk72/da1339/grids/0p1/bilnr_wgts_xesmf_bran_0p15_3600x1500_3600x2700_0p1_t.nc'
F_bran_u_wgt = '/g/data/jk72/da1339/grids/0p1/bilnr_wgts_xesmf_bran_0p15_3600x1500_3600x2700_0p1_u.nc'
F_cice_t_grd = '/g/data/jk72/da1339/grids/CICE/0p1_tgrid.nc'
F_cice_u_grd = '/g/data/jk72/da1339/grids/CICE/0p1_ugrid.nc'
F_era5       = '{var:s}_era5_oper_sfc_{stt:s}-{stp:s}.nc'.format(var=V_era5,stt=dt_stt,stp=dt_stp)
F_atm        = 'JRA55_03hr_forcing_{yr:s}.nc'.format(yr=yr_str)
F_ocn        = 'BRAN_daily_forcing_{yr:s}.nc'.format(yr=yr_str)
F_bran       = 'ocean_{var:s}_{yr:s}_{mo:s}.nc'.format(var=V_bran,yr=yr_str,mo=mo_str) #ocean_eta_t_2005_01.nc
P_bran       = os.path.join(D_bran,F_bran)
D_atm        = os.path.join(D_rgrd,'ERA5')
D_ocn        = os.path.join(D_rgrd,'BRAN')
P_atm        = os.path.join(D_atm,F_atm)
P_ocn        = os.path.join(D_ocn,F_ocn)

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
    return cp * (sst-T_mld)

def compute_diffuse_sfc_em(em_sfc, em_toa):
    '''
    a simple function to compute the diffuse short/long wave EM at surface
    '''
    k_t = em_sfc / em_toa
    return 0.952 - 1.041 * np.exp(np.exp((2.3 - 4.702*k_t)))

def compute_u2_from_u10(u10):
    '''
    a simple function transfer/compute wind speed at 2-metres from reported wind speed at 10-metres
    '''
    return (u10 * 4.87) / np.log((67.8 * 10) - 5.42)

def compute_sfc_airrho(t2m, d2m, sp):
    '''
    compute atmospheric air density at the surface based on temperature, dewpoint and surface pressure using metpy calculation
    '''
    RH  = mpc.relative_humidity_from_dewpoint(t2m,d2m)
    r   = mpc.mixing_ratio_from_relative_humidity(sp, t2m, RH)
    return mpc.density(sp, t2m, r)

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

def ocean_temp_from_theta(p0, p):
    '''
    compute ocean temperature from potential temperature
    p0 :: reference pressure [db]
    p  :: pressure at level [db]

    '''
    dp = p0 - p

G_cice        = xr.open_dataset(F_grd,chunks={'ny':100,'nx':100}).rename({'tlon':'lon','tlat':'lat'})
G_cice        = xr.open_dataset(F_grd).rename({'tlon':'lon','tlat':'lat'})
G_cice = G_cice.drop(('clon_t','clat_t','clon_u','clat_u','ulat','ulon','htn','hte','angle','angleT','tarea','uarea'))
ln_tmp = G_cice.tlon.values[0,:]
lt_tmp = G_cice.tlat.values[:,0]
G_cice = G_cice.drop(('tlon','tlat'))
G_cice['lon']                = ln_tmp
G_cice['lat']                = lt_tmp
G_cice['lon']                = G_cice.lon * 180/np.pi
G_cice['lat']                = G_cice.lat * 180/np.pi
G_cice['lon'].attrs['units'] = 'degrees_east'
G_cice['lat'].attrs['units'] = 'degrees_north'
G_cice.to_netcdf(F_cice_t_grd)

# P_era5 = os.path.join(D_era5,V_era5,yr_str,F_era5)
# G_era5 = xr.open_dataset(P_era5).rename({'longitude': 'lon', 'latitude': 'lat'})
# rg     = xe.Regridder(G_era5,G_cice,method="bilinear",periodic=True,filename='/g/data/jk72/da1339/grids/0p1/bilnr_wgts_xesmf_era5_0p25_1440x720_3600x7200_0p1.nc')
# t2m_na      = xr.open_mfdataset(os.path.join(D_era5,'2t',yr_str,'*.nc'), parallel=True, chunks={"time": 1})
# t2m         = rg(t2m_na)
# d2m_na      = xr.open_mfdataset(os.path.join(D_era5,'2d',yr_str,'*.nc'), parallel=True, chunks={"time": 1})
# d2m         = rg(d2m_na)
# u10_na      = xr.open_mfdataset(os.path.join(D_era5,'10u',yr_str,'*.nc'), parallel=True, chunks={"time": 1})
# u10         = rg(u10_na)
# v10_na      = xr.open_mfdataset(os.path.join(D_era5,'10v',yr_str,'*.nc'), parallel=True, chunks={"time": 1})
# v10         = rg(v10_na)
# metss_na    = xr.open_mfdataset(os.path.join(D_era5,'metss',yr_str,'*.nc'), parallel=True, chunks={"time": 1})
# metss       = rg(metss_na)
# mntss_na    = xr.open_mfdataset(os.path.join(D_era5,'mntss',yr_str,'*.nc'), parallel=True, chunks={"time": 1})
# mntss       = rg(mntss_na)
# mror_na     = xr.open_mfdataset(os.path.join(D_era5,'mror',yr_str,'*.nc'), parallel=True, chunks={"time": 1})
# mror        = rg(mror_na)
# msdrswrf_na = xr.open_mfdataset(os.path.join(D_era5,'msdrswrf',yr_str,'*.nc'), parallel=True, chunks={"time": 1})
# msdrswrf    = rg(msdrswrf_na)
# msdwlwrf_na = xr.open_mfdataset(os.path.join(D_era5,'msdwlwrf',yr_str,'*.nc'), parallel=True, chunks={"time": 1})
# msdwlwrf    = rg(msdwlwrf_na)
# msdwswrf_na = xr.open_mfdataset(os.path.join(D_era5,'msdwswrf',yr_str,'*.nc'), parallel=True, chunks={"time": 1})
# msdwswrf    = rg(msdwswrf_na)
# msl_na      = xr.open_mfdataset(os.path.join(D_era5,'msl',yr_str,'*.nc'), parallel=True, chunks={"time": 1})
# msl         = rg(msl_na)
# msr_na      = xr.open_mfdataset(os.path.join(D_era5,'msr',yr_str,'*.nc'), parallel=True, chunks={"time": 1})
# msr         = rg(msr_na)
# msror_na    = xr.open_mfdataset(os.path.join(D_era5,'msror',yr_str,'*.nc'), parallel=True, chunks={"time": 1})
# msror       = rg(msror_na)
# mtdwswrf_na = xr.open_mfdataset(os.path.join(D_era5,'mtdwswrf',yr_str,'*.nc'), parallel=True, chunks={"time": 1})
# mtdwswrf    = rg(mtdwswrf_na)
# mtnlwrf_na  = xr.open_mfdataset(os.path.join(D_era5,'mtnlwrf',yr_str,'*.nc'), parallel=True, chunks={"time": 1})
# mtnlwrf     = rg(mtnlwrf_na)
# mtnswrf_na  = xr.open_mfdataset(os.path.join(D_era5,'mtnswrf',yr_str,'*.nc'), parallel=True, chunks={"time": 1})
# mtnswrf     = rg(mtnswrf_na)
# mtpr_na     = xr.open_mfdataset(os.path.join(D_era5,'mtpr',yr_str,'*.nc'), parallel=True, chunks={"time": 1})
# mtpr        = rg(mtpr_na)
# sp_na       = xr.open_mfdataset(os.path.join(D_era5,'sp',yr_str,'*.nc'), parallel=True, chunks={"time": 1})
# sp          = rg(sp_na)
# z_na        = xr.open_mfdataset(os.path.join(D_era5,'z',yr_str,'*.nc'), parallel=True, chunks={"time": 1})
# z           = rg(z_na)
# qsat  = compute_sfc_qsat(d2m.d2m, sp.sp)
# ln    = t2m_na.longitude
# lt    = t2m_na.latitude
# LN,LT = np.meshgrid(ln,lt)
# LN_rg = rg(LN)
# LT_rg = rg(LT)
# time  = sp.time
# ATM   = xr.Dataset({ "airtmp" : t2m.t2m,
#                      "dlwsfc" : msdwlwrf.msdwlwrf,
#                      "glbrad" : msdwswrf.msdwswrf,
#                      "spchmd" : qsat,
#                      "ttlpcp" : mtpr.mtpr,
#                      "wndewd" : u10.u10,
#                      "wndnwd" : v10.v10 },
#                    coords = { "LON"  : (["ny","nx"],LON_rg),
#                               "LAT"  : (["ny","nx"],LAT_rg),
#                               "time" : time.values })
# ATM    = ATM.rename_dims({"ny":"nj","nx":"ni"})
# ATM.to_netcdf(P_atm)

G_bran = xr.open_dataset(F_Gbran).rename({'xt_ocean': 'lon', 'yt_ocean': 'lat'}).chunk({'lon':100,'lat':100})
G_bran = G_bran.drop(('st_edges_ocean','tmask','umask','kmu','kmt','hu','ht','xu_ocean','yu_ocean'))
rg      = xe.Regridder(G_bran,G_cice,method="bilinear",periodic=True,filename=F_bran_t_wgt,reuse_weights=True)
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
sst     = temp.sel(st_ocean=0,method='nearest')
sss     = salt.sel(st_ocean=0,method='nearest')
qdp     = compute_ocn_heat_flux(temp,salt,mld,lat_name='lat')
qdp.to_netcdf(os.path.join(D_ocn,'qdp_{yr:s}.nc'.format(yr=yr_str)))
dhdx  = compfute_ocn_sfc_slope(eta,G_cice,direction='x',grid_scale_factor=100)
dhdy  = compfute_ocn_sfc_slope(eta,G_cice,direction='y',grid_scale_factor=100)
ln    = sst.lon
lt    = sst.lat
LN,LT = np.meshgrid(ln,lt)
LN_rg = rg(LN)
LT_rg = rg(LT)
time  = sst.time
OCN   = xr.Dataset({ "T"    : sst,
                     "S"    : sss,
                     "hblt" : mld,
                     "u"    : uocn,
                     "v"    : vocn,
                     "dhdx" : dhdx,
                     "dhdy" : dhdy,
                     "qdp"  : qdp},
                   coords = { "LON"  : (["ny","nx"],LON_rg),
                              "LAT"  : (["ny","nx"],LAT_rg),
                              "time" : time.values })
OCN    = OCN.rename_dims({"ny":"nj","nx":"ni"})
OCN.to_netcdf(P_ocn)

