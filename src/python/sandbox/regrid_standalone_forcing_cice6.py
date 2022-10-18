#!/usr/bin/env python
# coding: utf-8
import os
import glob
import xesmf             as xe
import numpy             as np
import xarray            as xr
#import cosima_cookbook   as cc
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


G_cice        = xr.open_dataset(F_grd,chunks={'ny':100,'nx':100})
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

P_era5 = os.path.join(D_era5,V_era5,yr_str,F_era5)
G_era5 = xr.open_dataset(P_era5).rename({'longitude': 'lon', 'latitude': 'lat'})
rg     = xe.Regridder(G_era5,G_cice,method="bilinear",periodic=True,filename=F_era5_t_wgt,reuse_weights=True)
t2m_na = xr.open_mfdataset(os.path.join(D_era5,'2t',yr_str,'*.nc'), parallel=True, chunks={"time": 1})
t2m     = rg(t2m_na)
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
for i in range(90):
    F_png = f"./GRAPHICAL/animations/native_sequence/Python_Animation_03_frame_{i:04}.png"
    if not os.path.exists(F_png):
        print('plotting: '+str(t2m_na.coords['time'].values[i])[:13])
        t2m.isel(time=i).t2m.plot( figsize = (12,6),
                                   vmin=190, vmax=340,
                                   cmap='coolwarm',    # Change the colormap back to 'bwr'
                                   cbar_kwargs={ 'extend':'neither' } )
        plt.title("Time = " + str(t2m.coords['time'].values[i])[:13])
        plt.savefig(F_png)
        plt.close()
os.system("ffmpeg -r 60 -f image2 -s 1920x1080 -i ./GRAPHICAL/animations/native_sequence/Python_Animation_03_frame_%04d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p ./era5_native_t2m_2005_1st90days.mp4")
F_del = glob.glob('./GRAPHICAL/animations/native_sequence/*.png')
for F_delete in F_del: os.remove(F_delete)
for i in range(90*24):
    F_png = f"./GRAPHICAL/animations/regrid_sequence/Python_Animation_03_frame_{i:04}.png"
    if not os.path.exists(F_png):
        print('plotting: '+str(t2m.coords['time'].values[i])[:13])
        t2m.isel(time=i).t2m.plot( figsize = (12,6),
                                   vmin=190, vmax=340,
                                   cmap='coolwarm',    # Change the colormap back to 'bwr'
                                   cbar_kwargs={ 'extend':'neither' } )
        plt.title("Time = " + str(t2m.coords['time'].values[i])[:13])
        plt.savefig(F_png)
        plt.close()
os.system("ffmpeg -r 60 -f image2 -s 1920x1080 -i ../GRAPHICAL/animations/regrid_sequence/Python_Animation_03_frame_%04d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p ./era5_regrid_t2m_2005_1st90days.mp4")
F_del = glob.glob('./GRAPHICAL/animations/regrid_sequence/*.png')
for F_delete in F_del: os.remove(F_delete)
sss     = salt.sel(st_ocean=0,method='nearest')
print("I made it to here ..")
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

