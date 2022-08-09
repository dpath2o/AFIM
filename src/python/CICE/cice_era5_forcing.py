
import os
import pdb
import xarray     as xr
import numpy      as np
import pandas     as pd
import metpy.calc as mpc
from afim import regrid_netcdf_with_cdo

# Directories, files and switches
compute_sw_diffuse = False
compute_u2v2       = False
compute_rho        = False
compute_qsat       = False
Gtype  = '0p25'
regrd  = False
sttdt  = '2010-01-01'
months = 12
years  = 1
D02    = os.path.join('/','Volumes','ioa02')
D_ERA5 = os.path.join(D02,'reanalysis','ERA5')
Ftgrd  = os.path.join(D02,'model_input','grids',Gtype,'g{:s}_cice_tgrid.nc'.format(Gtype))
Fugrd  = os.path.join(D02,'model_input','grids',Gtype,'g{:s}_cice_ugrid.nc'.format(Gtype))
D_reG  = os.path.join(D02,'model_input','grids',Gtype,'regrid','ERA5')
D_frce = os.path.join('/','Users','dpath2o','cice-dirs','input','CICE_data','forcing',Gtype,'hourly')

# Load grid files
Gu = xr.open_dataset(os.path.join(D02,'model_input','grids',Gtype,'g{:s}_cice_ugrid.nc'.format(Gtype)))
Gt = xr.open_dataset(os.path.join(D02,'model_input','grids',Gtype,'g{:s}_cice_tgrid.nc'.format(Gtype)))

# regridding
tvars = ['2d','2t','msdrswrf','msdwlwrf','sp','mtpr','mtnswrf']
uvars = ['10u','10v']
varss = tvars+uvars
stdts = pd.date_range(sttdt, freq='MS', periods=months*years)
spdts = pd.date_range(sttdt, freq='M', periods=months*years)
Fstrs = 'era5_oper_sfc'
if regrd:
    for i in tvars:
        for j in np.arange(0,len(stdts),1):
            if spdts[j].is_year_end:
                DS = xr.open_mfdataset( os.path.join(D_reG,stdts[i].strftime('%Y'),i,'*.nc') )
                DS.to_netcdf( os.path.join(D_reG,stdts[i].strftime('%Y'),'g{Gtype:s}_{field:s}.nc'.format(Gtype=Gtype,field=i)) )
            Fin = '{field:s}_{Fstrs:s}_{stt:s}-{stp:s}.nc'.format(field=j,
                                                                  Fstrs=Fstrs,
                                                                  stt=stdts[j].strftime('%Y%M%d'),
                                                                  stp=spdts[j].strftime('%Y%M%d'))
            Pin  = os.path.join(D_ERA5,j,stdts[i].strftime('%Y'),Fin)
            Dout = os.path.join(D_reG,stdts[i].strftime('%Y'),i)
            if not(os.path.exists(Dout)): os.mkdir(Dout)
            Fout = 'g{Gtype:s}_{field:s}_{yrmo:s}'.format(Gtype=Gtype,field=i,yrmo=stdts[j].strftime('%Y%M'))
            Pout = os.path.join(Dout,Fout)
            regrid_netcdf_with_cdo(Pin, Pout, Ftgrd, cdo_options='-f nc -b F64')
    for i in uvars:
        for j in np.arange(0,len(stdts),1):
            if spdts[j].is_year_end:
                DS = xr.open_mfdataset( os.path.join(D_reG,stdts[i].strftime('%Y'),i,'*.nc') )
                DS.to_netcdf( os.path.join(D_reG,stdts[i].strftime('%Y'),'g{Gtype:s}_{field:s}.nc'.format(Gtype=Gtype,field=i)) )
            Fin = '{field:s}_{Fstrs:s}_{stt:s}-{stp:s}.nc'.format(field=j,
                                                                  Fstrs=Fstrs,
                                                                  stt=stdts[j].strftime('%Y%M%d'),
                                                                  stp=spdts[j].strftime('%Y%M%d'))
            Pin  = os.path.join(D_ERA5,j,stdts[i].strftime('%Y'),Fin)
            Dout = os.path.join(D_reG,stdts[i].strftime('%Y'),i)
            if not(os.path.exists(Dout)): os.mkdir(Dout)
            Fout = 'g{Gtype:s}_{field:s}_{yrmo:s}'.format(Gtype=Gtype,field=i,yrmo=stdts[j].strftime('%Y%M'))
            Pout = os.path.join(Dout,Fout)
            regrid_netcdf_with_cdo(Pin, Pout, Fugrd, cdo_options='-f nc -b F64')


####################################################
# shortwave diffuse calculation
if compute_sw_diffuse:
    msdrswrf = xr.open_dataset(os.path.join(D_reG,'g0p25_msdrswrf.nc'))
    mtnswrf  = xr.open_dataset(os.path.join(D_reG,'g0p25_mtnswrf.nc'))
    k_t      = msdrswrf.msdrswrf / mtnswrf.mtnswrf
    SW_diff  = 0.952 - 1.041 * np.exp(np.exp((2.3 - 4.702*k_t)))
    SW_diff.to_netcdf(os.path.join(D_reG,'g0p25_swvisdf.nc'))
    SW_diff  = ''
    msdrswrf = ''
    mtnswrf  = ''
    k_t      = ''

# U10 to U2
if compute_u2v2:
    u10 = xr.open_dataset(os.path.join(D_reG,'g0p25_u10.nc'))
    u2  = (u10.VAR_10U * 4.87) / np.log((67.8 * 10) - 5.42)
    u2.to_netcdf(os.path.join(D_reG,'g0p25_u2.nc'))
    u10 = ''
    u2  = ''
    v10 = xr.open_dataset(os.path.join(D_reG,'g0p25_v10.nc'))
    v2  = (v10.VAR_10V * 4.87) / np.log((67.8 * 10) - 5.42)
    v2.to_netcdf(os.path.join(D_reG,'g0p25_v2.nc'))
    v10 = ''
    v2  = ''

# density of air at 2-metres
if compute_rho:
    t2m = xr.open_dataset(os.path.join(D_reG,'g0p25_t2m.nc'))
    d2m = xr.open_dataset(os.path.join(D_reG,'g0p25_d2m.nc'))
    sp  = xr.open_dataset(os.path.join(D_reG,'g0p25_sp.nc'))
    RH  = mpc.relative_humidity_from_dewpoint(t2m.t2m,d2m.d2m)
    r   = mpc.mixing_ratio_from_relative_humidity(sp.sp, t2m.t2m, RH)
    rho = mpc.density(sp.sp, t2m.t2m, r)
    rho.to_netcdf(os.path.join(D_reG,'g0p25_rho.nc'))
    rho = ''
    RH  = ''
    r   = ''
    d2m = ''
    t2m = ''
    sp  = ''

# specific humidity at 2-metres
if compute_qsat:
    d2m = xr.open_dataset(os.path.join(D_reG,'g0p25_d2m.nc'))
    sp  = xr.open_dataset(os.path.join(D_reG,'g0p25_sp.nc'))
    Rdry = 287.0597
    Rvap = 461.5250
    a1   = 611.21
    a3   = 17.502
    a4   = 32.19
    T0   = 273.16
    E    = a1 * np.exp(a3 * (d2m.d2m-T0) / (d2m.d2m-a4) )
    qsat = (Rdry/Rvap) * E / (sp.sp - ( (1-Rdry/Rvap) * E) )
    qsat.to_netcdf(os.path.join(D_reG,'g0p25_qsat.nc'))
    qsat = ''
    E    = ''
    t2m  = ''
    d2m  = ''
    sp   = ''

#################################################################
# 
glbrad = xr.open_dataset(os.path.join(D_reG,'g0p25_msdrswrf.nc')).msdrswrf.values
dlwsfc = xr.open_dataset(os.path.join(D_reG,'g0p25_msdwlwrf.nc')).msdwlwrf.values
wndewd = xr.open_dataset(os.path.join(D_reG,'g0p25_u2.nc')).VAR_10U.values
wndnwd = xr.open_dataset(os.path.join(D_reG,'g0p25_v2.nc')).VAR_10V.values
airtmp = xr.open_dataset(os.path.join(D_reG,'g0p25_t2m.nc')).t2m.values
spchmd = xr.open_dataset(os.path.join(D_reG,'g0p25_qsat.nc')).__xarray_dataarray_variable__.values
ttlpcp = xr.open_dataset(os.path.join(D_reG,'g0p25_mtpr.nc')).mtpr.values
airrho = xr.open_dataset(os.path.join(D_reG,'g0p25_rho.nc')).__xarray_dataarray_variable__.values

# 
ERA5 = ''
ERA5 = xr.Dataset({'glbrad' : (['time','nj','ni'],glbrad),
                   'dlwsfc' : (['time','nj','ni'],dlwsfc),
                   'wndewd' : (['time','nj','ni'],wndewd),
                   'wndnwd' : (['time','nj','ni'],wndnwd),
                   'airtmp' : (['time','nj','ni'],airtmp),
                   'spchmd' : (['time','nj','ni'],spchmd),
                   'ttlpcp' : (['time','nj','ni'],ttlpcp),
                   'airrho' : (['time','nj','ni'],airrho)},
                   coords = {'ulon' : (['ni'],Gu.lon.data[0,:]),
                             'ulat' : (['nj'],Gu.lat.data[:,0]),
                             'tlon' : (['ni'],Gt.lon.data[0,:]),
                             'tlat' : (['nj'],Gt.lat.data[:,0]),
                             'time' : pd.date_range('2010-01-01', freq='H', periods=31*24)})
ERA5.to_netcdf(os.path.join(D_frce,'ERA5_g0p25_2010.nc'))
