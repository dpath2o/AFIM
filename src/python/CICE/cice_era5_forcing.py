
import os
import xarray as xr
import numpy  as np
import pandas as pd
from afim import regrid_netcdf_with_cdo

#
D02    = os.path.join('/','Volumes','ioa02')
D_ERA5 = os.path.join(D02,'reanalysis','ERA5')
Ftgrd  = os.path.join(D02,'model_input','grids','0p25','g0p25_cice_tgrid.nc')
Fugrd  = os.path.join(D02,'model_input','grids','0p25','g0p25_cice_ugrid.nc')
D_reG  = os.path.join(D02,'model_input','grids','0p25','regrid','ERA5')
D_frce = os.path.join('/','Users','dpath2o','cice-dirs','input','CICE_data','forcing','0p25','hourly')

#
G025u = xr.open_dataset(os.path.join(D02,'model_input','grids','0p25','g0p25_cice_ugrid.nc'))
G025t = xr.open_dataset(os.path.join(D02,'model_input','grids','0p25','g0p25_cice_tgrid.nc'))

#
Fin   = os.path.join(D_ERA5, '2d', '2010', '2d_era5_oper_sfc_20100101-20100131.nc')
Fout  = os.path.join(D_reG, 'g0p25_d2m.nc')
regrid_netcdf_with_cdo(Fin, Fout, Ftgrd, cdo_options='-f nc -b F64')

#
Fin   = os.path.join(D_ERA5,'mror','2010','mror_era5_oper_sfc_20100101-20100131.nc')
Fout  = os.path.join(D_reG,'g0p25_mror.nc')
regrid_netcdf_with_cdo(Fin, Fout, Ftgrd, cdo_options='-f nc -b F64')

#
Fin   = os.path.join(D_ERA5,'2t','2010','2t_era5_oper_sfc_20100101-20100131.nc')
Fout  = os.path.join(D_reG,'g0p25_t2m.nc')
regrid_netcdf_with_cdo(Fin, Fout, Ftgrd, cdo_options='-f nc -b F64')

#
Fin   = os.path.join(D_ERA5,'msdrswrf','2010','msdrswrf_era5_oper_sfc_20100101-20100131.nc')
Fout  = os.path.join(D_reG,'g0p25_msdrswrf.nc')
regrid_netcdf_with_cdo(Fin, Fout, Ftgrd, cdo_options='-f nc -b F64')

#
Fin   = os.path.join(D_ERA5,'msdwlwrf','2010','msdwlwrf_era5_oper_sfc_20100101-20100131.nc')
Fout  = os.path.join(D_reG,'g0p25_msdwlwrf.nc')
regrid_netcdf_with_cdo(Fin, Fout, Ftgrd, cdo_options='-f nc -b F64')

#
Fin   = os.path.join(D_ERA5,'msdwswrf','2010','msdwswrf_era5_oper_sfc_20100101-20100131.nc')
Fout  = os.path.join(D_reG,'g0p25_msdwswrf.nc')
regrid_netcdf_with_cdo(Fin, Fout, Ftgrd, cdo_options='-f nc -b F64')

#
Fin   = os.path.join(D_ERA5,'msl','2010','msl_era5_oper_sfc_20100101-20100131.nc')
Fout  = os.path.join(D_reG,'g0p25_msl.nc')
regrid_netcdf_with_cdo(Fin, Fout, Ftgrd, cdo_options='-f nc -b F64')

#
Fin   = os.path.join(D_ERA5,'msnlwrf','2010','msnlwrf_era5_oper_sfc_20100101-20100131.nc')
Fout  = os.path.join(D_reG,'g0p25_msnlwrf.nc')
regrid_netcdf_with_cdo(Fin, Fout, Ftgrd, cdo_options='-f nc -b F64')

#
Fin   = os.path.join(D_ERA5,'msnswrf','2010','msnswrf_era5_oper_sfc_20100101-20100131.nc')
Fout  = os.path.join(D_reG,'g0p25_msnswrf.nc')
regrid_netcdf_with_cdo(Fin, Fout, Ftgrd, cdo_options='-f nc -b F64')

#
Fin   = os.path.join(D_ERA5,'msr','2010','msr_era5_oper_sfc_20100101-20100131.nc')
Fout  = os.path.join(D_reG,'g0p25_msr.nc')
regrid_netcdf_with_cdo(Fin, Fout, Ftgrd, cdo_options='-f nc -b F64')

#
Fin   = os.path.join(D_ERA5,'msror','e5.oper.fc.sfc.meanflux.235_020_msror.ll025sc.2009121606_2010010106.nc')
Fout  = os.path.join(D_reG,'g0p25_msror.nc')
regrid_netcdf_with_cdo(Fin, Fout, Ftgrd, cdo_options='-f nc -b F64')

#
Fin   = os.path.join(D_ERA5,'mtdwswrf','2010','mtdwswrf_era5_oper_sfc_20100101-20100131.nc')
Fout  = os.path.join(D_reG,'g0p25_mtdwswrf.nc')
regrid_netcdf_with_cdo(Fin, Fout, Ftgrd, cdo_options='-f nc -b F64')

#
Fin   = os.path.join(D_ERA5,'mtnlwrf','2010','mtnlwrf_era5_oper_sfc_20100101-20100131.nc')
Fout  = os.path.join(D_reG,'g0p25_mtnlwrf.nc')
regrid_netcdf_with_cdo(Fin, Fout, Ftgrd, cdo_options='-f nc -b F64')

#
Fin   = os.path.join(D_ERA5,'mtnswrf','2010','mtnswrf_era5_oper_sfc_20100101-20100131.nc')
Fout  = os.path.join(D_reG,'g0p25_mtnswrf.nc')
regrid_netcdf_with_cdo(Fin, Fout, Ftgrd, cdo_options='-f nc -b F64')

#
Fin   = os.path.join(D_ERA5,'mtpr','2010','mtpr_era5_oper_sfc_20100101-20100131.nc')
Fout  = os.path.join(D_reG,'g0p25_mtpr.nc')
regrid_netcdf_with_cdo(Fin, Fout, Ftgrd, cdo_options='-f nc -b F64')

#
Fin   = os.path.join(D_ERA5,'sp','2010','sp_era5_oper_sfc_20100101-20100131.nc')
Fout  = os.path.join(D_reG,'g0p25_sp.nc')
regrid_netcdf_with_cdo(Fin, Fout, Ftgrd, cdo_options='-f nc -b F64')

# ERA5 u-grid re-gridding
#
Fin   = os.path.join(D_ERA5,'metss','2010','metss_era5_oper_sfc_20100101-20100131.nc')
Fout  = os.path.join(D_reG,'g0p25_metss.nc')
regrid_netcdf_with_cdo(Fin, Fout, Fugrd, cdo_options='-f nc -b F64')

#
Fin   = os.path.join(D_ERA5,'10v','e5.oper.an.sfc.128_166_10v.ll025sc.2010010100_2010013123.nc')
Fout  = os.path.join(D_reG,'g0p25_v10.nc')
regrid_netcdf_with_cdo(Fin, Fout, Fugrd, cdo_options='-f nc -b F64')

#
Fin   = os.path.join(D_ERA5,'mntss','2010','mntss_era5_oper_sfc_20100101-20100131.nc')
Fout  = os.path.join(D_reG,'g0p25_mntss.nc')
regrid_netcdf_with_cdo(Fin, Fout, Fugrd, cdo_options='-f nc -b F64')

#
Fin   = os.path.join(D_ERA5,'10u','e5.oper.an.sfc.128_165_10u.ll025sc.2010010100_2010013123.nc')
Fout  = os.path.join(D_reG,'g0p25_u10.nc')
regrid_netcdf_with_cdo(Fin, Fout, Fugrd, cdo_options='-f nc -b F64')

##################
# 
msdrswrf = xr.open_dataset(os.path.join(D_reG,'g0p25_msdrswrf.nc'))
mtnswrf  = xr.open_dataset(os.path.join(D_reG,'g0p25_mtnswrf.nc'))
k_t      = msdrswrf.msdrswrf / mtnswrf.mtnswrf
SW_diff  = 0.952 - 1.041 * np.exp(np.exp((2.3 - 4.702*k_t)))
SW_diff.to_netcdf(os.path.join(D_reG,'g0p25_swvisdf.nc'))


#
u10 = xr.open_dataset(os.path.join(D_reG,'g0p25_u10.nc'))
u2  = (u10.VAR_10U * 4.87) / np.log((67.8 * 10) - 5.42)
u2.to_netcdf(os.path.join(D_reG,'g0p25_u2.nc'))

#
v10 = xr.open_dataset(os.path.join(D_reG,'g0p25_v10.nc'))
v2  = (v10.VAR_10V * 4.87) / np.log((67.8 * 10) - 5.42)
v2.to_netcdf(os.path.join(D_reG,'g0p25_v2.nc'))

#
#density of air at 2-metres
t2m = xr.open_dataset(os.path.join(D_reG,'g0p25_t2m.nc'))
d2m = xr.open_dataset(os.path.join(D_reG,'g0p25_d2m.nc'))
sp  = xr.open_dataset(os.path.join(D_reG,'g0p25_sp.nc'))
RH  = mpc.relative_humidity_from_dewpoint(t2m.t2m,d2m.d2m)
r   = mpc.mixing_ratio_from_relative_humidity(sp.sp, t2m.t2m, RH)
rho = mpc.density(sp.sp, t2m.t2m, r)
rho.to_netcdf(os.path.join(D_reG,'g0p25_rho.nc'))

# specific humidity at 2-metres
Rdry = 287.0597
Rvap = 461.5250
a1   = 611.21
a3   = 17.502
a4   = 32.19
T0   = 273.16
E    = a1 * np.exp(a3 * (d2m.d2m-T0) / (d2m.d2m-a4) )
qsat = (Rdry/Rvap) * E / (sp.sp - ( (1-Rdry/Rvap) * E) )
qsat.to_netcdf(os.path.join(D_reG,'g0p25_qsat.nc'))

#
glbrad = xr.open_dataset(os.path.join(D_reG,'g0p25_msdrswrf.nc')).msdrswrf.values
dlwsfc = xr.open_dataset(os.path.join(D_reG,'g0p25_msdwlwrf.nc')).msdwlwrf.values
wndewd = xr.open_dataset(os.path.join(D_reG,'g0p25_u2.nc')).VAR_10U.values
wndnwd = xr.open_dataset(os.path.join(D_reG,'g0p25_v2.nc')).VAR_10V.values
airtmp = xr.open_dataset(os.path.join(D_reG,'g0p25_t2m.nc')).t2m.values
spchmd = xr.open_dataset(os.path.join(D_reG,'g0p25_qsat.nc')).values
ttlpcp = xr.open_dataset(os.path.join(D_reG,'g0p25_mtpr.nc')).mtpr.values

# 
ERA5 = ''
ERA5 = xr.Dataset({'glbrad' : (['time','nj','ni'],glbrad),
                   'dlwsfc' : (['time','nj','ni'],dlwsfc),
                   'wndewd' : (['time','nj','ni'],wndewd),
                   'wndnwd' : (['time','nj','ni'],wndnwd),
                   'airtmp' : (['time','nj','ni'],airtmp),
                   'spchmd' : (['time','nj','ni'],spchmd),
                   'ttlpcp' : (['time','nj','ni'],ttlpcp),},
                   coords = {'ulon' : (['nj'],G025u.lon.data[0,:]),
                             'ulat' : (['ni'],G025u.lat.data[:,0]),
                             'tlon' : (['nj'],G025t.lon.data[0,:]),
                             'tlat' : (['ni'],G025t.lat.data[:,0]),
                             'time' : pd.date_range('2010-01-01', freq='H', periods=31*24)})
ERA5.to_netcdf(os.path.join(D_frce,'ERA5_g0p25_2010.nc'))
