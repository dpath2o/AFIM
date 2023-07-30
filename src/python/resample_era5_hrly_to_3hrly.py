import os
import xarray as xr
ds = xr.open_dataset('/scratch/jk72/da1339/cice-dirs/input/AFIM/forcing/0p25/ERA5/8XDAILY/ERA5_2005_with_jra55do_var_names_reG_nearest_s2d_0p25_aom2_01hr.nc')
ds_resampled = ds.resample(time='3H').mean()
ds_resampled.to_netcdf('/scratch/jk72/da1339/cice-dirs/input/AFIM/forcing/0p25/ERA5/8XDAILY/JRA55_03hr_forcing_tx1_2005.nc')