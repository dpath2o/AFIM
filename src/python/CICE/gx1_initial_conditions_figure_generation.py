import os
import xarray            as xr
import matplotlib.pyplot as plt
import cartopy.crs       as ccrs

#

##################################################################################
# initial condition NON-tripolar grid provided from CICE6
F_gx1ic = '/Users/dpath2o/cice-dirs/input/CICE_data/ic/gx1/iced_gx1_v6.2005-01-01.nc'

fig = plt.figure(figsize=(12,10))
ax = fig.add_subplot(1,1,1, projection=ccrs.Robinson())
