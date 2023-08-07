import os
import re
import fnmatch
import glob
import sys
#!pip install proplot
sys.path.insert(0, '/home/581/da1339/.local/lib/python3.9/site-packages')
sys.path.insert(0, './src/python')
import afim
import cmocean
import xesmf             as xe
import numpy             as np
import pandas            as pd
import xarray            as xr
import cosima_cookbook   as cc
import matplotlib.pyplot as plt
import proplot           as plot
import cartopy.crs       as ccrs
import cartopy.crs       as ccrs
import cartopy.feature   as cfeature
import matplotlib.path   as mpath
from matplotlib          import pyplot as plt, animation
from dask.distributed    import Client
D_ani      = '/g/data/jk72/da1339/GRAPHICAL/animations'
D_afim     = '/g/data/jk72/da1339/afim_output'
D_ices     = ['20230703_jra55_aom2_ndte_240_mld_ice_dynamic',
              '20230706_jra55_aom2_ndte_3710_mld_ice_constant',
              '20230708_jra55_aom2_ndte_3710_mld_ice_dynamic',
              '20230708_jra55_bran_ndte_3710_mld_ice_dynamic']
var_names  = ['aice','frazil','hi','uvel','congel']
cmap       = cmocean.cm.ice
G_res      = '025'
figsize    = (12,9)
label_cbar = ['sea ice concentration (%)','frazil (cm/day)','ice volume per unit grid cell area','ice speed (m/s)','congelation ice growth']
vmin       = 0
vmax       = 1
frames     = 365
interval   = '-delay 1x30'
land_10m   = cfeature.NaturalEarthFeature('physical', 'land', '10m', edgecolor='black', facecolor='papayawhip', linewidth=0.5)
geolon     = xr.open_dataset('/g/data/ik11/grids/ocean_grid_{:s}.nc'.format(G_res)).geolon_t
geolat     = xr.open_dataset('/g/data/ik11/grids/ocean_grid_{:s}.nc'.format(G_res)).geolat_t
#########################################
os.system(f"rm {D_ani}/tmp/*.png")
for j in D_ices:
    ds = xr.open_mfdataset(f'{D_ani}/{j}/history/*.nc',chunks=None)
    ds = ds.assign_coords({'geolat': geolat, 'geolon_t': geolon})
    for k in var_names:
        var_name = k
        if var_name=='uvel':
            var_name = 'spd'
            ds = xr.Dataset(data_vars = {'spd'    : (('time','nj','ni'), np.sqrt(ds.uvel.data**2 + ds.vvel.data**2))},
                            coords    = {'geolon' : (['nj','ni'],ds.geolon.data,{'units':'degrees_east'}),
                                         'geolat' : (['nj','ni'],ds.geolat.data,{'units':'degrees_north'}),
                                         'time'   : (['time'],ds.time.data)}
        F_ani = '{:s}_{:s}_{:04d}.gif'.format(j,var_name,frames)
        for i in range(frames):
            f,axs = plot.subplots(ncols   = 2,
                                  axwidth = 2.2, 
                                  figsize = figsize,
                                  proj    = {1:'splaea',
                                             2:'npaeqd'}, 
                                  basemap = {1:False,
                                             2:False})
            axs.format(land     = True,
                       latlines = 10,
                       latmax   = 80,
                       suptitle = "{:s}".format(str(ds.coords['time'].values[i])[:13]))
            axs[0].format(boundinglat=-50, title='')
            axs[1].format(boundinglat=40 , title='')
            axs[0].add_feature(land_10m,color=[0.8, 0.8, 0.8])
            axs[0].coastlines(resolution='10m')
            axs[1].add_feature(land_10m,color=[0.8, 0.8, 0.8])
            axs[1].coastlines(resolution='10m')
            ds.isel(time=i)[var_name].plot(ax           = axs[0],
                                           x            = 'geolon_t',
                                           y            = 'geolat_t', 
                                           transform    = ccrs.PlateCarree(),
                                           vmin         = vmin,
                                           vmax         = vmax,
                                           cmap         = cmap,
                                           add_colorbar = False)
            ds.isel(time=i)[var_name].plot(ax          = axs[1],
                                           x           = 'geolon_t',
                                           y           = 'geolat_t',
                                           transform   = ccrs.PlateCarree(),
                                           vmin        = vmin,
                                           vmax        = vmax,
                                           cmap        = cmap,
                                           cbar_kwargs = {'label'    : label_cbar,
                                                          'fraction' : 0.03,
                                                          'aspect'   : 15,
                                                          'shrink'   : 0.7})
            plt.savefig(f"{D_ani}/tmp/{i:04}.png")
            plt.close(f)
        os.system(f"convert {interval} {D_ani}/tmp/*.png {D_ani}/{F_ani}")