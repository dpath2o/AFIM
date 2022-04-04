#!/usr/bin/env python

import os
import numpy as np
import xarray as xr
from affine import Affine
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cmocean
import cmocean.cm as cmo
from osgeo import gdal, osr

gdal.UseExceptions()

Ddat  = os.path.join(os.path.expanduser('~'),'SEASONG','PHD','data','bathy_topo','Cape_Darnley_Bathymetry_Grid_v1','03_GEOTIFF')
Fdat  = 'cdbg_v1.tif'

# method one
da = xr.open_rasterio(os.path.join(Ddat,Fdat))
xfrm  = Affine.from_gdal(*da.variable.attrs['transform'])
nx,ny = da.sizes['x'], da.sizes['y']
x,y   = np.meshgrid(np.arange(nx), np.arange(ny)) * xfrm

# method two
ds = gdal.Open(os.path.join(Ddat,Fdat))
data = ds.ReadAsArray()
gt = ds.GetGeoTransform()
proj = ds.GetProjection()
inproj = osr.SpatialReference()
inproj.ImportFromWkt(proj)
projcs = inproj.GetAuthorityCode('PROJCS')
projection = ccrs.epsg(int(projcs))
subplot_kw = dict(projection=projection)

clevs = np.arange(-3000, 0, 200)
fig, ax = plt.subplots(figsize=(9, 9), subplot_kw=subplot_kw)
extent = (gt[0], gt[0] + ds.RasterXSize * gt[1],
          gt[3] + ds.RasterYSize * gt[5], gt[3])
img = ax.imshow(da.variable.data[0], extent=extent, origin='upper', cmap=cmo.ice , vmin=-3000, vmax=0,)
plt.colorbar(img)
ax.coastlines(resolution='10m')#, alpha=0.1)
plt.show()
exit()

# import plotly.graph_objects as go
# import pandas as pd
# import numpy as np
# # Read data from a csv
# z_data = pd.read_csv('bathy_bedford.csv')
# z = z_data.values
# sh_0, sh_1 = z.shape
# x, y = np.linspace(44.66875, 44.74791667, sh_0), np.linspace(-63.69791667, -63.52708333, sh_1)
# fig = go.Figure(data=[go.Surface(z=z, x=x, y=y)])
# fig.update_traces(contours_z=dict(show=True, usecolormap=True,
#                                   highlightcolor="limegreen", project_z=True))
# fig.update_layout(title='Bedford Basin Elevation',xaxis_title="Latitude",
#                   yaxis_title="Longitude",autosize=False,
#                   width=900, height=900,
#                   margin=dict(l=65, r=50, b=65, t=90))
# fig.update_layout(scene = dict(
#                     xaxis_title='Latitude',
#                     yaxis_title='Longitude',
#                     zaxis_title='Elevation'),
#                     margin=dict(r=20, b=10, l=10, t=10))
# # fig.update_layout(color='Elevation')
# fig.update_layout(coloraxis_colorbar=dict(
#     title="Elevation",
#     thicknessmode="pixels", thickness=50,
#     lenmode="pixels", len=200,
#     yanchor="top", y=1,
#     ticks="outside", ticksuffix="",
#     dtick=5
# ))
# fig.show()

# PLOT
#crs = ccrs.PlateCarree()
crs = ccrs.UTM(zone="42D")
fig = plt.figure()#figsize=(10,10))
ax  = fig.add_subplot(111, projection=crs)
ax.coastlines(resolution='10m', alpha=0.1)
#ax.contourf(x, y, da.variable.data[0], cmap='Greys')
#ax.set_extent([gt[0], gt[0] + ds.RasterXSize * gt[1], gt[3] + ds.RasterYSize * gt[5], gt[3]])
#ax.set_extent([lon_min, lon_max, lat_min, lat_max])
#gl = ax.gridlines(crs=crs, draw_labels=True, alpha=0.5)
#gl.xlabels_top = None
#gl.ylabels_right = None
#xgrid = np.arange(lon_min-0.5, lon_max+0.5, 1.)
#ygrid = np.arange(lat_min, lat_max+1, 1.)
#gl.xlocator = mticker.FixedLocator(xgrid.tolist())
#gl.ylocator = mticker.FixedLocator(ygrid.tolist())
#gl.xformatter = LONGITUDE_FORMATTER
#gl.yformatter = LATITUDE_FORMATTER
#gl.xlabel_style = {'size': 14, 'color': 'black'}
#gl.ylabel_style = {'size': 14, 'color': 'black'}
plt.show()

# fig = plt.figure(figsize=(12,12))
# ax = fig.add_subplot(111)
# ax.imshow(da.variable.data[0])
# plt.show()
