#!/usr/bin/env python3

import os
import numpy as np
from pyproj import Proj
from pyproj.crs import ProjectedCRS
from pyproj.crs.coordinate_operation import UTMConversion
import rioxarray
import xarray
from affine import Affine
import matplotlib.pyplot as plt

Ddat   = os.path.join(os.path.expanduser('~'),'PhD','data','bathy_topo','Cape_Darnley_Bathymetry_Grid_v1','01_ESRI_ASCII')
Fdat   = open(os.path.join(Ddat,'cape_darnley_bathymetry_grid_v1.asc'), 'r')
Fprj   = open(os.path.join(Ddat,'cape_darnley_bathymetry_grid_v1.prj'), 'r')
strHdr = Fdat.readlines()[:5]
strPrj = Fprj.readlines()
Fdat.close()
Fprj.close()
Ncols  = int(strHdr[0].split()[1])
Nrows  = int(strHdr[1].split()[1])
xlcrn  = int(strHdr[2].split()[1])
ylcrn  = int(strHdr[3].split()[1])
Nclsz  = int(strHdr[4].split()[1])
prjct  = strPrj[0].split()[1]
zone   = strPrj[1].split()[1]
datum  = strPrj[2].split()[1]
sphrd  = strPrj[3].split()[1]
units  = strPrj[4].split()[1]
zunit  = strPrj[5].split()[1]
yshift = strPrj[6].split()[1]
Fdat   = open(os.path.join(Ddat,'cape_darnley_bathymetry_grid_v1.asc'), 'r')
strDat = Fdat.readlines()[7:]
Fdat.close()
Zdat = [np.fromstring(ln, dtype=float, sep=' ') for ln in strDat]
Z    = np.array(Zdat)
#Z    = np.transpose(Z)
x    = np.arange(xlcrn,xlcrn+(Z.shape[0]*Nclsz),Nclsz)
y    = np.arange(ylcrn,ylcrn+(Z.shape[1]*Nclsz),Nclsz)
X,Y  = np.meshgrid(x,y)

# plt.figure()
# lims = dict(cmap='RdBu_r', vmin=-3000, vmax=0)
# plt.pcolormesh(X,Y,Z,shading='auto',**lims)
# plt.colorbar()
# plt.show()

# Step 1: Build the transform
# GDAL link providing some insight on geo-transformation of raster data
# https://gdal.org/user/raster_data_model.html#affine-geotransform
#top = bottom - resolution_y * height
top = ylcrn + (Nclsz*Nrows)
transform  = Affine.translation(xlcrn,top) * Affine.scale(Nclsz,-Nclsz)

# Step 2: Get the CRS from UTM Zone
crs = ProjectedCRS(conversion=UTMConversion(-42))

# Step 3: Assign the transform/affine & CRS to raster
rds = xarray.DataArray(Z,  dims=("y", "x"))
rds.rio.write_transform(transform, inplace=True)
rds.rio.write_crs(crs, inplace=True)

# Step 4: Re-project to geographic
rds = rds.rio.reproject("EPSG:4326")
