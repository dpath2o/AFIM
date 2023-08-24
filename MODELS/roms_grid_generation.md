# test correct anaconda environment

``` python
import sys
return sys.exec_prefix
```

# import required libraries

``` python
import warnings
warnings.simplefilter('ignore')
import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Polygon, Point
from pyproj import Proj, transform
import pandas as pd
import geopandas as gpd
import xarray as xr
import pygridgen as pgg
import pygridtools as pgt
import pdb
```

# projections

This projection could possibly be used as well *epsg:3857*

``` python
in_projection_string = 'epsg:4326'
out_projection_string = 'epsg:3577'
outProj = Proj(init=out_projection_string)
inProj = Proj(init=in_projection_string)
```

# boundary points

``` python
blocs  = np.array([(150.8 , -34.7),
                   (150.4 , -35.4),
                   (151.5 , -36),
                   (152.5 , -34.5),
                   (150.8 , -34.7)] )
bbetas = np.array([1,
                   1,
                   1,
                   1,
                   0])
bndry = pd.DataFrame({'beta':bbetas,
                      'longitudes':blocs[:,0],
                      'latitudes':blocs[:,1]})
bndry_pts = gpd.GeoDataFrame(bndry,
                             geometry=gpd.points_from_xy(bndry.longitudes,bndry.latitudes),
                             crs=in_projection_string)
bndry_ply = gpd.GeoDataFrame(geometry=gpd.GeoSeries([Polygon(blocs)]),
                             crs=in_projection_string)
AusAlbers_bndry_pts = bndry_pts.to_crs(out_projection_string)
AusAlbers_bndry_ply = bndry_ply.to_crs(out_projection_string)
```

# grid resolutions

``` python
nx = 100
ny = nx
```

# plot the domain

``` python
fig,ax = plt.subplots(figsize=(5, 5), subplot_kw={'aspect':'equal'})
fig    = pgt.viz.plot_domain(AusAlbers_bndry_pts, betacol='beta', ax=ax)
plt.savefig('boundary.png')
```

# generate the grid with gridgen

``` python
grid   = pgg.Gridgen(AusAlbers_bndry_pts.geometry.x,
                     AusAlbers_bndry_pts.geometry.y,
                     AusAlbers_bndry_pts.beta,
                     shape=(nx,ny),
                     ul_idx=2)
fig,ax = plt.subplots(figsize=(7, 7), subplot_kw={'aspect':'equal'})
fig    = pgt.viz.plot_cells(grid.x, grid.y, ax=ax)
plt.savefig('grid_and_boundary.png')
```

# coast mask

``` python
cst     = gpd.read_file('./AusCst/DEACoastlines_annualcoastlines_v1.0.0.shp')
jb_cst  = cst.clip(AusAlbers_bndry_ply)
fig, ax = plt.subplots(figsize=(12, 8))
jb_cst.plot(ax=ax, color="purple")
AusAlbers_bndry_ply.boundary.plot(ax=ax, color="red")
ax.set_axis_off()
plt.savefig('./coast.png')
```

# import bathymetry

``` python
bath = xr.open_dataset('./GEBCO_2021_AustraliaStation.nc')
```
