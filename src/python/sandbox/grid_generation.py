import warnings
warnings.simplefilter('ignore')
import numpy as np
import matplotlib as mpl
mpl.use('cairo')
import matplotlib.pyplot as plt
from shapely.geometry import Polygon, Point, mapping
from pyproj import Proj, transform
import pandas as pd
import geopandas as gpd
import xarray as xr
import pygridgen as pgg
import pygridtools as pgt
import pdb

# projections
# epsg:3857
in_projection_string = 'epsg:4326'
out_projection_string = 'epsg:3577'
outProj = Proj(init=out_projection_string)
inProj = Proj(init=in_projection_string)

# close coast with a polygon over land
coast_extra_locs = np.array([(150.44 , -35.4),
                             (150    , -35.4),
                             (150    , -34.6),
                             (150.87 , -34.59),
                             (150.83 , -34.64),
                             (150.7  , -34.9),
                             (150.75 , -34.97),
                             (150.65 , -35),
                             (150.68 , -35.13),
                             (150.51 , -35.21),
                             (150.47 , -35.3),
                             (150.44 , -35.4)])
coast_ply        = gpd.GeoDataFrame(geometry=gpd.GeoSeries([Polygon(coast_extra_locs)]),
                                    crs=in_projection_string)
AusAlbers_coast_ply = coast_ply.to_crs(out_projection_string)

# boundary points
blocs  = np.array([(150.87 , -34.59),
                   (150.83 , -34.64),
                   (150.7  , -34.9),
                   (150.75 , -34.97),
                   (150.65 , -35),
                   (150.68 , -35.13),
                   (150.51 , -35.21),
                   (150.47 , -35.3),
                   (150.44 , -35.4),
                   (151.5 , -36),
                   (152.5 , -34.5),
                   (150.87 , -34.59)] )
bbetas = np.array([1,
                   0,
                   0,
                   0,
                   0,
                   0,
                   0,
                   0,
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

# grid resolutions
nx = 100
ny = nx

# plot the domain
fig,ax = plt.subplots(figsize=(10,12))
fig    = pgt.viz.plot_domain(AusAlbers_bndry_pts, betacol='beta', ax=ax)
plt.savefig('boundary.png')

# generate the grid with gridgen
grid   = pgg.Gridgen(AusAlbers_bndry_pts.geometry.x,
                     AusAlbers_bndry_pts.geometry.y,
                     AusAlbers_bndry_pts.beta,
                     shape=(nx,ny),
                     ul_idx=2)
fig,ax = plt.subplots(figsize=(7, 7), subplot_kw={'aspect':'equal'})
fig    = pgt.viz.plot_cells(grid.x, grid.y, ax=ax)
plt.savefig('grid_and_boundary.png')

# coastline
grid    = pgt.make_grid(domain=AusAlbers_bndry_pts,nx=nx,ny=ny,checksimplepoly=False,rawgrid=False)
# https://cmi.ga.gov.au/data-products/dea/581/dea-coastlines#access
cst     = gpd.read_file('./AusCst/DEACoastlines_annualcoastlines_v1.0.0.shp')
jb_cst  = cst.clip(AusAlbers_bndry_ply)
gshhg   = xr.open_dataset('./gshhg-gmt-2.3.7/binned_GSHHS_f.nc')
pdb.set_trace()
fig, ax = plt.subplots(figsize=(10,12))
jb_cst.plot(ax=ax, color="purple")
AusAlbers_bndry_ply.boundary.plot(ax=ax, color="red")
ax.set_axis_off()
plt.savefig('./coast.png')

# coastline masking
fig,ax = plt.subplots(figsize=(10,12))
common_opts = dict(use_existing=False)
tmp_list = []
geom = [x for x in AusAlbers_coast_ply.geometry]
tmp_list.append({'geometry':geom[0]})
geom = [x for x in jb_cst.geometry]
for g in geom:
    COORDS = mapping(g)['coordinates']
    for c in COORDS:
        c = np.array(c)
        if c.size < 6:
            continue
        lons = c[:,0]
        lats = c[:,1]
        plyg = Polygon(zip(lons, lats))
        tmp_list.append({'geometry': plyg})
mask = gpd.GeoDataFrame(tmp_list,crs=out_projection_string)
_ = ( grid.mask_centroids(inside=mask, **common_opts).plot_cells(ax=ax) )
jb_cst.plot(ax=ax, alpha=0.5, color='C0')
plt.savefig('./coast_mask_centroids.png')

# masking nodes of each grid cell
fig, ax = plt.subplots(figsize=(10,12))
common_opts = dict(use_existing=False)
_ = ( grid.mask_nodes(inside=mask, min_nodes=1, **common_opts).plot_cells(ax=ax) )
mask.plot(ax=ax, alpha=0.5)
plt.savefig('./coast_mask_nodes.png')

# bathymetry
bath = xr.open_dataset('./GEBCO_2021_AustraliaStation.nc')

