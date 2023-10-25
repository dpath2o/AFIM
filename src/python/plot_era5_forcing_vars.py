import sys
sys.path.insert(0, './AFIM/src/python')
sys.path.insert(0, '/home/581/da1339/.local/lib/python3.9/site-packages')
import os
import numpy              as np
import xarray             as xr
import matplotlib.pyplot  as plt
import cartopy.crs        as ccrs
import cmocean            as cm
import cartopy.feature    as cft
import matplotlib.path    as mpath
import metpy.calc         as mpc
from metpy.units          import units

def plot_variable(da, long_name, units, lon_name, lat_name, D_graph_base, F_name, vmin=None, vmax=None, cmap=None):
    print(f"attempting to plot: {D_graph_base}/{F_name}")
    fig, ax = plt.subplots(1,1, figsize=(15,9), dpi=150, facecolor="w", subplot_kw=dict(projection=ccrs.SouthPolarStereo()))
    theta = np.linspace(0, 2*np.pi, 100)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpath.Path(verts * radius + center)
    ax.set_extent([0,360,-90,-50], ccrs.PlateCarree())
    ax.set_boundary(circle, transform=ax.transAxes)
    ax.add_feature(cft.LAND, color='darkgrey')
    ax.add_feature(cft.COASTLINE, linewidth=.5)
    pcm = ax.pcolormesh(da[lon_name], da[lat_name], da, vmin=vmin, vmax=vmax, cmap=cmap, transform=ccrs.PlateCarree())
    plt.colorbar(pcm, label=units)
    plt.title(long_name)
    ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='gray', alpha=0.25, linestyle='--')
    if not os.path.exists(D_graph_base): os.makedirs(D_graph_base)
    plt.savefig(os.path.join(D_graph_base,F_name))
    plt.close()

era5         = xr.open_dataset("/scratch/jk72/da1339/cice-dirs/input/AFIM/forcing/0p25/ERA5/24XDAILY/era5_for_cice6_2005.CREATED_20230709.nc")
D_graph_base = "/g/data/jk72/da1339/GRAPHICAL/ERA5/0p25_CICE6_FRCG_VARS_REGRIDDED_20230709"

for var_name in era5.data_vars:
    var       = era5[var_name]
    long_name = var.attrs['long_name']
    units     = var.attrs['units']
    if var_name=='airtmp':
        vmin = 240
        vmax = 280
        cmap = cm.cm.thermal
    if var_name=='dlwsfc':
        vmin = 0
        vmax = 450
        cmap = cm.cm.amp
    if var_name=='glbrad':
        vmin = 0
        vmax = 1000
        cmap = cm.cm.amp
    if var_name=='spchmd':
        vmin = 0.005
        vmax = 0.02
        cmap = cm.cm.amp
    if var_name=='ttlpcp':
        vmin = 0
        vmax = 0.01
        cmap = cm.cm.amp
    if var_name=='wndnwd':
        continue
    D_save = os.path.join(D_graph_base,var_name)
    if not os.path.exists(D_save): os.makedirs(D_save)
    for i,t in enumerate(var["time"].values):
        if var_name=='wndewd':
            da        = mpc.wind_speed( era5['wndewd'].isel(time=i) , era5['wndnwd'].isel(time=i) )
            var_name  = 'wspd'
            D_save    = os.path.join(D_graph_base,var_name)
            if not os.path.exists(D_save): os.makedirs(D_save)
            long_name = "absolute wind speed"
            units     = "m/s"
            vmin      = 0
            vmax      = 25
            cmap      = cm.cm.speed
        else:
            da = var.isel(time=i)
        dt_str = t.astype('M8[m]').tolist().strftime('%Y-%m-%d-%H')
        F_name = f"{dt_str}_{var_name}_era5_SH.png"
        if not os.path.exists(os.path.join(D_save,F_name)):
            plot_variable(da, long_name, units, 'LON', 'LAT', D_save, F_name, vmin=vmin, vmax=vmax, cmap=cmap)