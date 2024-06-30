import os
import numpy             as np
import xarray            as xr
import seaborn           as sns
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.path   as mpath
import cmocean           as cm
import cartopy.crs       as ccrs
import cartopy.feature   as cfeature

#############################################################################################################
def southern_hemisphere_dataset_animation(self, ds, var_name='',cmap=None, G_res='025', figsize=(15,12),
                                          label_cbar='', vmin=0, vmax=1, frames=12, interval=1000):
    '''
    '''
    def animate(frame):
        cax.set_array(ds[var_name].isel(time=frame).values.flatten())
        ax.set_title("Time = " + str(ds.coords['time'].values[frame])[:13])
    land_10m = cft.NaturalEarthFeature('physical', 'land', '10m', edgecolor='black', facecolor='papayawhip', linewidth=0.5)
    geolon_t = xr.open_dataset('/g/data/ik11/grids/ocean_grid_{:s}.nc'.format(G_res)).geolon_t
    geolat_t = xr.open_dataset('/g/data/ik11/grids/ocean_grid_{:s}.nc'.format(G_res)).geolat_t
    ds       = ds.assign_coords({'geolat_t': geolat_t, 'geolon_t': geolon_t})
    ds['geolat_t'] = ds.geolat_t.swap_dims({'yt_ocean':'nj','xt_ocean':'ni'})
    ds['geolon_t'] = ds.geolon_t.swap_dims({'yt_ocean':'nj','xt_ocean':'ni'})
    plt.figure(figsize=figsize)
    fig = plt.gcf()
    ax = plt.axes(projection=ccrs.SouthPolarStereo())
    ax.set_extent([-280, 80, -80, -35], crs=ccrs.PlateCarree())
    ax.add_feature(land_10m, color=[0.8, 0.8, 0.8])
    ax.coastlines(resolution='10m')
    theta = np.linspace(0, 2*np.pi, 100)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpath.Path(verts * radius + center)
    ax.set_boundary(circle, transform=ax.transAxes)
    cax = ds.isel(time=0)[var_name].plot(x='geolon_t', y='geolat_t',
                                         transform=ccrs.PlateCarree(),
                                         vmin=vmin, vmax=vmax,
                                         cmap=cmap,
                                         cbar_kwargs = {'label': label_cbar,
                                                        'fraction': 0.03,
                                                        'aspect': 15,
                                                        'shrink': 0.7})
    ani = animation.FuncAnimation(fig,             # figure
                                  animate,         # name of the function above
                                  frames=frames,   # Could also be iterable or list
                                  interval=interval)
    HTML(ani.to_jshtml())
