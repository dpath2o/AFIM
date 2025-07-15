import os, time, imageio, pygmt
import xarray            as xr
import pandas            as pd
import geopandas         as gpd
import numpy             as np
import matplotlib.pyplot as plt
import matplotlib.dates  as mdates
from PIL                 import Image
from pathlib             import Path
from datetime            import datetime

class SeaIcePlotter:

    def __init__(self,**kwargs):
        return

    def load_ice_shelves(self):
        """
        Load and preprocess Antarctic ice shelf polygons for PyGMT plotting.

        This method reads a shapefile containing Antarctic coastal geometries,
        filters for polygons classified as ice shelves (`POLY_TYPE == 'S'`),
        ensures valid geometry, reprojects them to WGS84 (EPSG:4326),
        and applies a zero-width buffer to clean topology issues.

        Returns
        -------
        geopandas.GeoSeries
            Cleaned and reprojected geometries of Antarctic ice shelves.

        Notes
        -----
        - The input shapefile path is read from `self.config['pygmt_dict']['P_coast_shape']`.
        - This method is typically used to overlay ice shelf boundaries in PyGMT plots.
        - The returned geometries can be passed directly to `pygmt.Figure.plot()`.

        See Also
        --------
        - self.plot_FIA_FIP_faceted : Uses this method to overlay ice shelves in map panels.
        - geopandas.read_file : For reading shapefiles.
        """
        gdf                     = gpd.read_file(self.config['pygmt_dict']['P_coast_shape'])
        shelves                 = gdf[gdf['POLY_TYPE'] == 'S']
        shelves                 = shelves[~shelves.geometry.is_empty & shelves.geometry.notnull()]
        shelves                 = shelves.to_crs("EPSG:4326")
        shelves.geometry        = shelves.geometry.buffer(0)
        return shelves.geometry

    def create_IBCSO_bath(self):
        """
        Extract and save a masked IBCSO bathymetry layer for Antarctic plotting.

        This method loads the IBCSO v2.0 dataset (as a NetCDF raster), extracts
        the seafloor depth (negative elevations), masks out land areas (positive or zero),
        and saves a cleaned version to NetCDF for use in plotting.

        The result is saved at the path specified in
        `self.config['pygmt_dict']['P_IBCSO_bath']`.

        Returns
        -------
        None

        Notes
        -----
        - Input file is assumed to be a NetCDF raster with variable `band_data`.
        - Only `band=0` is used (i.e., the main bathymetric band).
        - Output variable is named `bath` and excludes attributes and encodings for simplicity.
        - This method is typically called once to prepare the bathymetry layer before reuse.

        See Also
        --------
        - self.load_IBCSO_bath : Loads the pre-saved masked bathymetry layer.
        - https://www.ibcso.org/ : International Bathymetric Chart of the Southern Ocean.
        """
        ds               = xr.open_dataset(self.config['pygmt_dict']['P_IBCSO_bed'])
        bed              = ds.band_data.isel(band=0)
        bed_masked       = bed.where(bed < 0)
        bed_masked.name  = "bath" 
        bed_masked.attrs = {}    
        ds_out           = bed_masked.to_dataset()
        ds_out.attrs     = {}        
        ds_out.encoding  = {}     
        ds_out.to_netcdf(self.config['pygmt_dict']["P_IBCSO_bath"])

    def load_IBCSO_bath(self):
        """
        Load masked IBCSO bathymetry dataset prepared by `create_IBCSO_bath()`.

        Returns the `bath` variable from the NetCDF file specified in
        `self.config['pygmt_dict']["P_IBCSO_bath"]`.

        Returns
        -------
        xarray.DataArray
            2D array of bathymetry values (only ocean depths, in meters).
            Values are negative below sea level, NaN over land.

        Notes
        -----
        - This method assumes that `create_IBCSO_bath()` has already been run.
        - Designed for use in PyGMT background plotting or masking.

        See Also
        --------
        - self.create_IBCSO_bath : Creates this file if it doesn't exist.
        - xarray.open_dataset : Loads the bathymetry layer.
        """
        return xr.open_dataset(self.config['pygmt_dict']["P_IBCSO_bath"]).bath

    def prepare_data_for_pygmt_plot(self, da, lon_coord_name=None, lat_coord_name=None, diff_plot=False):
        """
        Prepare gridded data for PyGMT scatter-style plotting by flattening and masking arrays.

        This method extracts 2D arrays of longitude, latitude, and data values from an input
        `xarray.DataArray`, applies appropriate masking, and returns flattened arrays suitable
        for `pygmt.Figure.plot(..., style="s0.2c", ...)`.

        Parameters
        ----------
        da : xarray.DataArray
            Input data array containing 2D or 3D fields (typically persistence or difference maps).
            Must contain latitude and longitude coordinates.
        lon_coord_name : str, optional
            Name of the longitude coordinate in `da`. Defaults to `self.pygmt_dict['lon_coord_name']`.
        lat_coord_name : str, optional
            Name of the latitude coordinate in `da`. Defaults to `self.pygmt_dict['lat_coord_name']`.
        diff_plot : bool, default=False
            If True, relaxes the value filter to allow valid difference values in range [-1, 1].
            If False, restricts to strictly positive values (e.g., valid FIP from 0 to 1).

        Returns
        -------
        dict
            Dictionary containing:
            - 'data'  : 1D array of filtered values to be plotted
            - 'lon'   : 1D array of longitudes (same shape as 'data')
            - 'lat'   : 1D array of latitudes  (same shape as 'data')
            - 'lon2d' : Full 2D longitude array from the original grid
            - 'lat2d' : Full 2D latitude array from the original grid

        Notes
        -----
        - Masking ensures only finite, valid values are plotted:
            * `diff_plot=True` accepts values between -1 and 1 (inclusive)
            * `diff_plot=False` accepts strictly positive values
        - This function is designed for `pygmt.Figure.plot()` with `style="s..."` scatter symbols,
          not for gridded plotting methods like `pygmt.Figure.grdimage()`.

        See Also
        --------
        - pygmt.Figure.plot : Used to plot flattened (x, y, value) tuples
        - pygmt.makecpt     : For setting the colormap before plotting
        - self.plot_FIA_FIP_faceted : Calls this method when plotting FIP panels
        """
        lon_coord_name = lon_coord_name if lon_coord_name is not None else self.pygmt_dict['lon_coord_name']
        lat_coord_name = lat_coord_name if lat_coord_name is not None else self.pygmt_dict['lat_coord_name']
        data_dict      = {}
        data2d         = np.asarray(da.data).astype(float)
        lon2d          = da[lon_coord_name].data
        lat2d          = da[lat_coord_name].data
        if diff_plot:
            mask = (data2d >= -1) & (data2d <= 1) & np.isfinite(data2d)
        else:
            mask = (data2d > 0) & np.isfinite(data2d)
        data_dict['data']  = data2d[mask].ravel()
        data_dict['lon']   = lon2d[mask].ravel()
        data_dict['lat']   = lat2d[mask].ravel()
        data_dict['lon2d'] = lon2d
        data_dict['lat2d'] = lat2d
        return data_dict

    def load_GI_lon_lats(self, data_dict):
        """
        Extract longitude and latitude positions of grounded iceberg (GI) grid cells.

        This method identifies cells where the landmask has been modified to add grounded icebergs
        (i.e., locations where the original landmask had ocean (`kmt_org == 1`) and the modified
        landmask has land (`kmt_mod == 0`)), and returns their corresponding geographic coordinates.

        Parameters
        ----------
        data_dict : dict
            Dictionary containing 'lon2d' and 'lat2d' 2D arrays (e.g., from `prepare_data_for_pygmt_plot`).

        Returns
        -------
        dict
            Dictionary with:
            - 'lon': 1D array of longitudes of grounded iceberg grid cells
            - 'lat': 1D array of latitudes  of grounded iceberg grid cells

        Notes
        -----
        - `self.P_KMT_org` and `self.P_KMT_mod` must be paths to NetCDF files with the `kmt` landmask variable.
        - `self.hemisphere_dict['nj_slice']` is used to subset the hemisphere-specific grid region.
        - The output is suitable for symbol plotting (e.g., `pygmt.Figure.plot(...)` with `style="c0.05c"`).
        """
        GI_loc_dict        = {}
        kmt_mod            = xr.open_dataset(self.P_KMT_mod).isel(nj=self.hemisphere_dict['nj_slice']).kmt.data
        kmt_org            = xr.open_dataset(self.P_KMT_org).isel(nj=self.hemisphere_dict['nj_slice']).kmt.data
        GI_mask            = (kmt_org == 1) & (kmt_mod == 0)
        GI_loc_dict['lon'] = data_dict['lon2d'][GI_mask].ravel()
        GI_loc_dict['lat'] = data_dict['lat2d'][GI_mask].ravel()
        return GI_loc_dict

    def create_cbar_frame(self, series, label, units=None, extend_cbar=False, max_ann_steps=10):
        """
        Construct a GMT-style colorbar annotation frame string for PyGMT plotting.

        This utility generates a clean, readable colorbar frame string using adaptive
        step sizing based on the data range. An optional second axis label (e.g., for units)
        and extension arrows can be included.

        Parameters
        ----------
        series : list or tuple of float
            Data range for the colorbar as [vmin, vmax].
        label : str
            Label text for the colorbar (e.g., "Fast Ice Persistence").
        units : str, optional
            Units label to be shown along the secondary (y) axis (e.g., "1/100").
        extend_cbar : bool, optional
            Whether to append extension arrows to the colorbar (+e). Default is False.
        max_ann_steps : int, default=10
            Desired maximum number of major annotations (controls tick spacing).

        Returns
        -------
        str or list of str
            GMT-format colorbar annotation string (e.g., "a0.1f0.02+lLabel"), or
            a list of two strings if `units` is provided.

        Notes
        -----
        - Tick spacing is determined using a scaled logarithmic rounding to ensure clean steps (e.g., 0.1, 0.2, 0.5).
        - If `extend_cbar` is True, `+e` is added to indicate out-of-bounds extension arrows.
        - Compatible with `pygmt.Figure.colorbar(frame=...)`.

        Examples
        --------
        >>> self.create_cbar_frame([0.01, 1.0], "Persistence", units="1/100")
        ['a0.1f0.02+lPersistence', 'y+l 1/100']
        """
        vmin, vmax = series[0], series[1]
        vrange = vmax - vmin
        raw_step = vrange / max_ann_steps
        exp = np.floor(np.log10(raw_step))
        base = 10 ** exp
        mult = raw_step / base
        if mult < 1.5:
            ann_step = base * 1
        elif mult < 3:
            ann_step = base * 2
        elif mult < 7:
            ann_step = base * 5
        else:
            ann_step = base * 10
        tick_step = ann_step / 5
        ann_str  = f"{ann_step:.3f}".rstrip("0").rstrip(".")
        tick_str = f"{tick_step:.3f}".rstrip("0").rstrip(".")
        # Build annotation string
        frame = f"a{ann_str}f{tick_str}+l{label}"
        if extend_cbar:
            frame += "+e"  # or use "+eU" / "+eL" for one-sided arrows
        if units is not None: 
            return [frame, f"y+l {units}"]
        else:
            return frame

    def get_meridian_center_from_geographic_extent(self, geographic_extent):
        """
        Determine the optimal central meridian for PyGMT polar stereographic projections.

        Given a geographic extent in longitude/latitude format, this method calculates
        the central meridian (longitude) to use in 'S<lon>/<lat>/<width>' PyGMT projection
        strings. It accounts for dateline wrapping and ensures the plot is centered visually.

        Parameters
        ----------
        geographic_extent : list of float
            Geographic region as [min_lon, max_lon, min_lat, max_lat]. Accepts longitudes in
            either [-180, 180] or [0, 360] and gracefully handles dateline crossing.

        Returns
        -------
        float
            Central meridian (longitude) in the range [-180, 180].

        Notes
        -----
        - If the computed center falls outside the intended range (e.g., due to dateline wrapping),
          the method rotates the meridian 180° to better align the figure.
        - The result is saved to `self.plot_meridian_center` for reuse.
        - Used in PyGMT stereographic projections like `'S{lon}/-90/30c'`.

        See Also
        --------
        - https://docs.generic-mapping-tools.org/latest/cookbook/proj.html
        - pygmt.Figure.basemap : For setting map projections using meridian centers.
        """
        lon_min, lon_max = geographic_extent[0], geographic_extent[1]
        lon_min_360 = lon_min % 360
        lon_max_360 = lon_max % 360
        if (lon_max_360 - lon_min_360) % 360 > 180:
            center = ((lon_min_360 + lon_max_360 + 360) / 2) % 360
        else:
            center = (lon_min_360 + lon_max_360) / 2
        if center > 180:
            center -= 360
        # Edge case fix: ensure center is visually aligned with geographic_extent
        # If the computed center is 180° out of phase (i.e., upside-down plots)
        if not (geographic_extent[0] <= center <= geographic_extent[1]):
            # Flip 180°
            center = (center + 180) % 360
            if center > 180:
                center -= 360
        #print(f"meridian center computed as {center:.2f}°")
        self.plot_meridian_center = center
        return center
        
    def plot_FIA_FIP_faceted(self, FIA_dict, FIP_DA,
                              sim_name       = None,
                              dt_range_str   = None,
                              P_png          = None,
                              enable_FIA     = True,
                              enable_FIP     = True,
                              plot_GI        = False,
                              GI_fill_color  = "red",
                              GI_sq_size     = 0.05,
                              diff_plot      = False,
                              FIA_ylim       = (0,1000),
                              roll_days      = 0,
                              lon_coord_name = None,
                              lat_coord_name = None,
                              cmap           = None,
                              series         = None,
                              cbar_frame     = None,
                              overwrite_fig  = None,
                              show_fig       = None):
        """
        Generate a composite figure showing seasonal fast ice area climatology (FIA) and fast ice persistence (FIP) maps.

        This method produces a multi-panel PyGMT figure containing:
        - A top panel showing the daily climatology of modeled and observed FIA across a full seasonal cycle.
        - Lower panels showing faceted FIP maps by Antarctic sector using symbol plots or difference plots.
        - Optional overlays of grounded iceberg locations.

        Parameters
        ----------
        FIA_dict : dict of xarray.DataArray or Dataset
            Dictionary containing fast ice area time series from one or more simulations or methods. Keys define line labels.
        FIP_DA : xarray.DataArray
            2D map of fast ice persistence (0–1), typically from seasonal averaging over time.
        sim_name : str, optional
            Simulation name used for labeling and figure output. Defaults to `self.sim_name`.
        dt_range_str : str, optional
            String describing the date range (e.g., "1993–1999") for filename or annotation. Auto-handled if None.
        P_png : pathlib.Path, optional
            Path to save the PNG figure. If not provided, figure is not saved.
        enable_FIA : bool, default=True
            Whether to plot the FIA climatology time series panel.
        enable_FIP : bool, default=True
            Whether to plot FIP spatial map panels.
        plot_GI : bool, default=False
            Whether to overlay grounded iceberg regions on the map panels.
        GI_fill_color : str, default="red"
            Color to use for grounded iceberg overlay.
        GI_sq_size : float, default=0.05
            Size of the circle or square marker used to plot grounded iceberg locations (in cm).
        diff_plot : bool, default=False
            If True, plot differences between model and observational FIA rather than raw values.
        FIA_ylim : tuple of float, default=(100, 1000)
            Y-axis range for the FIA climatology subplot (in 1000 km²).
        roll_days : int, optional
            Rolling window applied before plotting (currently unused but included for extensibility).
        lon_coord_name : str, optional
            Name of longitude coordinate in FIP_DA. If None, inferred.
        lat_coord_name : str, optional
            Name of latitude coordinate in FIP_DA. If None, inferred.
        cmap : str or PyGMT colormap, optional
            Color palette table (CPT) used for the FIP map. Defaults to `self.pygmt_dict["FIP_CPT"]`.
        series : list, optional
            Range and increment for color scaling (e.g., [0.01, 1.0, 0.01]).
        cbar_frame : list of str, optional
            List defining colorbar annotations. Default includes axis labels for FIP units.
        overwrite_fig : bool, optional
            If True, overwrite existing PNG if `P_png` exists. Defaults to `self.ow_fig`.
        show_fig : bool, optional
            If True, show the PyGMT figure interactively. Defaults to `self.show_fig`.

        Returns
        -------
        None
            This method creates and saves/shows a figure but does not return an object.

        Notes
        -----
        - FIA panel shows observed AF2020 climatology as dashed blue line, and modeled climatologies as filled bands (min–max) and solid mean lines.
        - If `diff_plot=True`, plots the difference between model mean and AF2020 climatology.
        - FIP panels are faceted by Antarctic regions as defined in `self.Ant_2sectors`, using small symbol plots.
        - Grounded iceberg data is loaded via `load_GI_lon_lats()` if enabled.
        - All PyGMT styling options are drawn from `self.pygmt_dict` for consistency across figures.
        - A helper text file is generated (`xannots.txt`) to customize month-based x-axis ticks.

        See Also
        --------
        - self.load_AF2020_FIA_summary
        - self.prepare_data_for_pygmt_plot
        - self.load_GI_lon_lats
        - self.Ant_2sectors
        """
        sim_name     = sim_name      if sim_name      is not None else self.sim_name
        show_fig     = show_fig      if show_fig      is not None else self.show_fig
        ow_fig       = overwrite_fig if overwrite_fig is not None else self.ow_fig
        cmap         = cmap          if cmap          is not None else self.pygmt_dict.get("FIP_CPT")
        series       = series        if series        is not None else [0.01, 1.0, 0.01]
        cbar_frame   = cbar_frame    if cbar_frame    is not None else ["x+lFast Ice Persistence", "y+l1/100"]
        AF2020_CSV   = self.load_AF2020_FIA_summary( start=self.dt0_str, end=self.dtN_str )
        obs_FIA      = AF2020_CSV['FIA_clim_repeat'].sel(region='circumpolar')
        xticks       = [1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 360]
        month_labels = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec", "Dec"]
        xannot_path  = Path("xannots.txt")
        xannot_lines = [f"{tick}\tig\t{label}" for tick, label in zip(xticks, month_labels)]
        xannot_path.write_text("\n".join(xannot_lines) + "\n")
        default_line_colors = ["orange", "green", "blue", "red", "magenta", "cyan"]
        fig = pygmt.Figure()
        pygmt.config(FORMAT_GEO_MAP = "ddd.x",
                     MAP_GRID_PEN   = "0p,white")
        if diff_plot:
            yaxis_lab = f"ya10f+lFast Ice Area Difference (obs-sim; 1000 km\u00b2)"
        else:
            yaxis_lab = f"ya50f+lFast Ice Area (1000 km\u00b2)"
        fig.basemap(region     = [1, 360, FIA_ylim[0], FIA_ylim[1]],
                    projection = "X30c/15c",
                    **{"frame" : [f"sxc{xannot_path}",
                                  f"ya100f+lFast Ice Area (1000 km\u00b2)",
                                  "WSne"]})
        if enable_FIA:
            for i, (name, da) in enumerate(FIA_dict.items()):
                if isinstance(da, xr.Dataset):
                    da = da.to_array().squeeze()
                da         = da.load()
                df         = pd.DataFrame({"time": pd.to_datetime(da["time"].values),
                                           "area": da.values})
                df["doy"]  = df["time"].dt.dayofyear
                df["year"] = df["time"].dt.year
                grouped    = df.groupby("doy")["area"]
                area_min   = grouped.min()
                area_max   = grouped.max()
                area_mean  = grouped.mean()
                if i == 0 and not diff_plot:
                    obs_df        = pd.DataFrame({"time": pd.to_datetime(obs_FIA["time"].values),
                                                  "area": obs_FIA.values})
                    obs_df["doy"] = obs_df["time"].dt.dayofyear
                    fig.plot(x=obs_df["doy"], y=obs_df["area"], pen="1.5p,blue,--", label="AF2020 Climatology (2000–2018)")
                if name in self.plot_ice_area_dict:
                    leg_lab    = self.plot_ice_area_dict[name]["label"].format(min_win=self.bool_min_days, win=self.bool_window)
                    line_color = self.plot_ice_area_dict[name]["color"]
                else:
                    leg_lab    = name
                    line_color = default_line_colors[i % len(default_line_colors)]
                if diff_plot:
                    obs_interp     = np.interp(area_mean.index, obs_df["doy"], obs_df["area"])
                    area_mean_diff = obs_interp - area_mean.values
                    fig.plot(x     = area_mean.index,
                             y     = area_mean_diff,
                             pen   = f"2p,{line_color}",
                             label = f"{leg_lab} diff")
                else:
                    fig.plot(x            = np.concatenate([area_min.index, area_max.index[::-1]]),
                             y            = np.concatenate([area_min.values, area_max.values[::-1]]),
                             fill         = f"{line_color}@70",
                             close        = True,
                             transparency = 60)
                    fig.plot(x     = area_mean.index,
                             y     = area_mean.values,
                             pen   = f"2p,{line_color}",
                             label = f"{leg_lab} climatology")
            fig.legend(position="JTL+jTL+o0.2c", box="+gwhite+p.5p")
        if enable_FIP:
            fig.shift_origin(yshift="-6c")
            ANT_IS         = self.load_ice_shelves()
            plot_data_dict = self.prepare_data_for_pygmt_plot(FIP_DA, lon_coord_name=lon_coord_name, lat_coord_name=lat_coord_name, diff_plot=diff_plot)
            if plot_GI:
                plot_GI_dict = self.load_GI_lon_lats(plot_data_dict)
            pygmt.makecpt(cmap=cmap, series=series)
            for i, (reg_name, reg_vals) in enumerate(self.Ant_2sectors.items()):
                region     = reg_vals['plot_region']
                projection = reg_vals['projection'].format(fig_width=30)
                if i>0:
                    fig.shift_origin(yshift="-10.25c")
                fig.basemap(region=region, projection=projection, frame=['xa10g10f','ya5g5f'])
                fig.coast(land='gray', water='white', shorelines="1/0.5p,gray30")
                fig.plot(x=plot_data_dict['lon'], y=plot_data_dict['lat'], fill=plot_data_dict['data'], style="s0.2c", cmap=True)
                if plot_GI:
                    fig.plot(x=plot_GI_dict['lon'], y=plot_GI_dict['lat'], fill=GI_fill_color, style=f"c{GI_sq_size}c")
                fig.plot(data=ANT_IS, pen="0.2p,gray", fill="lightgray")
            fig.colorbar(position="JBC+w20c/1c+mc+h", frame=cbar_frame)
        if P_png:
            if not P_png.exists():
                P_png.parent.mkdir(parents=True, exist_ok=True)
                fig.savefig(P_png)  
                self.logger.info(f"Saved figure to {P_png}")
            else:
                if ow_fig:
                    fig.savefig(P_png)  
                    self.logger.info(f"Saved figure to {P_png}")
                else:
                    self.logger.info(f"{P_png} already exists and not overwriting")
        if show_fig:
            fig.show()
        pygmt.clib.Session.__exit__

    def pygmt_map_plot_one_var(self, da, var_name,
                               sim_name       = None,
                               plot_regions   = None,
                               regional_dict  = None,
                               hemisphere     = "south",
                               time_stamp     = None,
                               tit_str        = None,
                               plot_GI        = False,
                               diff_plot      = False,
                               cmap           = None,
                               series         = None,
                               reverse        = None,
                               cbar_label     = None,
                               cbar_units     = None,
                               extend_cbar    = False,
                               cbar_position  = None,
                               lon_coord_name = None,
                               lat_coord_name = None,
                               fig_size       = None,
                               var_sq_size    = 0.2,
                               GI_sq_size     = 0.1,
                               GI_fill_color  = "red",
                               plot_iceshelves= True,
                               plot_bathymetry= True,
                               add_stat_annot = False,
                               land_color     = None,
                               water_color    = None,
                               P_png          = None,
                               var_out        = None,
                               overwrite_fig  = None,
                               show_fig       = None):
        """
        Generate a PyGMT figure showing a variable (e.g., FIA, FIP, differences) as a spatial map
        over Antarctic regions or hemispheric view, optionally including grounded icebergs,
        ice shelf outlines, and bathymetry.

        This flexible mapping function supports multiple types of visualizations:
        - Scalar field plots (e.g., fast ice persistence, concentration)
        - Binary or categorical masks (e.g., agreement maps, simulation masks)
        - Difference plots (e.g., observation minus simulation)

        Parameters
        ----------
        da : xarray.DataArray
            Input 2D (or broadcastable) data array to be plotted, e.g., fast ice persistence.
        var_name : str
            Name of the variable to be plotted. Used to look up color map settings and labels.
        sim_name : str, optional
            Simulation name used for file naming and figure annotation. Defaults to `self.sim_name`.
        plot_regions : int or None, optional
            Number of regional plots:
            - 8 : Antarctic 8-sector view
            - 2 : East/West sectors
            - None : Full hemisphere plot (default)
        regional_dict : dict, optional
            Custom dictionary of plotting regions (overrides built-ins if provided).
        hemisphere : str, default="south"
            Hemisphere name for hemispheric plot. Typically "south".
        time_stamp : str, optional
            Timestamp string for file naming. Defaults to `self.dt0_str`.
        tit_str : str, optional
            Title string to display on the figure.
        plot_GI : bool, default=False
            Whether to overlay grounded iceberg locations using `load_GI_lon_lats()`.
        diff_plot : bool, default=False
            Whether this is a difference plot (used for labeling, masking, and symbology).
        cmap : str, optional
            Color map name for continuous scalar plots.
        series : list of float, optional
            Min, max, and increment for the colorbar (e.g., [0, 1, 0.1]).
        reverse : bool, optional
            Whether to reverse the colormap.
        cbar_label : str, optional
            Label to show on the colorbar.
        cbar_units : str, optional
            Units for the colorbar (shown on secondary axis).
        extend_cbar : bool, default=False
            Whether to add extension arrows to colorbar.
        cbar_position : str, optional
            PyGMT-compatible position string for placing the colorbar.
        lon_coord_name : str, optional
            Longitude coordinate name in `da`. Defaults to `self.pygmt_dict`.
        lat_coord_name : str, optional
            Latitude coordinate name in `da`. Defaults to `self.pygmt_dict`.
        fig_size : float, optional
            Figure width in centimeters for hemispheric plot.
        var_sq_size : float, default=0.2
            Marker size (in cm) for main plotted variable (for scatter plotting).
        GI_sq_size : float, default=0.1
            Marker size (in cm) for grounded iceberg overlay.
        GI_fill_color : str, default="red"
            Color to use for grounded iceberg markers.
        plot_iceshelves : bool, default=True
            Whether to overlay Antarctic ice shelf outlines.
        plot_bathymetry : bool, default=True
            Whether to plot IBCSO bathymetry using shaded relief (`grdimage`).
        add_stat_annot : bool, default=False
            Whether to annotate figure with basic regional statistics.
        land_color : str, optional
            Color for land. Defaults to `self.pygmt_dict['land_color']`.
        water_color : str, optional
            Color for ocean/water. Defaults to `self.pygmt_dict['water_color']`.
        P_png : pathlib.Path, optional
            File path to save the figure. If None, path is generated automatically.
        var_out : str, optional
            Output variable name used in file naming. Defaults to `var_name`.
        overwrite_fig : bool, optional
            Whether to overwrite existing figure if it exists.
        show_fig : bool, optional
            Whether to display the figure interactively.

        Returns
        -------
        None
            The method generates and optionally saves or displays a PyGMT figure.

        Notes
        -----
        - Supports 8-region, 2-region, or full hemisphere plotting via `plot_regions`.
        - For "diff" plots (where 'diff' in var_name), colors are assigned categorically (e.g., agreement/simulation/observation).
        - Bathymetry and ice shelf overlays are loaded from NetCDF and shapefiles respectively.
        - The method gracefully skips plotting if required coordinates or data are missing.

        See Also
        --------
        - self.prepare_data_for_pygmt_plot : Prepares data dictionary for plotting
        - self.create_cbar_frame           : Builds formatted colorbar strings
        - self.load_GI_lon_lats            : Loads grounded iceberg locations
        - self.load_ice_shelves            : Loads Antarctic ice shelf polygons
        - self.load_IBCSO_bath             : Loads bathymetry grid from IBCSO

        Examples
        --------
        >>> self.pygmt_map_plot_one_var(FIP_DA, "FIP", plot_regions=8, show_fig=True)
        >>> self.pygmt_map_plot_one_var(diff_DA, "FIA_diff", diff_plot=True, cmap="categorical", plot_GI=True)
        """
        sim_name    = sim_name      if sim_name      is not None else self.sim_name
        show_fig    = show_fig      if show_fig      is not None else self.show_fig
        ow_fig      = overwrite_fig if overwrite_fig is not None else self.ow_fig
        time_stamp  = time_stamp    if time_stamp    is not None else self.dt0_str
        lon_coord_name = lon_coord_name if lon_coord_name is not None else self.pygmt_dict.get("lon_coord_name", "TLON")
        lat_coord_name = lat_coord_name if lat_coord_name is not None else self.pygmt_dict.get("lat_coord_name", "TLAT")
        cmap        = cmap          if cmap          is not None else self.plot_var_dict[var_name]['cmap']
        series      = series        if series        is not None else self.plot_var_dict[var_name]['series']
        reverse     = reverse       if reverse       is not None else self.plot_var_dict[var_name]['reverse']
        cbar_lab    = cbar_label    if cbar_label    is not None else self.plot_var_dict[var_name]['name']
        cbar_units  = cbar_units    if cbar_units    is not None else self.plot_var_dict[var_name]['units']
        fig_size    = fig_size      if fig_size      is not None else self.pygmt_dict['fig_size']
        cbar_pos    = cbar_position if cbar_position is not None else self.pygmt_dict['cbar_pos'].format(width=fig_size*0.8,height=0.75)
        land_color  = land_color    if land_color    is not None else self.pygmt_dict['land_color']
        water_color = water_color   if water_color   is not None else self.pygmt_dict['water_color']
        if var_out is None:
            var_out = var_name
        if plot_iceshelves:
            ANT_IS = self.load_ice_shelves()
        if plot_bathymetry:
            SO_BATH = self.load_IBCSO_bath()
        plot_data_dict = self.prepare_data_for_pygmt_plot(da, lon_coord_name=lon_coord_name, lat_coord_name=lat_coord_name, diff_plot=diff_plot)
        required_keys = ['lon', 'lat', 'data']
        try:
            if not isinstance(plot_data_dict, dict):
                self.logger.warning("plot_data_dict is not a dictionary — skipping plot.")
                return
            for k in required_keys:
                if k not in plot_data_dict:
                    self.logger.warning(f"Missing key '{k}' in plot_data_dict — skipping plot.")
                    return
                v = plot_data_dict[k]
                if v is None:
                    self.logger.warning(f"plot_data_dict['{k}'] is None — skipping plot.")
                    return
                if hasattr(v, "size") and v.size == 0:
                    self.logger.warning(f"plot_data_dict['{k}'] is empty — skipping plot.")
                    return
        except Exception as e:
            self.logger.warning(f"Skipping plot due to error: {e}")
            return
        cbar_frame = self.create_cbar_frame(series, cbar_lab, units=cbar_units, extend_cbar=extend_cbar)
        hem_plot   = False
        if plot_GI:
                plot_GI_dict = self.load_GI_lon_lats(plot_data_dict)
        if plot_regions is not None and plot_regions==8:
            self.logger.info("method will plot eight Antarctic sectors regional dictionary")
            reg_dict = self.Ant_8sectors
        elif plot_regions is not None and plot_regions==2:
            self.logger.info("method will plot two Antarctic sectors regional dictionary")
            reg_dict = self.Ant_2sectors
        elif plot_regions is not None and regional_dict is not None:
            self.logger.info("method will plot regional dictionary passed to this method")
            reg_dict = regional_dict
        elif plot_regions is not None and regional_dict is None:
            self.logger.info("plot_regions argument not valid")
        else:
            self.logger.info("method will plot hemispheric data")
            hem_plot = True
            reg_dict = self.hemispheres_dict
        if tit_str is not None:
            basemap_frame = ["af", f"+t{tit_str}"]
        else:
            basemap_frame = ["af"]
        for i, (reg_name, reg_vals) in enumerate(reg_dict.items()):
            if hem_plot and reg_name!=hemisphere:
                continue
            if P_png is None and self.save_fig:
                P_png = Path(self.D_graph, sim_name, reg_name, var_out, f"{time_stamp}_{sim_name}_{reg_name}_{var_out}.png")
            region     = reg_vals['plot_region']
            projection = reg_vals['projection']               
            if hem_plot:
                projection = projection.format(fig_size=fig_size)
            elif reg_name in list(self.Ant_8sectors.keys()):
                MC         = self.get_meridian_center_from_geographic_extent(region)
                projection = projection.format(MC=MC, fig_size=fig_size)
            elif reg_name in list(self.Ant_2sectors.keys()):
                projection = projection.format(fig_width=fig_size)
            fig = pygmt.Figure()
            pygmt.config(FONT_TITLE         = "16p,Courier-Bold",
                         FONT_ANNOT_PRIMARY = "14p,Helvetica",
                         COLOR_FOREGROUND   = 'black')
            fig.basemap(region=region, projection=projection, frame=basemap_frame)
            if plot_bathymetry:
                fig.grdimage(grid=SO_BATH, cmap='geo')
            else:
                fig.coast(region=region, projection=projection, shorelines="1/0.5p,gray30", land=land_color, water=water_color)
            if "diff" in var_name.lower():
                lat        = da[lat_coord_name].values.flatten()
                lon        = da[lon_coord_name].values.flatten()
                val        = da.values.flatten()
                valid_mask = ~np.isnan(val)
                lat_valid  = lat[valid_mask]
                lon_valid  = lon[valid_mask]
                val_valid  = val[valid_mask].astype(int)
                label_map  = {1: "simulation", 0: "agreement", 2: "observation"}
                labels     = [label_map[v] for v in val_valid]
                df         = pd.DataFrame({"longitude" : lon_valid,
                                           "latitude"  : lat_valid,
                                           "z"         : val_valid.astype(int)})
                pygmt.makecpt(cmap="categorical", series=[0,2,1], color_model="+cagreement,simulation,observation")
                fig.plot(data=df, style=f"s{var_sq_size}c", cmap=True)
            elif "mask" in var_name.lower():
                fig.plot(x=plot_data_dict['lon'], y=plot_data_dict['lat'], fill='red', style=f"s{var_sq_size}c")
            else:
                pygmt.makecpt(cmap=cmap, reverse=reverse, series=series)
                fig.plot(x=plot_data_dict['lon'], y=plot_data_dict['lat'], fill=plot_data_dict['data'], style=f"s{var_sq_size}c", cmap=True)           
            if plot_bathymetry:
                fig.coast(region=region, projection=projection, shorelines="1/0.5p,gray30")
            if plot_GI:
                fig.plot(x=plot_GI_dict['lon'], y=plot_GI_dict['lat'], fill=GI_fill_color, style=f"c{GI_sq_size}c")
            if plot_iceshelves:
                fig.plot(data=ANT_IS, pen="0.2p,gray", fill="lightgray")
            if add_stat_annot:
                annot_text = self.generate_regional_annotation_stats(da, region, lon_coord_name, lat_coord_name)
                for i, line in enumerate(annot_text):
                    try:
                        fig.text(position="TR", text=line, font="12p,Helvetica-Bold,black", justify="LM", no_clip=True, offset=f"-1/{-0.5*i}")
                    except pygmt.exceptions.GMTCLibError as e:
                        self.logger.warning(f"Error in plotting anotation text {e} -- skipping annotation")
            if "diff" in var_name.lower():
                fig.colorbar(position=cbar_pos, frame=["x+l" + cbar_lab])
            elif "mask" not in var_name.lower():
                try:
                    fig.colorbar(position=cbar_pos, frame=cbar_frame)
                except pygmt.exceptions.GMTCLibError as e:
                    self.logger.warning(f"Error in adding colorbar: {e} — skipping colorbar.")
            if P_png:
                if not P_png.exists():
                    P_png.parent.mkdir(parents=True, exist_ok=True)
                    fig.savefig(P_png)  
                    self.logger.info(f"Saved figure to {P_png}")
                else:   
                    if ow_fig:
                        fig.savefig(P_png)  
                        self.logger.info(f"Saved figure to {P_png}")
                    else:
                        self.logger.info(f"{P_png} already exists and not overwriting")
                P_png = None
            if show_fig:
                fig.show()
            pygmt.clib.Session.__exit__

    def generate_regional_annotation_stats(self, da, region, lon_coord_name, lat_coord_name):
        """
        Generate summary statistics (mean, std, min, max) and their locations within a geographic region.

        This helper function extracts the portion of a 2D `xarray.DataArray` that falls within a given 
        geographic bounding box, computes spatial statistics on valid (non-NaN) values, and formats 
        the results for annotation in PyGMT plots.

        Parameters
        ----------
        da : xarray.DataArray
            2D gridded data array (e.g., fast ice persistence) with latitude and longitude coordinates.
        region : list or tuple of float
            Bounding box for the region in the form [lon_min, lon_max, lat_min, lat_max].
        lon_coord_name : str
            Name of the longitude coordinate in `da` (e.g., 'TLON').
        lat_coord_name : str
            Name of the latitude coordinate in `da` (e.g., 'TLAT').

        Returns
        -------
        list of str
            Formatted text lines with:
            - Mean
            - Standard deviation
            - Data extent (number of valid grid cells)
            - Minimum value and its lat/lon location
            - Maximum value and its lat/lon location

        Notes
        -----
        - Longitude and latitude slices are derived using `np.searchsorted`, assuming monotonic grid coordinates.
        - Only finite (non-NaN) values are included in statistics.
        - Output is intended to be used directly in `fig.text()` or `fig.legend()` in PyGMT plotting routines.

        Examples
        --------
        >>> stats = self.generate_regional_annotation_stats(FIP_DA, region=[0, 90, -75, -60], lon_coord_name='TLON', lat_coord_name='TLAT')
        >>> for line in stats:
        >>>     print(line)
        Mean: 0.73
        Std: 0.18
        Extent: 1246
        Min:  0.02 at (-68.41, 23.65)
        Max: 1.00 at (-66.82, 45.12)
        """
        # Get the latitude and longitude coordinates (TLAT, TLON)
        da_lats    = da[lat_coord_name].values
        da_lons    = da[lon_coord_name].values
        lon_min_ni = np.searchsorted(da_lons[0, :], region[0])
        lon_max_ni = np.searchsorted(da_lons[0, :], region[1])        
        lat_min_nj = np.searchsorted(da_lats[:, 0], region[2])
        lat_max_nj = np.searchsorted(da_lats[:, 0], region[3])
        da_sliced  = da.isel(nj=slice(lat_min_nj, lat_max_nj), ni=slice(lon_min_ni, lon_max_ni))
        lon_sliced = da[lon_coord_name].isel(nj=slice(lat_min_nj, lat_max_nj), ni=slice(lon_min_ni, lon_max_ni))
        lat_sliced = da[lat_coord_name].isel(nj=slice(lat_min_nj, lat_max_nj), ni=slice(lon_min_ni, lon_max_ni))
        # Get the data, latitudes, and longitudes
        data       = da_sliced.values.flatten()
        lat        = lat_sliced.values.flatten()
        lon        = lon_sliced.values.flatten()
        valid_mask = ~np.isnan(data)
        data_valid = data[valid_mask]
        lat_valid  = lat[valid_mask]
        lon_valid  = lon[valid_mask]
        # Calculate statistics
        spatial_mean   = np.mean(data_valid)
        spatial_std    = np.std(data_valid)
        spatial_extent = len(data_valid)
        data_min       = np.min(data_valid)
        data_max       = np.max(data_valid)
        # Get min/max locations (lat, lon)
        min_index    = np.argmin(data_valid)
        max_index    = np.argmax(data_valid)
        min_location = (lat_valid[min_index], lon_valid[min_index])
        max_location = (lat_valid[max_index], lon_valid[max_index])
        # Format text to display
        text  = [f"Mean: {spatial_mean:.2f}", f"Std: {spatial_std:.2f}", f"Extent: {spatial_extent}"]
        text += [f"Min:  {data_min:.2f} at {min_location}", f"Max: {data_max:.2f} at {max_location}"]
        return text

    def plot_monthly_ice_metric_by_year(self, area_dict,
                                        ice_type         = "FIA",
                                        roll_days        = 0,
                                        ylim             = None,
                                        figsize          = (18, 10),
                                        tit_str          = None,
                                        P_png            = None,
                                        tick_fontsize    = 12,
                                        label_fontsize   = 14,
                                        title_fontsize   = 18,
                                        legend_fontsize  = 12,
                                        plot_annotations = False):
        """
        Plot monthly climatology of an ice area or volume metric by year and simulation.

        This method produces a seasonal line plot of a sea ice metric (e.g., FIA, SIA, FIV) 
        for multiple simulations or observational datasets. It computes daily climatologies 
        across years and shows min/max shading and mean lines. Optional annotations mark 
        annual maxima and minima.

        Parameters
        ----------
        area_dict : dict of xarray.DataArray
            Dictionary of time series datasets keyed by source name (e.g., model or obs name).
            Each entry must include a 'time' coordinate.
        ice_type : str, default="FIA"
            Type of ice metric to plot. Options include:
            - "FIA": Fast Ice Area (10³ km²)
            - "FIV": Fast Ice Volume (km³)
            - "SIA": Sea Ice Area (10⁶ km²)
            - "SIV": Sea Ice Volume (10⁶ km³)
        roll_days : int, default=0
            Number of days for centered rolling mean smoothing. Set to 0 to disable smoothing.
        ylim : tuple of float, optional
            Y-axis limits (e.g., `(100, 800)`). Defaults vary by `ice_type`.
        figsize : tuple of float, default=(18, 10)
            Figure size in inches.
        tit_str : str, optional
            Title to display above the plot.
        P_png : pathlib.Path, optional
            File path to save the figure (PNG). Saved only if `self.save_fig` is True.
        tick_fontsize : int, default=12
            Font size for x/y ticks.
        label_fontsize : int, default=14
            Font size for axis labels.
        title_fontsize : int, default=18
            Font size for the plot title.
        legend_fontsize : int, default=12
            Font size for the legend.
        plot_annotations : bool, default=False
            If True, adds max/min DOY labels and markers for each simulation.

        Returns
        -------
        None
            A Matplotlib figure is shown and optionally saved.

        Notes
        -----
        - For each dataset, the function groups daily values by DOY across all years (except the first year),
          computes mean/min/max, and plots them with fill and line styles.
        - Observational climatologies ("AF2020db_cli", "NSIDC") are treated separately and drawn as dashed lines.
        - When `plot_annotations=True`, the DOY of annual maxima and minima is printed and plotted.
        - The x-axis is labeled with month names, and the y-axis label depends on the `ice_type`.

        Examples
        --------
        >>> self.plot_monthly_ice_metric_by_year(FIA_dict, ice_type="FIA", roll_days=15, tit_str="Fast Ice Area Climatology")
        >>> self.plot_monthly_ice_metric_by_year(SIA_dict, ice_type="SIA", plot_annotations=True)

        See Also
        --------
        - self.plot_ice_area_dict : Internal dictionary for line labels and colors by ice type
        - pd.DataFrame.groupby : Used to compute daily climatologies
        """
        plt.figure(figsize=figsize, constrained_layout=True)
        cmap         = plt.get_cmap("tab10")
        month_labels = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
        xticks       = [1,32,62,92,122,152,182,212,242,272,302,332]
        self.plot_ice_area_dict = {"FI_BT": {"label": "daily ice-speed mask", "color": "black"},
                                   "FI_BT_bool": {"label": "Binary-days (6 of 7) ice-speed mask", "color": "orange"},
                                   "FI_BT_roll": {"label": "15-day-avg ice-speed mask", "color": "green"},}
        # Y-axis labels and defaults
        ice_type = ice_type.upper()
        ylabels = {"SIA": "Sea Ice Area (10⁶ km²)",
                   "SIV": "Sea Ice Volume (10⁶ km³)",
                   "FIA": "Fast Ice Area (10³ km²)",
                   "FIV": "Fast Ice Volume (km³)"}
        if ylim is None:
            ylim_defaults = {"SIA": (0, 20), "FIA": (100, 800)}
            ylim = ylim_defaults.get(ice_type, None)
        sim_annots = []
        sim_annot_colors = []
        for i, (name, da) in enumerate(area_dict.items()):
            if isinstance(da, xr.Dataset):
                da = da.to_array().squeeze()
            if name == "AF2020db_cli" and ice_type == "FIA":
                doy = da["doy"].values
                values = da.values
                plt.plot(doy, values, label="AF2020 Climatology (2000–2018)", linestyle="--", color="blue", linewidth=1.5)
                if plot_annotations:
                    max_doy, min_doy = doy[np.argmax(values)], doy[np.argmin(values)]
                    plt.plot(max_doy, values[np.argmax(values)], "b^")
                    plt.plot(min_doy, values[np.argmin(values)], "bv")
                    sim_annots += [f"AF2020 Max: {values[np.argmax(values)]:.1f} @ DOY {max_doy}", f"AF2020 Min: {values[np.argmin(values)]:.1f} @ DOY {min_doy}"]
                    sim_annot_colors += ["blue"] * 2
                continue
            if name == "NSIDC" and ice_type == "SIA":
                time = pd.to_datetime(da["time"].values)
                area = da.rolling(time=roll_days, center=True, min_periods=1).mean().values if roll_days >= 1 else da.values
                df = pd.DataFrame({"area": area}, index=time)
                df["doy"] = df.index.dayofyear
                clim_mean = df.groupby("doy")["area"].mean()
                plt.plot(clim_mean.index, clim_mean.values, label="NSIDC Climatology", linestyle="--", color="black", linewidth=1.5)
                if plot_annotations:
                    max_doy, min_doy = clim_mean.idxmax(), clim_mean.idxmin()
                    plt.plot(max_doy, clim_mean[max_doy], "k^")
                    plt.plot(min_doy, clim_mean[min_doy], "kv")
                    sim_annots += [f"NSIDC Max: {clim_mean[max_doy]:.1f} @ DOY {max_doy}", f"NSIDC Min: {clim_mean[min_doy]:.1f} @ DOY {min_doy}"]
                    sim_annot_colors += ["black"] * 2
                continue  
            time = pd.to_datetime(da["time"].values)
            area = da.rolling(time=roll_days, center=True, min_periods=1).mean().values if roll_days >= 1 else da.values
            df = pd.DataFrame({"area": area}, index=time)
            df["doy"] = df.index.dayofyear
            df["year"] = df.index.year
            df = df[df["year"] > df["year"].min()]  # drop first year
            grouped = df.groupby("doy")["area"]
            area_min, area_max, area_mean = grouped.min(), grouped.max(), grouped.mean()
            style = self.plot_ice_area_dict.get(name, {"label": name, "color": cmap(i)})
            plt.fill_between(area_min.index, area_min, area_max, alpha=0.2, color=style["color"], label=f"{style['label']} min/max range")
            plt.plot(area_mean.index, area_mean, color=style["color"], linewidth=2, label=f"{style['label']} climatology")
            if plot_annotations:
                max_doy, min_doy = area_mean.idxmax(), area_mean.idxmin()
                plt.plot(max_doy, area_mean[max_doy], marker="^", color=style["color"])
                plt.plot(min_doy, area_mean[min_doy], marker="v", color=style["color"])
                sim_annots += [f"{name} Max: {area_mean[max_doy]:.1f} @ DOY {max_doy}", f"{name} Min: {area_mean[min_doy]:.1f} @ DOY {min_doy}"]
                sim_annot_colors += [style["color"]] * 2
        plt.xticks(xticks, labels=month_labels, fontsize=tick_fontsize)
        plt.yticks(fontsize=tick_fontsize)
        plt.xlim(1, 366)    
        if ylim: plt.ylim(ylim)
        plt.ylabel(ylabels.get(ice_type, "Ice Metric"), fontsize=label_fontsize)
        plt.xlabel("Month", fontsize=label_fontsize)
        if tit_str:
            plt.title(tit_str, fontsize=title_fontsize)
        plt.grid(axis="y", linestyle="--", alpha=0.5)
        plt.legend(loc="upper left", fontsize=legend_fontsize)
        if sim_annots and plot_annotations:
            ax = plt.gca()
            for i, (text, color) in enumerate(zip(sim_annots, sim_annot_colors)):
                ax.text(0.78, 0.02 + i * 0.03, text, transform=ax.transAxes, fontsize=11, va='bottom', ha='center', color=color)
        if P_png and self.save_fig:
            plt.savefig(P_png, dpi=100)
            print(f"📏 Saved plot to {P_png}")

    def plot_timeseries(self, ts_dict,
                        primary_key : str = "FIA",
                        roll_days : int   = None,
                        tit_str   : str   = None,
                        ylim      : tuple = None,
                        ylabel    : str   = "Fast Ice Area (1000-km²)",
                        ytick_inc : int   = 100,
                        xlabel    : str   = "Date",
                        xtick_inc : int   = 1,
                        P_png     : str   = None,
                        fig_width : str   = "30c",
                        fig_height: str   = "5c",
                        legend_pos: str   = "JBL+jBL+w30c",
                        legend_box: str   = "+gwhite+p0.5p",
                        pen_weight: str   = "1p",
                        time_coord: str   = "time",
                        keys2plot : list  = None,
                        show_fig  : bool  = None):
        """
        Plot one or more time series using PyGMT.

        Parameters
        ----------
        ts_dict : dict
            Dictionary of {label: xarray.DataArray or Dataset}. Each value must have a 'time' coordinate.
        roll_days : int, optional
            Number of days for centered rolling average. If None or < 3, no smoothing is applied.
        tit_str : str, optional
            Title for the plot.
        ylim : tuple, optional
            Y-axis limits (ymin, ymax). Inferred from data if None.
        ylabel : str
            Y-axis label.
        ytick_inc : int
            Interval for y-axis ticks.
        xlabel : str
            X-axis label.
        xtick_inc : int
            Interval (in months or years) for x-axis ticks.
        P_png : str, optional
            Path to save PNG output. If None, figure is not saved.
        fig_width : str
            Width of the figure for PyGMT (e.g., "15c").
        fig_height : str
            Height of the figure for PyGMT (e.g., "5c").
        legend_pos : str
            PyGMT legend position string (e.g., "JTR+jTR+o0.2c").
        legend_box : str
            PyGMT legend box style string.
        pen_weight : str
            Pen thickness for the lines (e.g., "1p", "2p").
        time_coord : str
            Name of the time coordinate in each xarray object.
        keys2plot : list, optional
            Subset of keys in ts_dict to plot. If None, all keys are used.
        show_fig : bool, optional
            If True, show the figure interactively. Defaults to self.show_fig.

        Returns
        -------
        fig : pygmt.Figure
            The generated PyGMT figure.
        """
        show_fig = show_fig if show_fig is not None else self.show_fig
        dfs = []
        for i,(dict_key,data) in enumerate(ts_dict.items()):
            if dict_key=='AF2020':
                da = data
            else:
                da = data[primary_key]
            df         = pd.DataFrame({"time": pd.to_datetime(da["time"].values), 
                                        "data": da.values})
            dfs.append(df)
        df_all = pd.concat(dfs, ignore_index=True).dropna()
        all_times = df_all["time"]
        tmin, tmax = all_times.min(), all_times.max()
        # Auto y-limits with padding
        if ylim is None:
            ymin = df_all["data"].min()
            ymax = df_all["data"].max()
            pad = 0.05 * (ymax - ymin)
            ylim = [round(ymin - pad, 2), round(ymax + pad, 2)]
        region = [tmin.strftime("%Y-%m-%d"), tmax.strftime("%Y-%m-%d"), ylim[0], ylim[1]]
        projection = f"X{fig_width}/{fig_height}"
        frame = ["pxa3Of1o","sxa1Y",f"pya{ytick_inc}f50+l{ylabel}","WSnew"]
        # PyGMT settings
        fig = pygmt.Figure()
        fig.basemap(region=region, projection=projection, frame=frame)
        with pygmt.config(MAP_FRAME_TYPE          = "plain",
                     FONT_TITLE              = "18p,Helvetica-Bold",
                     FONT_LABEL              = "14p,Helvetica",
                     FONT_ANNOT_PRIMARY      = "20p,Helvetica",
                     MAP_TICK_LENGTH_PRIMARY = "0.1c",
                     FORMAT_DATE_MAP         = "o",
                     FORMAT_TIME_PRIMARY_MAP = "Character"):
            for i,(dict_key,data) in enumerate(ts_dict.items()):
                if dict_key=='AF2020':
                    da = data
                else:
                    da = data[primary_key]
                df         = pd.DataFrame({"time": pd.to_datetime(da["time"].values), 
                                           "data": da.values})
                color = self.plot_var_dict.get(dict_key, {}).get("line_clr", f"C{i}")
                if i==0:
                    leg_lab = f"{dict_key}+N{len(ts_dict.keys())}"
                else:
                    leg_lab = f"{dict_key}"
                fig.plot(x=df["time"], y=df["data"], pen=f"{pen_weight},{color}", label=leg_lab)
        fig.legend(position=legend_pos, box=legend_box)
        if P_png and self.save_fig:
            fig.savefig(P_png, dpi=600)
            print(f"📏 Saved plot to {P_png}")
        if show_fig:
            fig.show()
        pygmt.clib.Session.__exit__

    def plot_timeseries_groupby_doy(self, ts_dict,
                                    primary_key : str = "FIA",
                                    roll_days   : int   = None,
                                    tit_str     : str   = None,
                                    ylim        : tuple = None,
                                    ylabel      : str   = "Fast Ice Area (1000-km²)",
                                    ytick_inc   : int   = 100,
                                    xlabel      : str   = "Date",
                                    xtick_inc   : int   = 1,
                                    P_png       : str   = None,
                                    fig_width   : str   = "30c",
                                    fig_height  : str   = "10c",
                                    legend_pos  : str   = "JBR+jBR+o0.2c",
                                    pen_weight  : str   = "1p",
                                    time_coord  : str   = "time",
                                    keys2plot   : list  = None,
                                    show_fig    : bool  = None):
        """
        """
        dfs = []
        fake_year = 1996
        region = [f"{fake_year}-01-01", f"{fake_year}-12-31", 0, 1000]
        projection = f"X{fig_width}/{fig_height}"
        frame = ["WS", "pxa1O", f"sy+l{ylabel}"]#, "pya1+ucm"]  "sxa1Of30D",
        fig = pygmt.Figure()
        with pygmt.config(MAP_FRAME_TYPE          = "plain",
                          FONT_TITLE              = "18p,Helvetica-Bold",
                          FONT_LABEL              = "20p,Helvetica",
                          FONT_ANNOT_PRIMARY      = "20p,Helvetica",
                          MAP_TICK_LENGTH_PRIMARY = "0.1c",
                          FORMAT_DATE_MAP         = "o",
                          FORMAT_DATE_OUT         = "o",
                          FORMAT_TIME_PRIMARY_MAP = "a"):
            fig.basemap(projection=projection, region=region, frame=frame)
            for i,(dict_key,data) in enumerate(ts_dict.items()):
                if dict_key=='AF2020':
                    da = data
                else:
                    da = data[primary_key]
                df         = pd.DataFrame({"time": pd.to_datetime(da["time"].values), 
                                        "data": da.values})
                df["doy"]  = df["time"].dt.dayofyear
                df["year"] = df["time"].dt.year
                doy_grp    = df.groupby("doy")["data"]
                clim_min   = doy_grp.min()
                clim_max   = doy_grp.max()
                clim_mean  = doy_grp.mean()
                # Convert DOY index to datetime for use with pygmt
                clim_min.index  = pd.to_datetime(clim_min.index - 1, unit="D", origin=pd.Timestamp(f"{fake_year}-01-01"))
                clim_max.index  = pd.to_datetime(clim_max.index - 1, unit="D", origin=pd.Timestamp(f"{fake_year}-01-01"))
                clim_mean.index = pd.to_datetime(clim_mean.index - 1, unit="D", origin=pd.Timestamp(f"{fake_year}-01-01"))
                line_color      = self.plot_var_dict.get(dict_key,{}).get("line_clr",{})#def_colors[i % len(def_colors])
                fig.plot(x           = np.concatenate([clim_min.index, clim_max.index[::-1]]),
                        y            = np.concatenate([clim_min.values, clim_max.values[::-1]]),
                        fill         = f"{line_color}@80",
                        close        = True,
                        transparency = 60,)
                fig.plot(x     = clim_mean.index,
                        y     = clim_mean.values,
                        pen   = f"2p,{line_color}",
                        label = f"{dict_key}")
        fig.legend(position=legend_pos, box=True)
        fig.show()

    def plot_taylor(self, stats_dict, out_path):
        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_subplot(111, polar=True)
        ax.set_theta_zero_location('E')
        ax.set_theta_direction(-1)

        for rms in [0.2, 0.5, 1.0, 1.5]:
            rs = np.linspace(0.5, 2.0, 500)
            theta = np.arccos(np.clip(1 - (rms**2 - 1)/(2 * rs), -1, 1))
            ax.plot(theta, rs, '--', color='gray', lw=0.6)
        ax.plot([0], [1], 'ko', label='Reference')

        for label, stat in stats_dict.items():
            angle = np.arccos(stat["corr"])
            r = stat["std_ratio"]
            ax.plot(angle, r, 'o', label=label)

        ax.set_rmax(2)
        ax.set_rticks([0.5, 1.0, 1.5, 2.0])
        ax.set_rlabel_position(135)
        ax.set_title("Taylor Diagram: Sea Ice Speed Comparison", fontsize=12)
        ax.legend(loc='upper right', bbox_to_anchor=(1.45, 1.1))
        plt.tight_layout()
        plt.savefig(out_path)
        plt.close()