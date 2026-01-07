import os, time, imageio, pygmt
import xarray            as xr
import pandas            as pd
import geopandas         as gpd
import numpy             as np
import matplotlib.pyplot as plt
import matplotlib.dates  as mdates
from tqdm                import tqdm
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

    def prepare_data_for_pygmt_plot(self, da, bcoords=False, tcoords=True,
                                    lon_coord_name=None, lat_coord_name=None,
                                    diff_plot=False):
        """
        Prepare gridded data for PyGMT plotting.
        
        Priority:
        1. If lon_coord_name and lat_coord_name are provided -> use them
        2. Else, use bcoords/tcoords (cannot both be True)
        """
        data_dict = {}
        self.load_cice_grid(slice_hem=True)
        self.logger.info("preparing the data for plotting")
        # Determine which coordinates to use
        use_own_coords = lon_coord_name is not None or lat_coord_name is not None
        if use_own_coords:
            if lon_coord_name is None or lat_coord_name is None:
                raise ValueError("Both lon_coord_name and lat_coord_name must be provided if using own coordinates")
            self.logger.info(f"   using own coordinates: {lon_coord_name}, {lat_coord_name}")
            lon2d, lat2d = np.meshgrid(da[lon_coord_name], da[lat_coord_name])
        else:
            if bcoords and tcoords:
                raise ValueError("Cannot set both bcoords and tcoords to True")
            if bcoords:
                self.logger.info("   using B-grid coordinates")
                lon2d = self.G_u['lon'].values
                lat2d = self.G_u['lat'].values
            elif tcoords:
                self.logger.info("   using T-grid coordinates")
                lon2d = self.G_t['lon'].values
                lat2d = self.G_t['lat'].values
            else:
                raise ValueError("Must specify either bcoords, tcoords, or provide explicit coordinates")
        data2d = np.asarray(da.data).astype('float32')
        if diff_plot:
            mask = (data2d >= -1) & (data2d <= 1) & np.isfinite(data2d)
        else:
            mask = np.isfinite(data2d)
        data_dict['data'] = data2d[mask].ravel()
        data_dict['lon']  = lon2d[mask].ravel()
        data_dict['lat']  = lat2d[mask].ravel()
        return data_dict

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

    def extract_min_max_dates(self, ts_dict, keys2plot=None, primary_key='FIA', time_coord='time'):
        """
        Extract the minimum and maximum datetime values from time series data.

        This method scans through a dictionary of time series objects (either xarray.DataArrays,
        xarray.Datasets, or nested dictionaries containing a DataArray under `primary_key`),
        and returns the earliest and latest valid time values found across all selected entries.

        Parameters
        ----------
        ts_dict : dict
            A dictionary where each value is either:
            - an xarray.DataArray with a `time_coord` coordinate, or
            - a dictionary containing a DataArray under the key `primary_key`.
            Keys typically represent simulation or experiment identifiers.

        keys2plot : list of str, optional
            If provided, restricts the operation to keys within this list.
            Keys not in `keys2plot` will be skipped.

        primary_key : str, default 'FIA'
            The key to use when accessing nested dictionaries within `ts_dict`.
            Ignored if the value is already a DataArray or Dataset.

        time_coord : str, default 'time'
            The name of the time coordinate to extract from each DataArray.

        Returns
        -------
        tmin : pandas.Timestamp
            The earliest datetime found across all valid time series.

        tmax : pandas.Timestamp
            The latest datetime found across all valid time series.

        Notes
        -----
        - Entries in `ts_dict` with missing or malformed time coordinates are skipped.
        - The key "AF2020" is always excluded.
        - If no valid entries remain after filtering, a warning is logged and `None` is returned.

        """
        df_dts = []
        for dict_key, data in ts_dict.items():
            if (keys2plot is not None and dict_key not in keys2plot) or (dict_key == "AF2020"):
                continue
            self.logger.info(f"{dict_key} simulation will be included in {self._method_name()}()")
            # Handle both nested dict and direct DataArray
            if isinstance(data, dict) and primary_key in data:
                da = data[primary_key]
            else:
                da = data  # assume data is already a DataArray
            df_dt = pd.DataFrame({"time": pd.to_datetime(da[time_coord].values)})
            df_dts.append(df_dt)
        if not df_dts:
            self.logger.warning("No data to plot after filtering with keys2plot.")
            return
        df_all = pd.concat(df_dts, ignore_index=True).dropna()
        all_times = df_all["time"]
        return all_times.min(), all_times.max()

    def pygmt_fastice_panel(self,
                            fast_ice_variable : str   = "FIA",   # "FIA"|"fia" or "FIT"|"fit" or "FIS|fis"
                            ice_class         : str   = None,    # "FI_BT" ... see SeaIceToolbox for more help
                            class_type        : str   = "bin",   # 'bin', 'roll' or None
                            sim_name          : str   = None,
                            roll_days         : int   = 0,
                            # Generic (can be overridden)
                            fig_width         : str   = None,
                            fig_height        : str   = None,
                            ylim              : tuple = None,
                            frame_bndy        : str   = None,
                            yaxis_pri         : str   = None,
                            xaxis_pri         : str   = "a1Of15Dg",
                            leg_pos           : str   = None,
                            leg_box           : str   = "+gwhite+p.5p",
                            spat_var_style    : str   = None,
                            GI_plot_style     : str   = "c0.05c",
                            GI_fill_color     : str   = "#BA561A",
                            plot_GI           : bool  = False,
                            min_max_trans_val : int   = 80,
                            yshift_top        : str   = None,
                            yshift_bot        : str   = None,
                            bottom_frame_bndy : str   = "WSne",
                            bottom_yaxis      : str   = None,
                            bottom_xaxis      : str   = None,
                            land_clr          : str   = '#D1DDE0',
                            water_clr         : str   = "#EDF2F5",
                            coast_pen         : str   = "1/0.5p,black",
                            cbar_pos          : str   = None,
                            lon_coord_name    : str   = None,
                            lat_coord_name    : str   = None,
                            cmap              : str   = None,
                            series            : list  = None,
                            cbar_frame        : str   = None,
                            ANT_IS_pen        : str   = "0.2p,black",
                            ANT_IS_color      : str   = "#C1CED6",
                            font_annot_pri    : str   = "24p,Times-Roman",
                            font_annot_sec    : str   = "16p,Times-Roman",
                            font_lab          : str   = "22p,Times-Bold",
                            line_pen         : str    = "2p",
                            grid_pen_pri      : str   = ".5p",
                            grid_pen_sec      : str   = ".25p",
                            fmt_geo_map       : str   = "D:mm",
                            P_png             : str   = None,
                            save_fig          : bool  = None,
                            overwrite_fig     : bool  = None,
                            show_fig          : bool  = None):
        """
        Unified PyGMT fast-ice panel plotter.

        Set `fast_ice_variable` to:
            - "FIA" (or "fia")  -> plots FIA time series + FIP maps
            - "FIT" (or "fit")  -> plots FIT time series + FIHI maps

        Everything else (loading, styling, legends, grounded iceberg overlay, etc.)
        follows the same code path with variable-specific defaults injected via
        a small configuration dictionary.

        The rest of the arguments let you override those defaults if you need to.
        """
        var = fast_ice_variable.lower()
        if var not in ("fia", "fit", "fis", "fimar", "fimvr", "fitar", "fitvr"):
            raise ValueError(f"`fast_ice_variable` must be one of ['FIA','FIT','FIS','FIMAR','FIMVR','FITAR','FITVR']; got {fast_ice_variable}")
        # -------------------------------------------------------------------------
        # Per-variable defaults (you can push more things in here if you like)
        # -------------------------------------------------------------------------
        cfg = self.pygmt_FI_panel[var]
        # -------------------------------------------------------------------------
        # Resolve user overrides or fall back to defaults
        # -------------------------------------------------------------------------
        sim_name       = sim_name       if sim_name       is not None else self.sim_name
        ice_class      = ice_class      if ice_class      is not None else self.ice_class
        show_fig       = show_fig       if show_fig       is not None else self.show_fig
        save_fig       = save_fig       if save_fig       is not None else self.save_fig
        ow_fig         = overwrite_fig  if overwrite_fig  is not None else self.ow_fig
        frame_bndy     = frame_bndy     if frame_bndy     is not None else "WS"
        fig_width      = fig_width      if fig_width      is not None else "30c"
        fig_height     = fig_height     if fig_height     is not None else "25c"
        spat_var_style = spat_var_style if spat_var_style is not None else "s0.2c"
        yshift_top     = yshift_top     if yshift_top     is not None else "-6.25c"
        yshift_bot     = yshift_bot     if yshift_bot     is not None else "-10.5c"
        bottom_yaxis   = bottom_yaxis   if bottom_yaxis   is not None else "a5f1g"
        bottom_xaxis   = bottom_xaxis   if bottom_xaxis   is not None else "a30f10g"
        cbar_pos       = cbar_pos       if cbar_pos       is not None else "JBC+w25c/1c+mc+h"
        ylim           = ylim           if ylim           is not None else cfg["panel_ylim"]
        yaxis_pri      = yaxis_pri      if yaxis_pri      is not None else cfg["yaxis_pri"]
        leg_pos        = leg_pos        if leg_pos        is not None else cfg["leg_pos"]
        cbar_frame     = cbar_frame     if cbar_frame     is not None else cfg["cbar_frame"]
        cmap           = cmap           if cmap           is not None else cfg["cmap"]
        series         = series         if series         is not None else cfg["series"]
        # -------------------------------------------------------------------------
        # Paths / loads
        # -------------------------------------------------------------------------
        ANT_IS      = self.load_ice_shelves()
        ice_classes = [ice_class, f"{ice_class}_roll", f"{ice_class}_bin"]
        ts_dict     = {}
        if var=='fia':
            ts_dict["AF2020"] = xr.open_dataset(self.AF_FI_dict['P_AF2020_FIA'])["AF2020"]
        for iclass in ice_classes:
            P_mets          = Path(self.D_ispd_thresh, f"{iclass}_mets.zarr")
            ts_dict[iclass] = xr.open_dataset(P_mets)[cfg["top_name"]]
        tmin, tmax = self.extract_min_max_dates(ts_dict)
        # prioritise binary-days classification, then rolling-mean
        P_spat = Path(self.D_ispd_thresh, f"{ice_class}_{class_type}_mets.zarr") if class_type is not None else Path(self.D_ispd_thresh, f"{ice_class}_mets.zarr")
        try:
            da_spat = xr.open_dataset(P_spat)[cfg["bottom_name"]]
        except Exception:
            raise KeyError(f"could not load {P_spat}: {e}")
        df_spat = self.prepare_data_for_pygmt_plot(da_spat)
        if plot_GI:
            plot_GI_dict = self.load_GI_lon_lats()
        # -------------------------------------------------------------------------
        # Plot
        # -------------------------------------------------------------------------
        plot_region     = [f"{self.leap_year}-01-01", f"{self.leap_year}-12-31", ylim[0], ylim[1]]
        plot_projection = f"X{fig_width}/{fig_height}"
        frame           = [frame_bndy, f"px{xaxis_pri}", f"py{yaxis_pri}"]
        fig = pygmt.Figure()
        with pygmt.config(FONT_ANNOT_PRIMARY      = font_annot_pri,
                          FONT_ANNOT_SECONDARY    = font_annot_sec,
                          FONT_LABEL              = font_lab,
                          MAP_GRID_PEN_PRIMARY    = grid_pen_pri,
                          MAP_GRID_PEN_SECONDARY  = grid_pen_sec,
                          FORMAT_GEO_MAP          = fmt_geo_map,
                          FORMAT_DATE_MAP         = "o",
                          FORMAT_TIME_PRIMARY_MAP = "Abbreviated"):
            # ---- time-series top panel ----
            fig.basemap(region=plot_region, projection=plot_projection, **{"frame": frame})
            for k, da in ts_dict.items():
                if hasattr(self, "pygmt_FIA_dict") and k in self.pygmt_FIA_dict:
                    leg_lab    = self.pygmt_FIA_dict[k]["label"]
                    line_pen   = self.pygmt_FIA_dict[k]["line_pen"]
                    line_color = self.pygmt_FIA_dict[k]["line_color"]
                else:
                    leg_lab, line_pen, line_color = k, "1.5p", "black"
                clim = self.compute_doy_climatology(da)
                fig.plot(x            = np.concatenate([clim["min"].index, clim["max"].index[::-1]]),
                         y            = np.concatenate([clim["min"].values, clim["max"].values[::-1]]),
                         fill         = f"{line_color}@{min_max_trans_val}",
                         close        = True,
                         transparency = min_max_trans_val)
                fig.plot(x     = clim["mean"].index,
                         y     = clim["mean"].values,
                         pen   = f"{line_pen},{line_color}",
                         label = leg_lab)
            fig.legend(position=leg_pos, box=leg_box)
            # ---- bottom panel(s) ----
            fig.shift_origin(yshift=yshift_top)
            pygmt.makecpt(cmap=cmap, series=series)
            for i, (reg_name, reg_vals) in enumerate(self.Ant_2sectors.items()):
                b_region     = reg_vals["plot_region"]
                b_projection = reg_vals["projection"].format(fig_width=fig_width)
                if i > 0:
                    fig.shift_origin(yshift=yshift_bot)
                fig.basemap(region     = b_region,
                            projection = b_projection,
                            frame      = [f"x{bottom_xaxis}", f"y{bottom_yaxis}"])
                fig.coast(water=water_clr)
                fig.plot(x     = df_spat["lon"],
                         y     = df_spat["lat"],
                         fill  = df_spat["data"],
                         style = spat_var_style,
                         cmap  = True)
                if plot_GI:
                    fig.plot(x     = plot_GI_dict["lon"],
                             y     = plot_GI_dict["lat"],
                             fill  = GI_fill_color,
                             style = GI_plot_style)
                fig.coast(land=land_clr, shorelines=coast_pen)
                fig.plot(data=ANT_IS, pen=ANT_IS_pen, fill=ANT_IS_color)
            fig.colorbar(position=cbar_pos, frame=cbar_frame)
        # Save / Show
        if save_fig:
            F_png = f"{cfg['top_name']}_{cfg['bottom_name']}_{sim_name}_ispd_thresh{self.ispd_thresh_str}_{tmin.strftime('%Y')}-{tmax.strftime('%Y')}.png"
            P_png = P_png if P_png is not None else Path(self.D_graph, sim_name, F_png)
            P_png.parent.mkdir(parents=True, exist_ok=True)
            if not P_png.exists() or ow_fig:
                fig.savefig(P_png, dpi=300)
                self.logger.info(f"saved figure to {P_png}")
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
                               use_bcoords    = False,
                               use_tcoords    = False,
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
        #lon_coord_name = lon_coord_name if lon_coord_name is not None else self.pygmt_dict.get("lon_coord_name", "TLON")
        #lat_coord_name = lat_coord_name if lat_coord_name is not None else self.pygmt_dict.get("lat_coord_name", "TLAT")
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
        if lon_coord_name is not None or lat_coord_name is not None:
            plot_data_dict = self.prepare_data_for_pygmt_plot(da,
                                                              bcoords        = False,
                                                              tcoords        = False,
                                                              lon_coord_name = lon_coord_name,
                                                              lat_coord_name = lat_coord_name,
                                                              diff_plot      = diff_plot)
        else:
            if use_bcoords and use_tcoords:
                raise ValueError("Cannot set both use_bcoords and use_tcoords to True")
            plot_data_dict = self.prepare_data_for_pygmt_plot(da,
                                                              bcoords        = use_bcoords,
                                                              tcoords        = use_tcoords,
                                                              lon_coord_name = None,
                                                              lat_coord_name = None,
                                                              diff_plot      = diff_plot)
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
            plot_GI_dict = self.load_GI_lon_lats()#plot_data_dict)
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
            with pygmt.config(FONT_TITLE         = "16p,Courier-Bold",
                              FONT_ANNOT_PRIMARY = "14p,Helvetica",
                              COLOR_FOREGROUND   = 'black'):
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

    def pygmt_2D_data_prep(self, da,
                           x_coord_name = None,
                           y_coord_name = None,
                           region       = None,
                           extra_mask   = None,
                           mask_zero    = None):
        da2 = da.squeeze(drop=True)
        if da2.ndim != 2:
            raise ValueError(f"Expected 2-D data after squeeze; got dims {da2.dims}")
        x_name = x_coord_name or self.CICE_dict["lon_coord_name"]
        y_name = y_coord_name or self.CICE_dict["lat_coord_name"]
        xcoord = da2.coords[x_name]
        ycoord = da2.coords[y_name]
        if xcoord.ndim == 1 and ycoord.ndim == 1:
            X2D, Y2D = np.meshgrid(xcoord.values, ycoord.values, indexing="xy")
        elif xcoord.ndim == 2 and ycoord.ndim == 2:
            if xcoord.dims != da2.dims or ycoord.dims != da2.dims:
                da2 = da2.transpose(*xcoord.dims)
            X2D, Y2D = np.asarray(xcoord.values), np.asarray(ycoord.values)
        else:
            raise ValueError("Mixed coord dimensionality.")
        Z2D = np.asarray(da2.values)
        # --- decide zero-masking automatically if not specified ---
        if mask_zero is None:
            name = (da2.name or "").lower()
            fv = str(da2.attrs.get("flag_values", "")).strip()
            is_categorical = ("diff_cat" in name) or (fv in {"0 1 2", "0,1,2"})
            mask_zero = not is_categorical
        mask = np.isfinite(Z2D)
        if mask_zero:
            mask &= ~np.isclose(Z2D, 0.0, atol=1e-8)
        if region is not None:
            xmin, xmax, ymin, ymax = region
            mask &= (X2D >= xmin) & (X2D <= xmax) & (Y2D >= ymin) & (Y2D <= ymax)
        if extra_mask is not None:
            em = np.asarray(extra_mask(da2))
            if em.shape != Z2D.shape:
                raise ValueError("extra_mask shape mismatch")
            mask &= em
        x = X2D[mask]
        y = Y2D[mask]
        z = Z2D[mask]
        return np.column_stack([x, y, z])

    def pygmt_FIP_figure(self, plot_data,
                        show_fig      = False,
                        P_png         = None,
                        region        = [0, 360, -90, -62],
                        projection    = "S0.0/-90.0/50/20c",
                        cmap          = "/g/data/gv90/da1339/GRAPHICAL/CPTs/AF2020_YlGnBu.cpt",
                        reverse       = False,
                        series        = [0, 1],
                        cbar_pos      = "JBC+w15c/1c+mc+h",
                        cbar_frame    = ["xa0.1f0.05+lFast Ice Persistence"],
                        basemap_frame = ["af", f"+tAF2020-reG 2000-2018"],
                        land_color    = "#666666",
                        water_color   = "#BABCDE",
                        G_pt_color    = "#FF8903",
                        G_pt_marker   = "c",       # circle; for squares use "s"
                        G_pt_size     = "0.05",    # point size (string, no unit)
                        G_pt_unit     = "c",       # "c" for cm
                        shoreline_pen = "1/0.25p",
                        plot_GI         = False,
                        GI_color        = "#E349D0",
                        GI_marker       = "s",
                        GI_size         = "0.01",
                        GI_pt_unit      = "c",
                        plot_bathymetry = False,
                        # --- extras to help with DataArray inputs / categorical plotting ---
                        var_name      = None,        # e.g., "diff_cat" to trigger categorical mode
                        lat_coord_name = "lat",
                        lon_coord_name = "lon",
                        cat_cmap       = "hawaii",
                        cat_labels    = ("union", "model dominate", "observation alone"),
                        cat_series    = (0, 1, 2),   # codes that match cat_labels order
                        cat_cbar_frame = None, #["+lDifference category"],
                        # NEW:
                        weight_da=None,                 # xr.DataArray matching plot_data shape (optional)
                        weight_to_transparency=True):    # map weight∈[0,1] -> transparency∈[100,0]):
        """
        Plot FIP fields with PyGMT. If `var_name` contains 'diff_cat', uses a categorical CPT
        with codes 0,1,2 => agreement, simulation, observation.

        plot_data:
        - xr.DataArray with 2-D 'lat' and 'lon' coords, or
        - pd.DataFrame with columns ['longitude','latitude','z'].
        """
        # Prep GI overlay / bathy
        if plot_GI:
            self.load_cice_grid(slice_hem=True)
        if plot_bathymetry:
            SO_BATH = self.load_IBCSO_bath()
        # Convert input to a plotting table
        if isinstance(plot_data, xr.DataArray):
            da = plot_data
            # Use provided var_name if given, else fall back to DataArray.name
            if var_name is None and hasattr(da, "name") and da.name is not None:
                var_name = da.name
            lat = da[lat_coord_name].values.ravel()
            lon = da[lon_coord_name].values.ravel()
            val = da.values.ravel()
            valid = np.isfinite(val)
            df = pd.DataFrame( {"longitude": lon[valid], "latitude": lat[valid], "z": val[valid]})
        elif isinstance(plot_data, pd.DataFrame):
            df = plot_data.copy()
        else:
            raise TypeError("plot_data must be an xarray.DataArray or a pandas.DataFrame")
        # Build the figure + base
        fig = pygmt.Figure()
        with pygmt.config(FONT_TITLE           = "22p,Bookman-Demi",
                          FONT_ANNOT_PRIMARY   = "16p,NewCenturySchlbk-Roman",
                          FONT_ANNOT_SECONDARY = "16p,NewCenturySchlbk-Bold",
                          FONT_LABEL           = "16p,NewCenturySchlbk-Bold",
                          COLOR_FOREGROUND     = "black",):
            fig.basemap(region=region, projection=projection, frame=basemap_frame)
            if plot_bathymetry:
                fig.grdimage(grid=SO_BATH, cmap="geo")
            else:
                fig.coast(region=region, projection=projection, land=land_color, water=water_color)
            # Categorical vs continuous color mapping
            is_categorical = (var_name is not None) and ("diff_cat" in var_name.lower())
            if is_categorical:
                df["z"] = df["z"].astype(int)
                # categorical CPT + legend
                codes = np.unique(df["z"].to_numpy())
                zmin, zmax = int(codes.min()), int(codes.max())
                if zmax <= zmin:
                    zmax = zmin + 1
                pygmt.makecpt(cmap        = cat_cmap or "categorical",
                              series      = [zmin, zmax, 1],             # e.g., 0..2 step 1
                              color_model = "+c" + ",".join(cat_labels)) # e.g., "agreement,simulation,observation"
                if (weight_da is not None) and weight_to_transparency:
                    # Align weight with the same mask used to build df
                    if isinstance(plot_data, xr.DataArray):
                        z_full = plot_data.values.ravel()
                        valid  = np.isfinite(z_full)
                        w_full = weight_da.squeeze(drop=True).values.ravel()
                        wv     = w_full[valid]
                    else:
                        # If caller passed a DataFrame, require a 'weight' column
                        if "weight" not in df.columns:
                            raise ValueError("With DataFrame input, include a 'weight' column for opacity weighting.")
                        wv = df["weight"].to_numpy()
                    # Drop fully transparent points (weight<=0)
                    keep = wv > 0
                    if not np.all(keep):
                        df = df.loc[keep].reset_index(drop=True)
                        wv = wv[keep]
                    # Bin weights → a few opacity levels (scalar transparency per batch)
                    # Bins in weight-space; 4 bins is a good balance of speed and fidelity
                    w_edges = np.array([0.00, 0.25, 0.50, 0.75, 1.01])
                    bin_idx = np.digitize(wv, w_edges) - 1  # 0..3
                    df["_bin"] = bin_idx
                    # Representative transparency for each bin (tau = (1 - w_center)*100)
                    tau_for_bin = np.array([87.5, 62.5, 37.5, 12.5], dtype=float)
                    # Draw by category then by opacity bin (keeps CPT colors + avoids overlap artifacts)
                    for k in sorted(codes):
                        sel_k = (df["z"].to_numpy() == int(k))
                        if not np.any(sel_k):
                            continue
                        for b in range(len(tau_for_bin)):
                            sel = sel_k & (df["_bin"].to_numpy() == b)
                            if not np.any(sel):
                                continue
                            fig.plot(data         = df.loc[sel],
                                     style        = f"{G_pt_marker}{G_pt_size}{G_pt_unit}",
                                     cmap         = True,
                                     transparency = float(tau_for_bin[b]))   # <-- SCALAR, not array
                    # Clean up helper column
                    df.drop(columns=["_bin"], inplace=True)
                else:
                    # No weighting: single draw
                    fig.plot(data  = df,
                             style = f"{G_pt_marker}{G_pt_size}{G_pt_unit}",
                             cmap  = True)
                fig.colorbar(position=cbar_pos, frame=cat_cbar_frame)
            else:
                # Continuous CPT (default behavior)
                pygmt.makecpt(cmap=cmap, reverse=reverse, series=series)
                fig.plot(data  = df,
                         style = f"{G_pt_marker}{G_pt_size}{G_pt_unit}",
                         cmap  = True)
                fig.colorbar(position=cbar_pos, frame=cbar_frame)
            if plot_GI:
                fig.plot(x     = self.G_GI['lon'].values.ravel(), 
                         y     = self.G_GI['lat'].values.ravel(),
                         fill  = GI_color, 
                         style = f"{GI_marker}{GI_size}{GI_pt_unit}")
            fig.coast(region=region, projection=projection, shorelines=shoreline_pen)
        if P_png:
            fig.savefig(P_png)
        if show_fig:
            fig.show()
        return fig

    def _smooth_df_time(self, df: pd.DataFrame, window, center=True, min_periods=1):
        """
        Smooth a time series in df with columns ['time','data'].
        'window' can be an int (# of samples) or a time offset string like '15D'.
        Returns a new df with smoothed 'data'.
        """
        if window is None:
            return df
        s = df.set_index("time")["data"]
        # Pandas rolling supports both integer and time-based windows on a DatetimeIndex
        s_sm = s.rolling(window=window, center=center, min_periods=min_periods).mean()
        out = df.copy()
        out["data"] = s_sm.values
        return out

    def _repeat_doy_over_range(self, df_time_data, full_range, time_coord="time"):
        src         = df_time_data.copy()
        src["doy"]  = src[time_coord].dt.dayofyear
        clim_lookup = src.groupby("doy")["data"].mean()  # index = 1..365 or 366
        rep         = pd.DataFrame({"time": full_range})
        rep["doy"]  = rep["time"].dt.dayofyear
        # handle leap day: if 366 not in lookup, remap DOY=366 -> 365
        if 366 not in clim_lookup.index:
            rep["doy_eff"] = rep["doy"].where(rep["doy"] != 366, 365)
        else:
            rep["doy_eff"] = rep["doy"]
        rep["rep"] = rep["doy_eff"].map(clim_lookup)
        return rep[["time", "rep"]]

    def pygmt_timeseries(self, ts_dict,
                        comp_name        : str              = "test",
                        primary_key      : str              = "FIA",  # "FIA"
                        smooth           : str|int | None   = None,   # e.g., 15, "15D"
                        clim_smooth      : int | None       = None,   # 15
                        climatology      : bool             = False,
                        ylabel           : str              = None,   # "@[Fast Ice Area (1\\times10^3 km^2)@[",
                        ylim             : tuple            = [0,1000],
                        yaxis_pri        : int              = None,
                        ytick_pri        : int              = 100,
                        ytick_sec        : int              = 50,
                        projection       : str              = None,
                        fig_width        : str              = None,
                        fig_height       : str              = None,
                        xaxis_pri        : str              = None,
                        xaxis_sec        : str              = None,
                        frame_bndy       : str              = "WS",
                        legend_pos       : str              = None,
                        legend_box       : str              = "+gwhite+p0.5p",
                        fmt_dt_pri       : str              = None,
                        fmt_dt_sec       : str              = None,
                        fmt_dt_map       : str              = None,
                        fnt_type         : str              = "Helvetica",
                        fnt_wght_lab     : str              = "20p",
                        fnt_wght_ax      : str              = "18p",
                        line_pen         : str              = "1p",
                        grid_wght_pri    : str              = ".25p",
                        grid_wght_sec    : str              = ".1p",
                        P_png            : str              = None,
                        time_coord       : str              = "time",
                        time_coord_alt   : str              = "date",
                        keys2plot        : list             = None,
                        repeat_keys      : list[str] | None = None,
                        repeat_policy    : str              = "inside_others", # "inside_others" | "outside_others" | "fill_gaps" | "always"
                        repeat_ref_keys  : list[str] | None = None,            # which keys define the "others" window; default: all except current
                        clip_x_axis      : bool             = False,
                        zero_line        : bool             = False,
                        zero_line_level  : float            = 0.0,
                        zero_line_pen    : str              = "2p,black",
                        save_fig         : bool             = None,
                        show_fig         : bool             = None):
        """
        Plot time series of a primary variable (e.g., FIA) for a set of simulations or observations.

        Parameters
        ----------
        ts_dict : dict
            Dictionary of xarray DataArrays keyed by simulation or dataset name.
        comp_name : str
            Name for the comparison (used in figure title and filename).
        primary_key : str
            Key used to extract the variable from each dataset (except 'AF2020').
        climatology : bool
            If True, plot daily climatology with fill and mean lines.
        ylim : tuple or None
            Y-axis limits. If None, inferred from data with 5% padding.
        ylabel : str
            Label for the Y-axis.
        ytick_inc : int
            Interval between Y-axis ticks.
        xaxis_pri, xaxis_sec : str
            GMT frame settings for primary and secondary axes (when climatology is False).
        P_png : str or None
            Optional full path to save figure. If None, default filename is constructed.
        legend_box : str
            GMT legend box styling.
        line_pen : str
            Line thickness for plotted time series.
        time_coord : str
            Name of time coordinate in each DataArray.
        keys2plot : list or None
            If provided, only datasets with keys in this list are plotted.
        show_fig : bool or None
            If True, show figure interactively. Defaults to self.show_fig.
        repeat_keys : list[str] or None
            Keys in `ts_dict` to plot as a repeated day-of-year climatology when
            `climatology=False`.
        repeat_policy : {"outside_others","fill_gaps","always"}
            - "outside_others": keep original values where *other* series exist in time,
            but replace values outside that union with day-of-year climatology of
            the nominated series (good for fair comparison).
            - "fill_gaps": use climatology only where the nominated series has no data,
            keep original where it does.
            - "always": ignore original values and plot the repeated climatology across
            the full x-range.
        repeat_ref_keys : list[str] or None
            If provided, these keys define the "others" time span for the
            "outside_others" policy. By default, it's all plotted series except the
            current one.
        """
        show_fig   = show_fig   if show_fig   is not None else self.show_fig
        save_fig   = save_fig   if save_fig   is not None else self.save_fig
        fmt_dt_pri = fmt_dt_pri if fmt_dt_pri is not None else "Character"
        fmt_dt_sec = fmt_dt_sec if fmt_dt_sec is not None else "Abbreviated"
        fmt_dt_map = fmt_dt_map if fmt_dt_map is not None else "o"
        # primary_key = primary_key if primary_key is not None else keys2plot[0]
        # if keys2plot is None:
        #     raise("either 'primary_key' or 'keys2plot' must be defined")
        # need to get out the maximum times for plot boundaries
        tmin, tmax = self.extract_min_max_dates(ts_dict, keys2plot=keys2plot, primary_key=primary_key, time_coord=time_coord)
        # --- normalize new params ---
        if isinstance(repeat_keys, (str, bytes)):
            repeat_keys = [repeat_keys]
        repeat_keys = set(repeat_keys or [])
        if isinstance(repeat_ref_keys, (str, bytes)):
            repeat_ref_keys = [repeat_ref_keys]
        x_start, x_end = tmin, tmax
        if (not climatology) and clip_x_axis and repeat_keys:
            # reference keys = all plotted except the repeated ones, unless user specified
            if repeat_ref_keys:
                ref_keys = set(repeat_ref_keys)
            else:
                ref_keys = {k for k in ts_dict.keys()
                            if (keys2plot is None or k in (keys2plot or []))
                            and k not in repeat_keys}

            other_mins, other_maxs = [], []
            for k in ref_keys:
                da2 = ts_dict[k][primary_key]
                tt2 = pd.to_datetime(da2[time_coord].values)
                if len(tt2) > 0:
                    other_mins.append(tt2.min()); other_maxs.append(tt2.max())

            if other_mins:  # only override if we actually have refs
                x_start = pd.to_datetime(min(other_mins)).normalize()
                x_end   = pd.to_datetime(max(other_maxs)).normalize()
        # there are differences in the projection and x-axis for the two types of figures
        if climatology:
            fake_year  = 1996
            fig_width  = fig_width  if fig_width  is not None else "20c"
            fig_height = fig_height if fig_height is not None else "15c"
            xaxis_sec  = xaxis_sec  if xaxis_sec  is not None else None
            xaxis_pri  = xaxis_pri  if xaxis_pri  is not None else "a1Og"
            region     = [f"{fake_year}-01-01", f"{fake_year}-12-31", ylim[0], ylim[1]]
        else:
            fig_width  = fig_width  if fig_width  is not None else "50c"
            fig_height = fig_height if fig_height is not None else "15c"
            xaxis_pri  = xaxis_pri  if xaxis_pri  is not None else "a2Of30D"#g30D"
            xaxis_sec  = xaxis_sec  if xaxis_sec  is not None else "a1Y"
            region     = [x_start.strftime("%Y-%m-%d"), x_end.strftime("%Y-%m-%d"), ylim[0], ylim[1]]
        # define the projection and frame of the figure
        legend_pos = legend_pos if legend_pos is not None else f"JTL+jTL+o0.2c+w{fig_width}"
        projection = f"X{fig_width}/{fig_height}"
        yaxis_pri  = yaxis_pri if yaxis_pri is not None else f"a{ytick_pri}f{ytick_sec}+l{ylabel}"
        if xaxis_sec is not None:
            frame = [frame_bndy, f"sx{xaxis_sec}", f"px{xaxis_pri}", f"py{yaxis_pri}"]
        else:
            frame = [frame_bndy, f"px{xaxis_pri}", f"py{yaxis_pri}"]
        # make sure the figure has the same configuration by wrapping it in a with condition
        fig = pygmt.Figure()
        with pygmt.config(FONT_LABEL                = f"{fnt_wght_lab},{fnt_type}",
                          FONT                      = f"{fnt_wght_ax},{fnt_type}",
                          MAP_GRID_PEN_PRIMARY      = grid_wght_pri,
                          MAP_GRID_PEN_SECONDARY    = grid_wght_sec,
                          FORMAT_TIME_PRIMARY_MAP   = fmt_dt_pri,
                          FORMAT_TIME_SECONDARY_MAP = fmt_dt_sec,
                          FORMAT_DATE_MAP           = fmt_dt_map):
            fig.basemap(projection=projection, region=region, frame=frame)
            # loop over each key in the dictionary and if keys2plot is defined only plot those dictionaries
            cnt=0
            for i, (dict_key, data) in enumerate(ts_dict.items()):
                if data.get("line_clr", None) is not None:
                    line_color = data['line_clr']
                else:
                    line_color = self.plot_var_dict.get(dict_key, {}).get("line_clr", f"C{i}")
                if data.get("leg_lab", None) is not None:
                    leg_lab = data['leg_lab']
                else:
                    leg_lab = self.plot_var_dict.get(dict_key, {}).get("leg_lab", dict_key )
                da = data[primary_key]
                self.logger.info(f"pulling out data array for {dict_key} and putting into dataframe")
                self.logger.info(f"legend label: {leg_lab}")
                self.logger.info(f"line color  : {line_color}")
                if dict_key=="AF2020" and primary_key=="FIA":
                    df = pd.DataFrame({"time": pd.to_datetime(da[time_coord_alt].values), "data": da.values})
                else:
                    df = pd.DataFrame({"time": pd.to_datetime(da[time_coord].values), "data": da.values})
                if climatology:
                    if dict_key=="AF2020" and primary_key=="FIA":
                        clim = self.compute_doy_climatology(da, time_coord=time_coord_alt)
                    else:
                        clim = self.compute_doy_climatology(da)
                    mean_x = clim['mean'].index
                    mean_y = clim['mean'].values
                    if clim_smooth is not None and clim_smooth > 1:
                        mean_y = (pd.Series(mean_y, index=mean_x)
                                    .rolling(window=clim_smooth, center=True, min_periods=1)
                                    .mean()
                                    .values)
                    fig.plot(x            = np.concatenate([clim['min'].index, clim['max'].index[::-1]]),
                             y            = np.concatenate([clim['min'].values, clim['max'].values[::-1]]),
                             fill         = f"{line_color}@80",
                             close        = True,
                             transparency = 80)
                    fig.plot(x     = mean_x,
                             y     = mean_y,
                             pen   = f"{line_pen},{line_color}",
                             label = leg_lab)
                else:       
                    if dict_key in repeat_keys:
                        # target x-range (union for figure)
                        t_start, t_end = tmin.normalize(), tmax.normalize()
                        full_range = pd.date_range(t_start, t_end, freq="D")
                        # repeated climatology over the full range
                        rep_df = self._repeat_doy_over_range(df, full_range)  # cols: time, rep
                        # original reindexed over full range
                        base = pd.DataFrame({"time": full_range}).merge(df.rename(columns={"data": "orig"}), on="time", how="left")
                        mix = base.merge(rep_df, on="time", how="left")
                        if repeat_policy == "always":
                            out_df = mix[["time", "rep"]].rename(columns={"rep": "data"})
                        elif repeat_policy == "fill_gaps":
                            mix["data"] = mix["orig"].fillna(mix["rep"])
                            out_df = mix[["time", "data"]]
                        elif repeat_policy in ("outside_others", "inside_others"):
                            # Determine the "others" window
                            if repeat_ref_keys is not None:
                                ref_keys = set(repeat_ref_keys)
                            else:
                                ref_keys = set(k for k in ts_dict.keys()
                                            if k != dict_key and (keys2plot is None or k in keys2plot))
                            if len(ref_keys) == 0:
                                # no refs -> behave like fill_gaps (but still allow inside_others to clip nothing)
                                mix["data"] = mix["orig"].fillna(mix["rep"])
                                if repeat_policy == "inside_others":
                                    mix.loc[:, "data"] = np.nan  # nothing to show if no reference window
                                out_df = mix[["time", "data"]]
                            else:
                                # compute min/max over "others"
                                other_mins, other_maxs = [], []
                                for k in ref_keys:
                                    d2 = ts_dict[k][primary_key]
                                    tt = pd.to_datetime(d2[time_coord].values)
                                    if len(tt) > 0:
                                        other_mins.append(tt.min()); other_maxs.append(tt.max())
                                other_min = pd.to_datetime(min(other_mins)).normalize()
                                other_max = pd.to_datetime(max(other_maxs)).normalize()
                                inside = (mix["time"] >= other_min) & (mix["time"] <= other_max)
                                # inside the others’ window: use original where present, else repeated
                                mix.loc[inside, "data"] = mix.loc[inside, "orig"].fillna(mix.loc[inside, "rep"])
                                if repeat_policy == "inside_others":
                                    # outside -> hide (clip)
                                    mix.loc[~inside, "data"] = np.nan
                                else:
                                    # outside -> use repeated climatology
                                    mix.loc[~inside, "data"] = mix.loc[~inside, "rep"]
                                out_df = mix[["time", "data"]]
                        else:
                            raise ValueError(f"Unknown repeat_policy: {repeat_policy}")
                        # optional smoothing after composition
                        out_df = self._smooth_df_time(out_df, window=smooth)
                        fig.plot(x=out_df["time"], y=out_df["data"], pen=f"{line_pen},{line_color}", label=leg_lab)
                    else:
                        df_sm = self._smooth_df_time(df, window=smooth)
                        fig.plot(x=df_sm["time"], y=df_sm["data"], pen=f"{line_pen},{line_color}", label=leg_lab)
            if zero_line:
                y0 = float(zero_line_level)
                y_min, y_max = ylim
                if (y_min <= y0 <= y_max):
                    x0, x1 = region[0], region[1]
                    fig.plot(x=[x0, x1], y=[y0, y0], pen=zero_line_pen)
            fig.legend(position=legend_pos, box=legend_box)
        if save_fig:
            F_png = f"{primary_key}_{comp_name}_{'climatology_' if climatology else ''}{tmin.strftime('%Y')}-{tmax.strftime('%Y')}.png"
            P_png = P_png if P_png is not None else Path(self.D_graph, "timeseries", F_png)
            fig.savefig(P_png, dpi=300)
            self.logger.info(f"saved figure to {P_png}")
        if show_fig:
            fig.show()
        pygmt.clib.Session.__exit__

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
        
    def pygmt_map_plot_multi_var_8sectors(self, das, var_names,
                                          sim_name        = None,
                                            panel_titles    = None,
                                            hemisphere      = "south",   # kept for API symmetry; not used (8-sector only)
                                            time_stamp      = None,
                                            tit_str         = None,
                                            plot_GI         = False,
                                            diff_plot       = False,     # can be bool or list[bool] length N
                                            cmaps           = None,      # can be str or list[str] length N
                                            series_list     = None,      # can be list[3] or list[list[3]] length N
                                            reverse_list    = None,      # can be bool or list[bool] length N
                                            cbar_labels     = None,      # can be str or list[str] length N
                                            cbar_units_list = None,      # can be str or list[str] length N
                                            extend_cbar     = False,
                                            cbar_position   = None,
                                            lon_coord_name  = None,      # can be str or list[str] length N
                                            lat_coord_name  = None,      # can be str or list[str] length N
                                            use_bcoords     = False,     # can be bool or list[bool] length N
                                            use_tcoords     = False,     # can be bool or list[bool] length N
                                            fig_size        = None,      # width (cm) of EACH panel
                                            panel_gap_cm    = 0.6,       # horizontal gap between panels in cm
                                            var_sq_size     = 0.2,
                                            GI_sq_size      = 0.1,
                                            GI_fill_color   = "red",
                                            plot_iceshelves = True,
                                            plot_bathymetry = False,
                                            add_stat_annot  = False,
                                            land_color      = None,
                                            water_color     = None,
                                            P_png           = None,      # base filename; sector name will be appended
                                            var_out         = None,
                                            overwrite_fig   = None,
                                            show_fig        = None):
        """
        Multi-panel (2–3 columns) plot for the fixed 8 Antarctic sectors using PyGMT.

        - Produces ONE figure per sector.
        - Within each sector figure, places N panels left-to-right using fig.shift_origin().
        - N must be 2 or 3 (two or three DataArrays).

        Notes
        -----
        - If P_png is provided, the sector name is appended before the suffix to avoid overwrite,
        e.g. myplot.png -> myplot_WS.png.
        - For "diff" categorical plots, the method follows the existing convention:
        if "diff" in var_name.lower(): treats the field as categorical (agreement/simulation/observation).
        """
        # -------------------------
        # helpers
        # -------------------------
        def _as_list(x, n, name, default=None):
            if x is None:
                return [default] * n
            if isinstance(x, (list, tuple)):
                if len(x) != n:
                    raise ValueError(f"{name} must have length {n}; got {len(x)}")
                return list(x)
            return [x] * n
        # -------------------------
        # validate inputs
        # -------------------------
        if das is None or not isinstance(das, (list, tuple)):
            raise TypeError("das must be a list/tuple of 2 or 3 xarray.DataArray objects")
        n_pan = len(das)
        if n_pan not in (2, 3):
            raise ValueError(f"das must contain 2 or 3 DataArrays; got {n_pan}")
        if var_names is None or not isinstance(var_names, (list, tuple)) or len(var_names) != n_pan:
            raise ValueError(f"var_names must be a list/tuple of length {n_pan}")
        # defaults (mirrors your single-var method)
        sim_name    = sim_name if sim_name is not None else self.sim_name
        show_fig    = show_fig if show_fig is not None else self.show_fig
        ow_fig      = overwrite_fig if overwrite_fig is not None else self.ow_fig
        time_stamp  = time_stamp if time_stamp is not None else self.dt0_str
        fig_size    = fig_size if fig_size is not None else self.pygmt_dict["fig_size"]
        land_color  = land_color if land_color is not None else self.pygmt_dict["land_color"]
        water_color = water_color if water_color is not None else self.pygmt_dict["water_color"]
        cbar_pos    = cbar_position if cbar_position is not None else self.pygmt_dict["cbar_pos"].format(width=fig_size * 0.8, height=0.75)
        # per-panel lists
        panel_titles    = _as_list(panel_titles,    n_pan, "panel_titles",    default=None)
        diff_plot_list  = _as_list(diff_plot,       n_pan, "diff_plot",       default=False)
        cmaps           = _as_list(cmaps,           n_pan, "cmaps",           default=None)
        series_list     = _as_list(series_list,     n_pan, "series_list",     default=None)
        reverse_list    = _as_list(reverse_list,    n_pan, "reverse_list",    default=None)
        cbar_labels     = _as_list(cbar_labels,     n_pan, "cbar_labels",     default=None)
        cbar_units_list = _as_list(cbar_units_list, n_pan, "cbar_units_list", default=None)
        lon_coord_name_list = _as_list(lon_coord_name, n_pan, "lon_coord_name", default=None)
        lat_coord_name_list = _as_list(lat_coord_name, n_pan, "lat_coord_name", default=None)
        use_bcoords_list    = _as_list(use_bcoords,    n_pan, "use_bcoords",    default=False)
        use_tcoords_list    = _as_list(use_tcoords,    n_pan, "use_tcoords",    default=False)
        # resolve per-panel plot defaults from plot_var_dict where possible
        cmap_eff, series_eff, reverse_eff, cbar_lab_eff, cbar_unit_eff = [], [], [], [], []
        for j in range(n_pan):
            vn = var_names[j]
            if vn in getattr(self, "plot_var_dict", {}):
                cmap_eff.append(cmaps[j] if cmaps[j] is not None else self.plot_var_dict[vn]["cmap"])
                series_eff.append(series_list[j] if series_list[j] is not None else self.plot_var_dict[vn]["series"])
                reverse_eff.append(reverse_list[j] if reverse_list[j] is not None else self.plot_var_dict[vn]["reverse"])
                cbar_lab_eff.append(cbar_labels[j] if cbar_labels[j] is not None else self.plot_var_dict[vn]["name"])
                cbar_unit_eff.append(cbar_units_list[j] if cbar_units_list[j] is not None else self.plot_var_dict[vn]["units"])
            else:
                # require explicit overrides if vn not in plot_var_dict
                if cmaps[j] is None or series_list[j] is None:
                    raise KeyError(f"var_name='{vn}' not found in plot_var_dict. Provide cmaps[{j}] and series_list[{j}] explicitly.")
                cmap_eff.append(cmaps[j])
                series_eff.append(series_list[j])
                reverse_eff.append(False if reverse_list[j] is None else reverse_list[j])
                cbar_lab_eff.append(vn if cbar_labels[j] is None else cbar_labels[j])
                cbar_unit_eff.append("" if cbar_units_list[j] is None else cbar_units_list[j])
        if plot_iceshelves:
            ANT_IS = self.load_ice_shelves()
        if plot_bathymetry:
            SO_BATH = self.load_IBCSO_bath()
        if plot_GI:
            plot_GI_dict = self.load_GI_lon_lats()
        # output naming
        if var_out is None:
            var_out = "multi"
        if tit_str is None:
            tit_str = ""
        # x-shift between panels
        xshift = f"{(fig_size + panel_gap_cm)}c"
        # -------------------------
        # loop sectors (8-sector only)
        # -------------------------
        reg_dict = self.Ant_8sectors
        for reg_name, reg_vals in reg_dict.items():
            region     = reg_vals["plot_region"]
            projection = reg_vals["projection"]
            # sector projections usually need MC
            MC = self.get_meridian_center_from_geographic_extent(region)
            projection = projection.format(MC=MC, fig_size=fig_size)
            # sector output path handling
            P_png_sector = None
            if P_png is None and self.save_fig:
                P_png_sector = Path(self.D_graph, sim_name, reg_name, var_out, f"{time_stamp}_{sim_name}_{reg_name}_{var_out}.png")
            elif P_png is not None:
                P_png = Path(P_png)
                # append sector to avoid overwriting
                P_png_sector = P_png.with_name(f"{P_png.stem}_{reg_name}{P_png.suffix}")
            fig = pygmt.Figure()
            with pygmt.config(
                FONT_TITLE         = "16p,Courier-Bold",
                FONT_ANNOT_PRIMARY = "14p,Helvetica",
                COLOR_FOREGROUND   = "black"):
                # build each panel left-to-right using shift_origin()
                for j in range(n_pan):
                    if j > 0:
                        fig.shift_origin(xshift=xshift, yshift="0c")
                    da = das[j]
                    vn = var_names[j]
                    pn_title = panel_titles[j] if panel_titles[j] is not None else vn
                    # Prepare plotting coordinates/data
                    if lon_coord_name_list[j] is not None or lat_coord_name_list[j] is not None:
                        plot_data_dict = self.prepare_data_for_pygmt_plot(da,
                                                                          bcoords=False,
                                                                          tcoords=False,
                                                                          lon_coord_name=lon_coord_name_list[j],
                                                                          lat_coord_name=lat_coord_name_list[j],
                                                                          diff_plot=diff_plot_list[j])
                        lon_nm = lon_coord_name_list[j]
                        lat_nm = lat_coord_name_list[j]
                    else:
                        if use_bcoords_list[j] and use_tcoords_list[j]:
                            raise ValueError("Cannot set both use_bcoords and use_tcoords to True (panel index {j})")
                        plot_data_dict = self.prepare_data_for_pygmt_plot(da,
                                                                          bcoords=use_bcoords_list[j],
                                                                          tcoords=use_tcoords_list[j],
                                                                          lon_coord_name=None,
                                                                          lat_coord_name=None,
                                                                          diff_plot=diff_plot_list[j])
                        # these names are only needed for the categorical 'diff' branch below
                        lon_nm = self.pygmt_dict.get("lon_coord_name", "TLON")
                        lat_nm = self.pygmt_dict.get("lat_coord_name", "TLAT")
                    # Guard clauses (mirrors your one-var method)
                    if (not isinstance(plot_data_dict, dict)) or any(k not in plot_data_dict for k in ("lon","lat","data")):
                        self.logger.warning(f"[{reg_name}] panel {j}: plot_data_dict invalid — skipping panel.")
                        continue
                    if plot_data_dict["lon"] is None or plot_data_dict["lat"] is None or plot_data_dict["data"] is None:
                        self.logger.warning(f"[{reg_name}] panel {j}: plot_data_dict has None entries — skipping panel.")
                        continue
                    # Title logic
                    if tit_str and j == 0:
                        # Put overall title on first panel only
                        basemap_frame = ["af", f"+t{tit_str}"]
                    else:
                        basemap_frame = ["af", f"+t{pn_title}"]
                    fig.basemap(region=region, projection=projection, frame=basemap_frame)
                    # Background
                    if plot_bathymetry:
                        fig.grdimage(grid=SO_BATH, cmap="geo")
                    else:
                        fig.coast(region     = region,
                                  projection = projection,
                                  shorelines = "1/0.5p,gray30",
                                  land       = land_color,
                                  water      = water_color)
                    # Main variable plotting
                    if "diff" in vn.lower():
                        # categorical diff plot (expects integer labels)
                        lat   = da[lat_nm].values.flatten()
                        lon   = da[lon_nm].values.flatten()
                        val   = da.values.flatten()
                        valid = ~np.isnan(val)
                        df    = pd.DataFrame({"longitude": lon[valid],
                                              "latitude" : lat[valid],
                                              "z"        : val[valid].astype(int)})
                        pygmt.makecpt(cmap="categorical", series=[0, 2, 1],
                                    color_model="+cagreement,simulation,observation")
                        fig.plot(data=df, style=f"s{var_sq_size}c", cmap=True)
                    elif "mask" in vn.lower():
                        fig.plot(x=plot_data_dict["lon"], y=plot_data_dict["lat"],
                                fill="red", style=f"s{var_sq_size}c")
                    else:
                        pygmt.makecpt(cmap=cmap_eff[j], reverse=reverse_eff[j], series=series_eff[j])
                        fig.plot(x=plot_data_dict["lon"], y=plot_data_dict["lat"],
                                fill=plot_data_dict["data"], style=f"s{var_sq_size}c", cmap=True)
                    # Overlays
                    if plot_bathymetry:
                        fig.coast(region=region, projection=projection, shorelines="1/0.5p,gray30")
                    if plot_GI:
                        fig.plot(x=plot_GI_dict["lon"], y=plot_GI_dict["lat"],
                                fill=GI_fill_color, style=f"c{GI_sq_size}c")
                    if plot_iceshelves:
                        fig.plot(data=ANT_IS, pen="0.2p,gray", fill="lightgray")
                    # Optional stats annotation
                    if add_stat_annot:
                        try:
                            annot_text = self.generate_regional_annotation_stats(da, region, lon_nm, lat_nm)
                            for ii, line in enumerate(annot_text):
                                fig.text(position="TR", text=line, font="12p,Helvetica-Bold,black",
                                        justify="LM", no_clip=True, offset=f"-1/{-0.5*ii}")
                        except Exception as e:
                            self.logger.warning(f"[{reg_name}] panel {j}: stats annotation failed: {e}")
                    # Colorbar (relative to current panel origin, so it works with shift_origin)
                    if "diff" in vn.lower():
                        fig.colorbar(position=cbar_pos, frame=["x+l" + cbar_lab_eff[j]])
                    elif "mask" not in vn.lower():
                        cbar_frame = self.create_cbar_frame(
                            series_eff[j],
                            cbar_lab_eff[j],
                            units=cbar_unit_eff[j],
                            extend_cbar=extend_cbar,
                        )
                        try:
                            fig.colorbar(position=cbar_pos, frame=cbar_frame)
                        except pygmt.exceptions.GMTCLibError as e:
                            self.logger.warning(f"[{reg_name}] panel {j}: colorbar failed: {e}")

            # Save / show
            if P_png_sector is not None:
                if not P_png_sector.exists():
                    P_png_sector.parent.mkdir(parents=True, exist_ok=True)
                    fig.savefig(P_png_sector)
                    self.logger.info(f"Saved figure to {P_png_sector}")
                else:
                    if ow_fig:
                        fig.savefig(P_png_sector)
                        self.logger.info(f"Saved figure to {P_png_sector}")
                    else:
                        self.logger.info(f"{P_png_sector} exists and overwrite_fig=False; skipping save.")

            if show_fig:
                fig.show()

        # Important: do not manually call Session.__exit__ here; pygmt handles it.


# import os, time, imageio, pygmt
# import xarray            as xr
# import pandas            as pd
# import geopandas         as gpd
# import numpy             as np
# import matplotlib.pyplot as plt
# import matplotlib.dates  as mdates
# from tqdm                import tqdm
# from PIL                 import Image
# from pathlib             import Path
# from datetime            import datetime

# class SeaIcePlotter:

#     def __init__(self,**kwargs):
#         return

#     def load_ice_shelves(self):
#         """
#         Load and preprocess Antarctic ice shelf polygons for PyGMT plotting.

#         This method reads a shapefile containing Antarctic coastal geometries,
#         filters for polygons classified as ice shelves (`POLY_TYPE == 'S'`),
#         ensures valid geometry, reprojects them to WGS84 (EPSG:4326),
#         and applies a zero-width buffer to clean topology issues.

#         Returns
#         -------
#         geopandas.GeoSeries
#             Cleaned and reprojected geometries of Antarctic ice shelves.

#         Notes
#         -----
#         - The input shapefile path is read from `self.config['pygmt_dict']['P_coast_shape']`.
#         - This method is typically used to overlay ice shelf boundaries in PyGMT plots.
#         - The returned geometries can be passed directly to `pygmt.Figure.plot()`.

#         See Also
#         --------
#         - self.plot_FIA_FIP_faceted : Uses this method to overlay ice shelves in map panels.
#         - geopandas.read_file : For reading shapefiles.
#         """
#         gdf                     = gpd.read_file(self.config['pygmt_dict']['P_coast_shape'])
#         shelves                 = gdf[gdf['POLY_TYPE'] == 'S']
#         shelves                 = shelves[~shelves.geometry.is_empty & shelves.geometry.notnull()]
#         shelves                 = shelves.to_crs("EPSG:4326")
#         shelves.geometry        = shelves.geometry.buffer(0)
#         return shelves.geometry

#     def create_IBCSO_bath(self):
#         """
#         Extract and save a masked IBCSO bathymetry layer for Antarctic plotting.

#         This method loads the IBCSO v2.0 dataset (as a NetCDF raster), extracts
#         the seafloor depth (negative elevations), masks out land areas (positive or zero),
#         and saves a cleaned version to NetCDF for use in plotting.

#         The result is saved at the path specified in
#         `self.config['pygmt_dict']['P_IBCSO_bath']`.

#         Returns
#         -------
#         None

#         Notes
#         -----
#         - Input file is assumed to be a NetCDF raster with variable `band_data`.
#         - Only `band=0` is used (i.e., the main bathymetric band).
#         - Output variable is named `bath` and excludes attributes and encodings for simplicity.
#         - This method is typically called once to prepare the bathymetry layer before reuse.

#         See Also
#         --------
#         - self.load_IBCSO_bath : Loads the pre-saved masked bathymetry layer.
#         - https://www.ibcso.org/ : International Bathymetric Chart of the Southern Ocean.
#         """
#         ds               = xr.open_dataset(self.config['pygmt_dict']['P_IBCSO_bed'])
#         bed              = ds.band_data.isel(band=0)
#         bed_masked       = bed.where(bed < 0)
#         bed_masked.name  = "bath" 
#         bed_masked.attrs = {}    
#         ds_out           = bed_masked.to_dataset()
#         ds_out.attrs     = {}        
#         ds_out.encoding  = {}     
#         ds_out.to_netcdf(self.config['pygmt_dict']["P_IBCSO_bath"])

#     def load_IBCSO_bath(self):
#         """
#         Load masked IBCSO bathymetry dataset prepared by `create_IBCSO_bath()`.

#         Returns the `bath` variable from the NetCDF file specified in
#         `self.config['pygmt_dict']["P_IBCSO_bath"]`.

#         Returns
#         -------
#         xarray.DataArray
#             2D array of bathymetry values (only ocean depths, in meters).
#             Values are negative below sea level, NaN over land.

#         Notes
#         -----
#         - This method assumes that `create_IBCSO_bath()` has already been run.
#         - Designed for use in PyGMT background plotting or masking.

#         See Also
#         --------
#         - self.create_IBCSO_bath : Creates this file if it doesn't exist.
#         - xarray.open_dataset : Loads the bathymetry layer.
#         """
#         return xr.open_dataset(self.config['pygmt_dict']["P_IBCSO_bath"]).bath

#     def prepare_data_for_pygmt_plot(self, da, bcoords=False, tcoords=True,
#                                     lon_coord_name=None, lat_coord_name=None,
#                                     diff_plot=False):
#         """
#         Prepare gridded data for PyGMT plotting.
        
#         Priority:
#         1. If lon_coord_name and lat_coord_name are provided -> use them
#         2. Else, use bcoords/tcoords (cannot both be True)
#         """
#         data_dict = {}
#         self.load_bgrid(slice_hem=True)
#         self.logger.info("preparing the data for plotting")
#         # Determine which coordinates to use
#         use_own_coords = lon_coord_name is not None or lat_coord_name is not None
#         if use_own_coords:
#             if lon_coord_name is None or lat_coord_name is None:
#                 raise ValueError("Both lon_coord_name and lat_coord_name must be provided if using own coordinates")
#             self.logger.info(f"   using own coordinates: {lon_coord_name}, {lat_coord_name}")
#             lon2d, lat2d = np.meshgrid(da[lon_coord_name], da[lat_coord_name])
#         else:
#             if bcoords and tcoords:
#                 raise ValueError("Cannot set both bcoords and tcoords to True")
#             if bcoords:
#                 self.logger.info("   using B-grid coordinates")
#                 lon2d = self.G_u['lon'].values
#                 lat2d = self.G_u['lat'].values
#             elif tcoords:
#                 self.logger.info("   using T-grid coordinates")
#                 lon2d = self.G_t['lon'].values
#                 lat2d = self.G_t['lat'].values
#             else:
#                 raise ValueError("Must specify either bcoords, tcoords, or provide explicit coordinates")
#         data2d = np.asarray(da.data).astype('float32')
#         if diff_plot:
#             mask = (data2d >= -1) & (data2d <= 1) & np.isfinite(data2d)
#         else:
#             mask = np.isfinite(data2d)
#         data_dict['data'] = data2d[mask].ravel()
#         data_dict['lon']  = lon2d[mask].ravel()
#         data_dict['lat']  = lat2d[mask].ravel()
#         return data_dict

#     def create_cbar_frame(self, series, label, units=None, extend_cbar=False, max_ann_steps=10):
#         """
#         Construct a GMT-style colorbar annotation frame string for PyGMT plotting.

#         This utility generates a clean, readable colorbar frame string using adaptive
#         step sizing based on the data range. An optional second axis label (e.g., for units)
#         and extension arrows can be included.

#         Parameters
#         ----------
#         series : list or tuple of float
#             Data range for the colorbar as [vmin, vmax].
#         label : str
#             Label text for the colorbar (e.g., "Fast Ice Persistence").
#         units : str, optional
#             Units label to be shown along the secondary (y) axis (e.g., "1/100").
#         extend_cbar : bool, optional
#             Whether to append extension arrows to the colorbar (+e). Default is False.
#         max_ann_steps : int, default=10
#             Desired maximum number of major annotations (controls tick spacing).

#         Returns
#         -------
#         str or list of str
#             GMT-format colorbar annotation string (e.g., "a0.1f0.02+lLabel"), or
#             a list of two strings if `units` is provided.

#         Notes
#         -----
#         - Tick spacing is determined using a scaled logarithmic rounding to ensure clean steps (e.g., 0.1, 0.2, 0.5).
#         - If `extend_cbar` is True, `+e` is added to indicate out-of-bounds extension arrows.
#         - Compatible with `pygmt.Figure.colorbar(frame=...)`.

#         Examples
#         --------
#         >>> self.create_cbar_frame([0.01, 1.0], "Persistence", units="1/100")
#         ['a0.1f0.02+lPersistence', 'y+l 1/100']
#         """
#         vmin, vmax = series[0], series[1]
#         vrange = vmax - vmin
#         raw_step = vrange / max_ann_steps
#         exp = np.floor(np.log10(raw_step))
#         base = 10 ** exp
#         mult = raw_step / base
#         if mult < 1.5:
#             ann_step = base * 1
#         elif mult < 3:
#             ann_step = base * 2
#         elif mult < 7:
#             ann_step = base * 5
#         else:
#             ann_step = base * 10
#         tick_step = ann_step / 5
#         ann_str  = f"{ann_step:.3f}".rstrip("0").rstrip(".")
#         tick_str = f"{tick_step:.3f}".rstrip("0").rstrip(".")
#         # Build annotation string
#         frame = f"a{ann_str}f{tick_str}+l{label}"
#         if extend_cbar:
#             frame += "+e"  # or use "+eU" / "+eL" for one-sided arrows
#         if units is not None: 
#             return [frame, f"y+l {units}"]
#         else:
#             return frame

#     def get_meridian_center_from_geographic_extent(self, geographic_extent):
#         """
#         Determine the optimal central meridian for PyGMT polar stereographic projections.

#         Given a geographic extent in longitude/latitude format, this method calculates
#         the central meridian (longitude) to use in 'S<lon>/<lat>/<width>' PyGMT projection
#         strings. It accounts for dateline wrapping and ensures the plot is centered visually.

#         Parameters
#         ----------
#         geographic_extent : list of float
#             Geographic region as [min_lon, max_lon, min_lat, max_lat]. Accepts longitudes in
#             either [-180, 180] or [0, 360] and gracefully handles dateline crossing.

#         Returns
#         -------
#         float
#             Central meridian (longitude) in the range [-180, 180].

#         Notes
#         -----
#         - If the computed center falls outside the intended range (e.g., due to dateline wrapping), 
#             the method rotates the meridian 180° to better align the figure.
#         - The result is saved to `self.plot_meridian_center` for reuse.
#         - Used in PyGMT stereographic projections like `'S{lon}/-90/30c'`.

#         See Also
#         --------
#         - https://docs.generic-mapping-tools.org/latest/cookbook/proj.html
#         - pygmt.Figure.basemap : For setting map projections using meridian centers.
#         """
#         lon_min, lon_max = geographic_extent[0], geographic_extent[1]
#         lon_min_360 = lon_min % 360
#         lon_max_360 = lon_max % 360
#         if (lon_max_360 - lon_min_360) % 360 > 180:
#             center = ((lon_min_360 + lon_max_360 + 360) / 2) % 360
#         else:
#             center = (lon_min_360 + lon_max_360) / 2
#         if center > 180:
#             center -= 360
#         # Edge case fix: ensure center is visually aligned with geographic_extent
#         # If the computed center is 180° out of phase (i.e., upside-down plots)
#         if not (geographic_extent[0] <= center <= geographic_extent[1]):
#             # Flip 180°
#             center = (center + 180) % 360
#             if center > 180:
#                 center -= 360
#         #print(f"meridian center computed as {center:.2f}°")
#         self.plot_meridian_center = center
#         return center

#     def generate_regional_annotation_stats(self, da, region, lon_coord_name, lat_coord_name):
#         """
#         Generate summary statistics (mean, std, min, max) and their locations within a geographic region.

#         This helper function extracts the portion of a 2D `xarray.DataArray` that falls within a given 
#         geographic bounding box, computes spatial statistics on valid (non-NaN) values, and formats 
#         the results for annotation in PyGMT plots.

#         Parameters
#         ----------
#         da : xarray.DataArray
#             2D gridded data array (e.g., fast ice persistence) with latitude and longitude coordinates.
#         region : list or tuple of float
#             Bounding box for the region in the form [lon_min, lon_max, lat_min, lat_max].
#         lon_coord_name : str
#             Name of the longitude coordinate in `da` (e.g., 'TLON').
#         lat_coord_name : str
#             Name of the latitude coordinate in `da` (e.g., 'TLAT').

#         Returns
#         -------
#         list of str
#             Formatted text lines with:
#             - Mean
#             - Standard deviation
#             - Data extent (number of valid grid cells)
#             - Minimum value and its lat/lon location
#             - Maximum value and its lat/lon location

#         Notes
#         -----
#         - Longitude and latitude slices are derived using `np.searchsorted`, assuming monotonic grid coordinates.
#         - Only finite (non-NaN) values are included in statistics.
#         - Output is intended to be used directly in `fig.text()` or `fig.legend()` in PyGMT plotting routines.

#         Examples
#         --------
#         >>> stats = self.generate_regional_annotation_stats(FIP_DA, region=[0, 90, -75, -60], lon_coord_name='TLON', lat_coord_name='TLAT')
#         >>> for line in stats:
#         >>>     print(line)
#         Mean: 0.73
#         Std: 0.18
#         Extent: 1246
#         Min:  0.02 at (-68.41, 23.65)
#         Max: 1.00 at (-66.82, 45.12)
#         """
#         # Get the latitude and longitude coordinates (TLAT, TLON)
#         da_lats    = da[lat_coord_name].values
#         da_lons    = da[lon_coord_name].values
#         lon_min_ni = np.searchsorted(da_lons[0, :], region[0])
#         lon_max_ni = np.searchsorted(da_lons[0, :], region[1])        
#         lat_min_nj = np.searchsorted(da_lats[:, 0], region[2])
#         lat_max_nj = np.searchsorted(da_lats[:, 0], region[3])
#         da_sliced  = da.isel(nj=slice(lat_min_nj, lat_max_nj), ni=slice(lon_min_ni, lon_max_ni))
#         lon_sliced = da[lon_coord_name].isel(nj=slice(lat_min_nj, lat_max_nj), ni=slice(lon_min_ni, lon_max_ni))
#         lat_sliced = da[lat_coord_name].isel(nj=slice(lat_min_nj, lat_max_nj), ni=slice(lon_min_ni, lon_max_ni))
#         # Get the data, latitudes, and longitudes
#         data       = da_sliced.values.flatten()
#         lat        = lat_sliced.values.flatten()
#         lon        = lon_sliced.values.flatten()
#         valid_mask = ~np.isnan(data)
#         data_valid = data[valid_mask]
#         lat_valid  = lat[valid_mask]
#         lon_valid  = lon[valid_mask]
#         # Calculate statistics
#         spatial_mean   = np.mean(data_valid)
#         spatial_std    = np.std(data_valid)
#         spatial_extent = len(data_valid)
#         data_min       = np.min(data_valid)
#         data_max       = np.max(data_valid)
#         # Get min/max locations (lat, lon)
#         min_index    = np.argmin(data_valid)
#         max_index    = np.argmax(data_valid)
#         min_location = (lat_valid[min_index], lon_valid[min_index])
#         max_location = (lat_valid[max_index], lon_valid[max_index])
#         # Format text to display
#         text  = [f"Mean: {spatial_mean:.2f}", f"Std: {spatial_std:.2f}", f"Extent: {spatial_extent}"]
#         text += [f"Min:  {data_min:.2f} at {min_location}", f"Max: {data_max:.2f} at {max_location}"]
#         return text

#     def extract_min_max_dates(self, ts_dict, keys2plot=None, primary_key='FIA', time_coord='time'):
#         """
#         Extract the minimum and maximum datetime values from time series data.

#         This method scans through a dictionary of time series objects (either xarray.DataArrays,
#         xarray.Datasets, or nested dictionaries containing a DataArray under `primary_key`),
#         and returns the earliest and latest valid time values found across all selected entries.

#         Parameters
#         ----------
#         ts_dict : dict
#             A dictionary where each value is either:
#             - an xarray.DataArray with a `time_coord` coordinate, or
#             - a dictionary containing a DataArray under the key `primary_key`.
#             Keys typically represent simulation or experiment identifiers.

#         keys2plot : list of str, optional
#             If provided, restricts the operation to keys within this list.
#             Keys not in `keys2plot` will be skipped.

#         primary_key : str, default 'FIA'
#             The key to use when accessing nested dictionaries within `ts_dict`.
#             Ignored if the value is already a DataArray or Dataset.

#         time_coord : str, default 'time'
#             The name of the time coordinate to extract from each DataArray.

#         Returns
#         -------
#         tmin : pandas.Timestamp
#             The earliest datetime found across all valid time series.

#         tmax : pandas.Timestamp
#             The latest datetime found across all valid time series.

#         Notes
#         -----
#         - Entries in `ts_dict` with missing or malformed time coordinates are skipped.
#         - The key "AF2020" is always excluded.
#         - If no valid entries remain after filtering, a warning is logged and `None` is returned.

#         """
#         df_dts = []
#         for dict_key, data in ts_dict.items():
#             if (keys2plot is not None and dict_key not in keys2plot) or (dict_key == "AF2020"):
#                 continue
#             self.logger.info(f"{dict_key} simulation will be included in {self._method_name()}()")
#             # Handle both nested dict and direct DataArray
#             if isinstance(data, dict) and primary_key in data:
#                 da = data[primary_key]
#             else:
#                 da = data  # assume data is already a DataArray
#             df_dt = pd.DataFrame({"time": pd.to_datetime(da[time_coord].values)})
#             df_dts.append(df_dt)
#         if not df_dts:
#             self.logger.warning("No data to plot after filtering with keys2plot.")
#             return
#         df_all = pd.concat(df_dts, ignore_index=True).dropna()
#         all_times = df_all["time"]
#         return all_times.min(), all_times.max()

#     def pygmt_fastice_panel(self,
#                             fast_ice_variable : str   = "FIA",   # "FIA"|"fia" or "FIT"|"fit" or "FIS|fis"
#                             ice_class         : str   = None,    # "FI_BT" ... see SeaIceToolbox for more help
#                             class_type        : str   = "bin",   # 'bin', 'roll' or None
#                             sim_name          : str   = None,
#                             roll_days         : int   = 0,
#                             # Generic (can be overridden)
#                             fig_width         : str   = None,
#                             fig_height        : str   = None,
#                             ylim              : tuple = None,
#                             frame_bndy        : str   = None,
#                             yaxis_pri         : str   = None,
#                             xaxis_pri         : str   = "a1Of15Dg",
#                             leg_pos           : str   = None,
#                             leg_box           : str   = "+gwhite+p.5p",
#                             spat_var_style    : str   = None,
#                             GI_plot_style     : str   = "c0.05c",
#                             GI_fill_color     : str   = "#BA561A",
#                             plot_GI           : bool  = False,
#                             min_max_trans_val : int   = 80,
#                             yshift_top        : str   = None,
#                             yshift_bot        : str   = None,
#                             bottom_frame_bndy : str   = "WSne",
#                             bottom_yaxis      : str   = None,
#                             bottom_xaxis      : str   = None,
#                             land_clr          : str   = '#D1DDE0',
#                             water_clr         : str   = "#EDF2F5",
#                             coast_pen         : str   = "1/0.5p,black",
#                             cbar_pos          : str   = None,
#                             lon_coord_name    : str   = None,
#                             lat_coord_name    : str   = None,
#                             cmap              : str   = None,
#                             series            : list  = None,
#                             cbar_frame        : str   = None,
#                             ANT_IS_pen        : str   = "0.2p,black",
#                             ANT_IS_color      : str   = "#C1CED6",
#                             font_annot_pri    : str   = "24p,Times-Roman",
#                             font_annot_sec    : str   = "16p,Times-Roman",
#                             font_lab          : str   = "22p,Times-Bold",
#                             line_pen         : str    = "2p",
#                             grid_pen_pri      : str   = ".5p",
#                             grid_pen_sec      : str   = ".25p",
#                             fmt_geo_map       : str   = "D:mm",
#                             P_png             : str   = None,
#                             save_fig          : bool  = None,
#                             overwrite_fig     : bool  = None,
#                             show_fig          : bool  = None):
#         """
#         Unified PyGMT fast-ice panel plotter.

#         Set `fast_ice_variable` to:
#             - "FIA" (or "fia")  -> plots FIA time series + FIP maps
#             - "FIT" (or "fit")  -> plots FIT time series + FIHI maps

#         Everything else (loading, styling, legends, grounded iceberg overlay, etc.)
#         follows the same code path with variable-specific defaults injected via
#         a small configuration dictionary.

#         The rest of the arguments let you override those defaults if you need to.
#         """
#         var = fast_ice_variable.lower()
#         if var not in ("fia", "fit", "fis", "fimar", "fimvr", "fitar", "fitvr"):
#             raise ValueError(f"`fast_ice_variable` must be one of ['FIA','FIT','FIS','FIMAR','FIMVR','FITAR','FITVR']; got {fast_ice_variable}")
#         # -------------------------------------------------------------------------
#         # Per-variable defaults (you can push more things in here if you like)
#         # -------------------------------------------------------------------------
#         cfg = self.pygmt_FI_panel[var]
#         # -------------------------------------------------------------------------
#         # Resolve user overrides or fall back to defaults
#         # -------------------------------------------------------------------------
#         sim_name       = sim_name       if sim_name       is not None else self.sim_name
#         ice_class      = ice_class      if ice_class      is not None else self.ice_class
#         show_fig       = show_fig       if show_fig       is not None else self.show_fig
#         save_fig       = save_fig       if save_fig       is not None else self.save_fig
#         ow_fig         = overwrite_fig  if overwrite_fig  is not None else self.ow_fig
#         frame_bndy     = frame_bndy     if frame_bndy     is not None else "WS"
#         fig_width      = fig_width      if fig_width      is not None else "30c"
#         fig_height     = fig_height     if fig_height     is not None else "25c"
#         spat_var_style = spat_var_style if spat_var_style is not None else "s0.2c"
#         yshift_top     = yshift_top     if yshift_top     is not None else "-6.25c"
#         yshift_bot     = yshift_bot     if yshift_bot     is not None else "-10.5c"
#         bottom_yaxis   = bottom_yaxis   if bottom_yaxis   is not None else "a5f1g"
#         bottom_xaxis   = bottom_xaxis   if bottom_xaxis   is not None else "a30f10g"
#         cbar_pos       = cbar_pos       if cbar_pos       is not None else "JBC+w25c/1c+mc+h"
#         ylim           = ylim           if ylim           is not None else cfg["panel_ylim"]
#         yaxis_pri      = yaxis_pri      if yaxis_pri      is not None else cfg["yaxis_pri"]
#         leg_pos        = leg_pos        if leg_pos        is not None else cfg["leg_pos"]
#         cbar_frame     = cbar_frame     if cbar_frame     is not None else cfg["cbar_frame"]
#         cmap           = cmap           if cmap           is not None else cfg["cmap"]
#         series         = series         if series         is not None else cfg["series"]
#         # -------------------------------------------------------------------------
#         # Paths / loads
#         # -------------------------------------------------------------------------
#         ANT_IS      = self.load_ice_shelves()
#         ice_classes = [ice_class, f"{ice_class}_roll", f"{ice_class}_bin"]
#         ts_dict     = {}
#         if var=='fia':
#             ts_dict["AF2020"] = xr.open_dataset(self.AF_FI_dict['P_AF2020_FIA'])["AF2020"]
#         for iclass in ice_classes:
#             P_mets          = Path(self.D_ispd_thresh, f"{iclass}_mets.zarr")
#             ts_dict[iclass] = xr.open_dataset(P_mets)[cfg["top_name"]]
#         tmin, tmax = self.extract_min_max_dates(ts_dict)
#         # prioritise binary-days classification, then rolling-mean
#         P_spat = Path(self.D_ispd_thresh, f"{ice_class}_{class_type}_mets.zarr") if class_type is not None else Path(self.D_ispd_thresh, f"{ice_class}_mets.zarr")
#         try:
#             da_spat = xr.open_dataset(P_spat)[cfg["bottom_name"]]
#         except Exception:
#             raise KeyError(f"could not load {P_spat}: {e}")
#         df_spat = self.prepare_data_for_pygmt_plot(da_spat)
#         if plot_GI:
#             plot_GI_dict = self.load_GI_lon_lats()
#         # -------------------------------------------------------------------------
#         # Plot
#         # -------------------------------------------------------------------------
#         plot_region     = [f"{self.leap_year}-01-01", f"{self.leap_year}-12-31", ylim[0], ylim[1]]
#         plot_projection = f"X{fig_width}/{fig_height}"
#         frame           = [frame_bndy, f"px{xaxis_pri}", f"py{yaxis_pri}"]
#         fig = pygmt.Figure()
#         with pygmt.config(FONT_ANNOT_PRIMARY      = font_annot_pri,
#                           FONT_ANNOT_SECONDARY    = font_annot_sec,
#                           FONT_LABEL              = font_lab,
#                           MAP_GRID_PEN_PRIMARY    = grid_pen_pri,
#                           MAP_GRID_PEN_SECONDARY  = grid_pen_sec,
#                           FORMAT_GEO_MAP          = fmt_geo_map,
#                           FORMAT_DATE_MAP         = "o",
#                           FORMAT_TIME_PRIMARY_MAP = "Abbreviated"):
#             # ---- time-series top panel ----
#             fig.basemap(region=plot_region, projection=plot_projection, **{"frame": frame})
#             for k, da in ts_dict.items():
#                 if hasattr(self, "pygmt_FIA_dict") and k in self.pygmt_FIA_dict:
#                     leg_lab    = self.pygmt_FIA_dict[k]["label"]
#                     line_pen   = self.pygmt_FIA_dict[k]["line_pen"]
#                     line_color = self.pygmt_FIA_dict[k]["line_color"]
#                 else:
#                     leg_lab, line_pen, line_color = k, "1.5p", "black"
#                 clim = self.compute_doy_climatology(da)
#                 fig.plot(x            = np.concatenate([clim["min"].index, clim["max"].index[::-1]]),
#                          y            = np.concatenate([clim["min"].values, clim["max"].values[::-1]]),
#                          fill         = f"{line_color}@{min_max_trans_val}",
#                          close        = True,
#                          transparency = min_max_trans_val)
#                 fig.plot(x     = clim["mean"].index,
#                          y     = clim["mean"].values,
#                          pen   = f"{line_pen},{line_color}",
#                          label = leg_lab)
#             fig.legend(position=leg_pos, box=leg_box)
#             # ---- bottom panel(s) ----
#             fig.shift_origin(yshift=yshift_top)
#             pygmt.makecpt(cmap=cmap, series=series)
#             for i, (reg_name, reg_vals) in enumerate(self.Ant_2sectors.items()):
#                 b_region     = reg_vals["plot_region"]
#                 b_projection = reg_vals["projection"].format(fig_width=fig_width)
#                 if i > 0:
#                     fig.shift_origin(yshift=yshift_bot)
#                 fig.basemap(region     = b_region,
#                             projection = b_projection,
#                             frame      = [f"x{bottom_xaxis}", f"y{bottom_yaxis}"])
#                 fig.coast(water=water_clr)
#                 fig.plot(x     = df_spat["lon"],
#                          y     = df_spat["lat"],
#                          fill  = df_spat["data"],
#                          style = spat_var_style,
#                          cmap  = True)
#                 if plot_GI:
#                     fig.plot(x     = plot_GI_dict["lon"],
#                              y     = plot_GI_dict["lat"],
#                              fill  = GI_fill_color,
#                              style = GI_plot_style)
#                 fig.coast(land=land_clr, shorelines=coast_pen)
#                 fig.plot(data=ANT_IS, pen=ANT_IS_pen, fill=ANT_IS_color)
#             fig.colorbar(position=cbar_pos, frame=cbar_frame)
#         # Save / Show
#         if save_fig:
#             F_png = f"{cfg['top_name']}_{cfg['bottom_name']}_{sim_name}_ispd_thresh{self.ispd_thresh_str}_{tmin.strftime('%Y')}-{tmax.strftime('%Y')}.png"
#             P_png = P_png if P_png is not None else Path(self.D_graph, sim_name, F_png)
#             P_png.parent.mkdir(parents=True, exist_ok=True)
#             if not P_png.exists() or ow_fig:
#                 fig.savefig(P_png, dpi=300)
#                 self.logger.info(f"saved figure to {P_png}")
#             else:
#                 self.logger.info(f"{P_png} already exists and not overwriting")
#         if show_fig:
#             fig.show()
#         pygmt.clib.Session.__exit__

#     def pygmt_map_plot_one_var(self, da, var_name,
#                                sim_name       = None,
#                                plot_regions   = None,
#                                regional_dict  = None,
#                                hemisphere     = "south",
#                                time_stamp     = None,
#                                tit_str        = None,
#                                plot_GI        = False,
#                                diff_plot      = False,
#                                cmap           = None,
#                                series         = None,
#                                reverse        = None,
#                                cbar_label     = None,
#                                cbar_units     = None,
#                                extend_cbar    = False,
#                                cbar_position  = None,
#                                lon_coord_name = None,
#                                lat_coord_name = None,
#                                use_bcoords    = False,
#                                use_tcoords    = False,
#                                fig_size       = None,
#                                var_sq_size    = 0.2,
#                                GI_sq_size     = 0.1,
#                                GI_fill_color  = "red",
#                                plot_iceshelves= True,
#                                plot_bathymetry= True,
#                                add_stat_annot = False,
#                                land_color     = None,
#                                water_color    = None,
#                                P_png          = None,
#                                var_out        = None,
#                                overwrite_fig  = None,
#                                show_fig       = None):
#         """
#         Generate a PyGMT figure showing a variable (e.g., FIA, FIP, differences) as a spatial map
#         over Antarctic regions or hemispheric view, optionally including grounded icebergs,
#         ice shelf outlines, and bathymetry.

#         This flexible mapping function supports multiple types of visualizations:
#         - Scalar field plots (e.g., fast ice persistence, concentration)
#         - Binary or categorical masks (e.g., agreement maps, simulation masks)
#         - Difference plots (e.g., observation minus simulation)

#         Parameters
#         ----------
#         da : xarray.DataArray
#             Input 2D (or broadcastable) data array to be plotted, e.g., fast ice persistence.
#         var_name : str
#             Name of the variable to be plotted. Used to look up color map settings and labels.
#         sim_name : str, optional
#             Simulation name used for file naming and figure annotation. Defaults to `self.sim_name`.
#         plot_regions : int or None, optional
#             Number of regional plots:
#             - 8 : Antarctic 8-sector view
#             - 2 : East/West sectors
#             - None : Full hemisphere plot (default)
#         regional_dict : dict, optional
#             Custom dictionary of plotting regions (overrides built-ins if provided).
#         hemisphere : str, default="south"
#             Hemisphere name for hemispheric plot. Typically "south".
#         time_stamp : str, optional
#             Timestamp string for file naming. Defaults to `self.dt0_str`.
#         tit_str : str, optional
#             Title string to display on the figure.
#         plot_GI : bool, default=False
#             Whether to overlay grounded iceberg locations using `load_GI_lon_lats()`.
#         diff_plot : bool, default=False
#             Whether this is a difference plot (used for labeling, masking, and symbology).
#         cmap : str, optional
#             Color map name for continuous scalar plots.
#         series : list of float, optional
#             Min, max, and increment for the colorbar (e.g., [0, 1, 0.1]).
#         reverse : bool, optional
#             Whether to reverse the colormap.
#         cbar_label : str, optional
#             Label to show on the colorbar.
#         cbar_units : str, optional
#             Units for the colorbar (shown on secondary axis).
#         extend_cbar : bool, default=False
#             Whether to add extension arrows to colorbar.
#         cbar_position : str, optional
#             PyGMT-compatible position string for placing the colorbar.
#         lon_coord_name : str, optional
#             Longitude coordinate name in `da`. Defaults to `self.pygmt_dict`.
#         lat_coord_name : str, optional
#             Latitude coordinate name in `da`. Defaults to `self.pygmt_dict`.
#         fig_size : float, optional
#             Figure width in centimeters for hemispheric plot.
#         var_sq_size : float, default=0.2
#             Marker size (in cm) for main plotted variable (for scatter plotting).
#         GI_sq_size : float, default=0.1
#             Marker size (in cm) for grounded iceberg overlay.
#         GI_fill_color : str, default="red"
#             Color to use for grounded iceberg markers.
#         plot_iceshelves : bool, default=True
#             Whether to overlay Antarctic ice shelf outlines.
#         plot_bathymetry : bool, default=True
#             Whether to plot IBCSO bathymetry using shaded relief (`grdimage`).
#         add_stat_annot : bool, default=False
#             Whether to annotate figure with basic regional statistics.
#         land_color : str, optional
#             Color for land. Defaults to `self.pygmt_dict['land_color']`.
#         water_color : str, optional
#             Color for ocean/water. Defaults to `self.pygmt_dict['water_color']`.
#         P_png : pathlib.Path, optional
#             File path to save the figure. If None, path is generated automatically.
#         var_out : str, optional
#             Output variable name used in file naming. Defaults to `var_name`.
#         overwrite_fig : bool, optional
#             Whether to overwrite existing figure if it exists.
#         show_fig : bool, optional
#             Whether to display the figure interactively.

#         Returns
#         -------
#         None
#             The method generates and optionally saves or displays a PyGMT figure.

#         Notes
#         -----
#         - Supports 8-region, 2-region, or full hemisphere plotting via `plot_regions`.
#         - For "diff" plots (where 'diff' in var_name), colors are assigned categorically (e.g., agreement/simulation/observation).
#         - Bathymetry and ice shelf overlays are loaded from NetCDF and shapefiles respectively.
#         - The method gracefully skips plotting if required coordinates or data are missing.

#         See Also
#         --------
#         - self.prepare_data_for_pygmt_plot : Prepares data dictionary for plotting
#         - self.create_cbar_frame           : Builds formatted colorbar strings
#         - self.load_GI_lon_lats            : Loads grounded iceberg locations
#         - self.load_ice_shelves            : Loads Antarctic ice shelf polygons
#         - self.load_IBCSO_bath             : Loads bathymetry grid from IBCSO

#         Examples
#         --------
#         >>> self.pygmt_map_plot_one_var(FIP_DA, "FIP", plot_regions=8, show_fig=True)
#         >>> self.pygmt_map_plot_one_var(diff_DA, "FIA_diff", diff_plot=True, cmap="categorical", plot_GI=True)
#         """
#         sim_name    = sim_name      if sim_name      is not None else self.sim_name
#         show_fig    = show_fig      if show_fig      is not None else self.show_fig
#         ow_fig      = overwrite_fig if overwrite_fig is not None else self.ow_fig
#         time_stamp  = time_stamp    if time_stamp    is not None else self.dt0_str
#         #lon_coord_name = lon_coord_name if lon_coord_name is not None else self.pygmt_dict.get("lon_coord_name", "TLON")
#         #lat_coord_name = lat_coord_name if lat_coord_name is not None else self.pygmt_dict.get("lat_coord_name", "TLAT")
#         cmap        = cmap          if cmap          is not None else self.plot_var_dict[var_name]['cmap']
#         series      = series        if series        is not None else self.plot_var_dict[var_name]['series']
#         reverse     = reverse       if reverse       is not None else self.plot_var_dict[var_name]['reverse']
#         cbar_lab    = cbar_label    if cbar_label    is not None else self.plot_var_dict[var_name]['name']
#         cbar_units  = cbar_units    if cbar_units    is not None else self.plot_var_dict[var_name]['units']
#         fig_size    = fig_size      if fig_size      is not None else self.pygmt_dict['fig_size']
#         cbar_pos    = cbar_position if cbar_position is not None else self.pygmt_dict['cbar_pos'].format(width=fig_size*0.8,height=0.75)
#         land_color  = land_color    if land_color    is not None else self.pygmt_dict['land_color']
#         water_color = water_color   if water_color   is not None else self.pygmt_dict['water_color']
#         if var_out is None:
#             var_out = var_name
#         if plot_iceshelves:
#             ANT_IS = self.load_ice_shelves()
#         if plot_bathymetry:
#             SO_BATH = self.load_IBCSO_bath()
#         if lon_coord_name is not None or lat_coord_name is not None:
#             plot_data_dict = self.prepare_data_for_pygmt_plot(da,
#                                                               bcoords        = False,
#                                                               tcoords        = False,
#                                                               lon_coord_name = lon_coord_name,
#                                                               lat_coord_name = lat_coord_name,
#                                                               diff_plot      = diff_plot)
#         else:
#             if use_bcoords and use_tcoords:
#                 raise ValueError("Cannot set both use_bcoords and use_tcoords to True")
#             plot_data_dict = self.prepare_data_for_pygmt_plot(da,
#                                                               bcoords        = use_bcoords,
#                                                               tcoords        = use_tcoords,
#                                                               lon_coord_name = None,
#                                                               lat_coord_name = None,
#                                                               diff_plot      = diff_plot)
#         required_keys = ['lon', 'lat', 'data']
#         try:
#             if not isinstance(plot_data_dict, dict):
#                 self.logger.warning("plot_data_dict is not a dictionary — skipping plot.")
#                 return
#             for k in required_keys:
#                 if k not in plot_data_dict:
#                     self.logger.warning(f"Missing key '{k}' in plot_data_dict — skipping plot.")
#                     return
#                 v = plot_data_dict[k]
#                 if v is None:
#                     self.logger.warning(f"plot_data_dict['{k}'] is None — skipping plot.")
#                     return
#                 if hasattr(v, "size") and v.size == 0:
#                     self.logger.warning(f"plot_data_dict['{k}'] is empty — skipping plot.")
#                     return
#         except Exception as e:
#             self.logger.warning(f"Skipping plot due to error: {e}")
#             return
#         cbar_frame = self.create_cbar_frame(series, cbar_lab, units=cbar_units, extend_cbar=extend_cbar)
#         hem_plot   = False
#         if plot_GI:
#             plot_GI_dict = self.load_GI_lon_lats()#plot_data_dict)
#         if plot_regions is not None and plot_regions==8:
#             self.logger.info("method will plot eight Antarctic sectors regional dictionary")
#             reg_dict = self.Ant_8sectors
#         elif plot_regions is not None and plot_regions==2:
#             self.logger.info("method will plot two Antarctic sectors regional dictionary")
#             reg_dict = self.Ant_2sectors
#         elif plot_regions is not None and regional_dict is not None:
#             self.logger.info("method will plot regional dictionary passed to this method")
#             reg_dict = regional_dict
#         elif plot_regions is not None and regional_dict is None:
#             self.logger.info("plot_regions argument not valid")
#         else:
#             self.logger.info("method will plot hemispheric data")
#             hem_plot = True
#             reg_dict = self.hemispheres_dict
#         if tit_str is not None:
#             basemap_frame = ["af", f"+t{tit_str}"]
#         else:
#             basemap_frame = ["af"]
#         for i, (reg_name, reg_vals) in enumerate(reg_dict.items()):
#             if hem_plot and reg_name!=hemisphere:
#                 continue
#             if P_png is None and self.save_fig:
#                 P_png = Path(self.D_graph, sim_name, reg_name, var_out, f"{time_stamp}_{sim_name}_{reg_name}_{var_out}.png")
#             region     = reg_vals['plot_region']
#             projection = reg_vals['projection']               
#             if hem_plot:
#                 projection = projection.format(fig_size=fig_size)
#             elif reg_name in list(self.Ant_8sectors.keys()):
#                 MC         = self.get_meridian_center_from_geographic_extent(region)
#                 projection = projection.format(MC=MC, fig_size=fig_size)
#             elif reg_name in list(self.Ant_2sectors.keys()):
#                 projection = projection.format(fig_width=fig_size)
#             fig = pygmt.Figure()
#             with pygmt.config(FONT_TITLE         = "16p,Courier-Bold",
#                               FONT_ANNOT_PRIMARY = "14p,Helvetica",
#                               COLOR_FOREGROUND   = 'black'):
#                 fig.basemap(region=region, projection=projection, frame=basemap_frame)
#                 if plot_bathymetry:
#                     fig.grdimage(grid=SO_BATH, cmap='geo')
#                 else:
#                     fig.coast(region=region, projection=projection, shorelines="1/0.5p,gray30", land=land_color, water=water_color)
#                 if "diff" in var_name.lower():
#                     lat        = da[lat_coord_name].values.flatten()
#                     lon        = da[lon_coord_name].values.flatten()
#                     val        = da.values.flatten()
#                     valid_mask = ~np.isnan(val)
#                     lat_valid  = lat[valid_mask]
#                     lon_valid  = lon[valid_mask]
#                     val_valid  = val[valid_mask].astype(int)
#                     label_map  = {1: "simulation", 0: "agreement", 2: "observation"}
#                     labels     = [label_map[v] for v in val_valid]
#                     df         = pd.DataFrame({"longitude" : lon_valid,
#                                             "latitude"  : lat_valid,
#                                             "z"         : val_valid.astype(int)})
#                     pygmt.makecpt(cmap="categorical", series=[0,2,1], color_model="+cagreement,simulation,observation")
#                     fig.plot(data=df, style=f"s{var_sq_size}c", cmap=True)
#                 elif "mask" in var_name.lower():
#                     fig.plot(x=plot_data_dict['lon'], y=plot_data_dict['lat'], fill='red', style=f"s{var_sq_size}c")
#                 else:
#                     pygmt.makecpt(cmap=cmap, reverse=reverse, series=series)
#                     fig.plot(x=plot_data_dict['lon'], y=plot_data_dict['lat'], fill=plot_data_dict['data'], style=f"s{var_sq_size}c", cmap=True)           
#                 if plot_bathymetry:
#                     fig.coast(region=region, projection=projection, shorelines="1/0.5p,gray30")
#                 if plot_GI:
#                     fig.plot(x=plot_GI_dict['lon'], y=plot_GI_dict['lat'], fill=GI_fill_color, style=f"c{GI_sq_size}c")
#                 if plot_iceshelves:
#                     fig.plot(data=ANT_IS, pen="0.2p,gray", fill="lightgray")
#                 if add_stat_annot:
#                     annot_text = self.generate_regional_annotation_stats(da, region, lon_coord_name, lat_coord_name)
#                     for i, line in enumerate(annot_text):
#                         try:
#                             fig.text(position="TR", text=line, font="12p,Helvetica-Bold,black", justify="LM", no_clip=True, offset=f"-1/{-0.5*i}")
#                         except pygmt.exceptions.GMTCLibError as e:
#                             self.logger.warning(f"Error in plotting anotation text {e} -- skipping annotation")
#                 if "diff" in var_name.lower():
#                     fig.colorbar(position=cbar_pos, frame=["x+l" + cbar_lab])
#                 elif "mask" not in var_name.lower():
#                     try:
#                         fig.colorbar(position=cbar_pos, frame=cbar_frame)
#                     except pygmt.exceptions.GMTCLibError as e:
#                         self.logger.warning(f"Error in adding colorbar: {e} — skipping colorbar.")
#             if P_png:
#                 if not P_png.exists():
#                     P_png.parent.mkdir(parents=True, exist_ok=True)
#                     fig.savefig(P_png)  
#                     self.logger.info(f"Saved figure to {P_png}")
#                 else:   
#                     if ow_fig:
#                         fig.savefig(P_png)  
#                         self.logger.info(f"Saved figure to {P_png}")
#                     else:
#                         self.logger.info(f"{P_png} already exists and not overwriting")
#                 P_png = None
#             if show_fig:
#                 fig.show()
#             pygmt.clib.Session.__exit__

#     def pygmt_2D_data_prep(self, da,
#                            x_coord_name = None,
#                            y_coord_name = None,
#                            region       = None,
#                            extra_mask   = None,
#                            mask_zero    = None):
#         da2 = da.squeeze(drop=True)
#         if da2.ndim != 2:
#             raise ValueError(f"Expected 2-D data after squeeze; got dims {da2.dims}")
#         x_name = x_coord_name or SI_tools.CICE_dict["lon_coord_name"]
#         y_name = y_coord_name or SI_tools.CICE_dict["lat_coord_name"]
#         xcoord = da2.coords[x_name]
#         ycoord = da2.coords[y_name]
#         if xcoord.ndim == 1 and ycoord.ndim == 1:
#             X2D, Y2D = np.meshgrid(xcoord.values, ycoord.values, indexing="xy")
#         elif xcoord.ndim == 2 and ycoord.ndim == 2:
#             if xcoord.dims != da2.dims or ycoord.dims != da2.dims:
#                 da2 = da2.transpose(*xcoord.dims)
#             X2D, Y2D = np.asarray(xcoord.values), np.asarray(ycoord.values)
#         else:
#             raise ValueError("Mixed coord dimensionality.")
#         Z2D = np.asarray(da2.values)
#         # --- decide zero-masking automatically if not specified ---
#         if mask_zero is None:
#             name = (da2.name or "").lower()
#             fv = str(da2.attrs.get("flag_values", "")).strip()
#             is_categorical = ("diff_cat" in name) or (fv in {"0 1 2", "0,1,2"})
#             mask_zero = not is_categorical
#         mask = np.isfinite(Z2D)
#         if mask_zero:
#             mask &= ~np.isclose(Z2D, 0.0, atol=1e-8)
#         if region is not None:
#             xmin, xmax, ymin, ymax = region
#             mask &= (X2D >= xmin) & (X2D <= xmax) & (Y2D >= ymin) & (Y2D <= ymax)
#         if extra_mask is not None:
#             em = np.asarray(extra_mask(da2))
#             if em.shape != Z2D.shape:
#                 raise ValueError("extra_mask shape mismatch")
#             mask &= em
#         x = X2D[mask]
#         y = Y2D[mask]
#         z = Z2D[mask]
#         return np.column_stack([x, y, z])

#     def pygmt_FIP_figure(self, plot_data,
#                         show_fig      = False,
#                         P_png         = None,
#                         region        = [0, 360, -90, -62],
#                         projection    = "S0.0/-90.0/50/20c",
#                         cmap          = "/g/data/gv90/da1339/GRAPHICAL/CPTs/AF2020_YlGnBu.cpt",
#                         reverse       = False,
#                         series        = [0, 1],
#                         cbar_pos      = "JBC+w15c/1c+mc+h",
#                         cbar_frame    = ["xa0.1f0.05+lFast Ice Persistence"],
#                         basemap_frame = ["af", f"+tAF2020-reG 2000-2018"],
#                         land_color    = "#666666",
#                         water_color   = "#BABCDE",
#                         G_pt_color    = "#FF8903",
#                         G_pt_marker   = "c",       # circle; for squares use "s"
#                         G_pt_size     = "0.05",    # point size (string, no unit)
#                         G_pt_unit     = "c",       # "c" for cm
#                         shoreline_pen = "1/0.25p",
#                         plot_GI         = False,
#                         GI_color        = "#E349D0",
#                         GI_marker       = "s",
#                         GI_size         = "0.01",
#                         GI_pt_unit      = "c",
#                         plot_bathymetry = False,
#                         # --- extras to help with DataArray inputs / categorical plotting ---
#                         var_name      = None,        # e.g., "diff_cat" to trigger categorical mode
#                         lat_coord_name = "lat",
#                         lon_coord_name = "lon",
#                         cat_cmap       = "hawaii",
#                         cat_labels    = ("union", "model dominate", "observation alone"),
#                         cat_series    = (0, 1, 2),   # codes that match cat_labels order
#                         cat_cbar_frame = None, #["+lDifference category"],
#                         # NEW:
#                         weight_da=None,                 # xr.DataArray matching plot_data shape (optional)
#                         weight_to_transparency=True):    # map weight∈[0,1] -> transparency∈[100,0]):
#         """
#         Plot FIP fields with PyGMT. If `var_name` contains 'diff_cat', uses a categorical CPT
#         with codes 0,1,2 => agreement, simulation, observation.

#         plot_data:
#         - xr.DataArray with 2-D 'lat' and 'lon' coords, or
#         - pd.DataFrame with columns ['longitude','latitude','z'].
#         """
#         # Prep GI overlay / bathy
#         if plot_GI:
#             SI_tools.load_bgrid(slice_hem=True)
#             GI_mask = (SI_tools.G_t["kmt_org"] == 1) & (SI_tools.G_t["kmt_mod"] == 0)
#             GI_lon  = SI_tools.normalise_longitudes(SI_tools.G_t["lon"].where(GI_mask.astype('float32')), "0-360").values.ravel()
#             GI_lat  = SI_tools.G_t["lat"].where(GI_mask.astype('float32')).values.ravel()
#         if plot_bathymetry:
#             SO_BATH = SI_tools.load_IBCSO_bath()
#         # Convert input to a plotting table
#         if isinstance(plot_data, xr.DataArray):
#             da = plot_data
#             # Use provided var_name if given, else fall back to DataArray.name
#             if var_name is None and hasattr(da, "name") and da.name is not None:
#                 var_name = da.name
#             lat = da[lat_coord_name].values.ravel()
#             lon = da[lon_coord_name].values.ravel()
#             val = da.values.ravel()
#             valid = np.isfinite(val)
#             df = pd.DataFrame( {"longitude": lon[valid], "latitude": lat[valid], "z": val[valid]})
#         elif isinstance(plot_data, pd.DataFrame):
#             df = plot_data.copy()
#         else:
#             raise TypeError("plot_data must be an xarray.DataArray or a pandas.DataFrame")
#         # Build the figure + base
#         fig = pygmt.Figure()
#         with pygmt.config(FONT_TITLE           = "22p,Bookman-Demi",
#                           FONT_ANNOT_PRIMARY   = "16p,NewCenturySchlbk-Roman",
#                           FONT_ANNOT_SECONDARY = "16p,NewCenturySchlbk-Bold",
#                           FONT_LABEL           = "16p,NewCenturySchlbk-Bold",
#                           COLOR_FOREGROUND     = "black",):
#             fig.basemap(region=region, projection=projection, frame=basemap_frame)
#             if plot_bathymetry:
#                 fig.grdimage(grid=SO_BATH, cmap="geo")
#             else:
#                 fig.coast(region=region, projection=projection, land=land_color, water=water_color)
#             # Categorical vs continuous color mapping
#             is_categorical = (var_name is not None) and ("diff_cat" in var_name.lower())
#             if is_categorical:
#                 df["z"] = df["z"].astype(int)
#                 # categorical CPT + legend
#                 codes = np.unique(df["z"].to_numpy())
#                 zmin, zmax = int(codes.min()), int(codes.max())
#                 if zmax <= zmin:
#                     zmax = zmin + 1
#                 pygmt.makecpt(cmap        = cat_cmap or "categorical",
#                               series      = [zmin, zmax, 1],             # e.g., 0..2 step 1
#                               color_model = "+c" + ",".join(cat_labels)) # e.g., "agreement,simulation,observation"
#                 if (weight_da is not None) and weight_to_transparency:
#                     # Align weight with the same mask used to build df
#                     if isinstance(plot_data, xr.DataArray):
#                         z_full = plot_data.values.ravel()
#                         valid  = np.isfinite(z_full)
#                         w_full = weight_da.squeeze(drop=True).values.ravel()
#                         wv     = w_full[valid]
#                     else:
#                         # If caller passed a DataFrame, require a 'weight' column
#                         if "weight" not in df.columns:
#                             raise ValueError("With DataFrame input, include a 'weight' column for opacity weighting.")
#                         wv = df["weight"].to_numpy()
#                     # Drop fully transparent points (weight<=0)
#                     keep = wv > 0
#                     if not np.all(keep):
#                         df = df.loc[keep].reset_index(drop=True)
#                         wv = wv[keep]
#                     # Bin weights → a few opacity levels (scalar transparency per batch)
#                     # Bins in weight-space; 4 bins is a good balance of speed and fidelity
#                     w_edges = np.array([0.00, 0.25, 0.50, 0.75, 1.01])
#                     bin_idx = np.digitize(wv, w_edges) - 1  # 0..3
#                     df["_bin"] = bin_idx
#                     # Representative transparency for each bin (tau = (1 - w_center)*100)
#                     tau_for_bin = np.array([87.5, 62.5, 37.5, 12.5], dtype=float)
#                     # Draw by category then by opacity bin (keeps CPT colors + avoids overlap artifacts)
#                     for k in sorted(codes):
#                         sel_k = (df["z"].to_numpy() == int(k))
#                         if not np.any(sel_k):
#                             continue
#                         for b in range(len(tau_for_bin)):
#                             sel = sel_k & (df["_bin"].to_numpy() == b)
#                             if not np.any(sel):
#                                 continue
#                             fig.plot(data         = df.loc[sel],
#                                      style        = f"{G_pt_marker}{G_pt_size}{G_pt_unit}",
#                                      cmap         = True,
#                                      transparency = float(tau_for_bin[b]))   # <-- SCALAR, not array
#                     # Clean up helper column
#                     df.drop(columns=["_bin"], inplace=True)
#                 else:
#                     # No weighting: single draw
#                     fig.plot(data  = df,
#                              style = f"{G_pt_marker}{G_pt_size}{G_pt_unit}",
#                              cmap  = True)
#                 fig.colorbar(position=cbar_pos, frame=cat_cbar_frame)
#             else:
#                 # Continuous CPT (default behavior)
#                 pygmt.makecpt(cmap=cmap, reverse=reverse, series=series)
#                 fig.plot(data  = df,
#                          style = f"{G_pt_marker}{G_pt_size}{G_pt_unit}",
#                          cmap  = True)
#                 fig.colorbar(position=cbar_pos, frame=cbar_frame)
#             if plot_GI:
#                 fig.plot(x=GI_lon, y=GI_lat, fill=GI_color, style=f"{GI_marker}{GI_size}{GI_pt_unit}")
#             fig.coast(region=region, projection=projection, shorelines=shoreline_pen)
#         if P_png:
#             fig.savefig(P_png)
#         if show_fig:
#             fig.show()
#         return fig

#     def _smooth_df_time(self, df: pd.DataFrame, window, center=True, min_periods=1):
#         """
#         Smooth a time series in df with columns ['time','data'].
#         'window' can be an int (# of samples) or a time offset string like '15D'.
#         Returns a new df with smoothed 'data'.
#         """
#         if window is None:
#             return df
#         s = df.set_index("time")["data"]
#         # Pandas rolling supports both integer and time-based windows on a DatetimeIndex
#         s_sm = s.rolling(window=window, center=center, min_periods=min_periods).mean()
#         out = df.copy()
#         out["data"] = s_sm.values
#         return out

#     def _repeat_doy_over_range(self, df_time_data, full_range, time_coord="time"):
#         src         = df_time_data.copy()
#         src["doy"]  = src[time_coord].dt.dayofyear
#         clim_lookup = src.groupby("doy")["data"].mean()  # index = 1..365 or 366
#         rep         = pd.DataFrame({"time": full_range})
#         rep["doy"]  = rep["time"].dt.dayofyear
#         # handle leap day: if 366 not in lookup, remap DOY=366 -> 365
#         if 366 not in clim_lookup.index:
#             rep["doy_eff"] = rep["doy"].where(rep["doy"] != 366, 365)
#         else:
#             rep["doy_eff"] = rep["doy"]
#         rep["rep"] = rep["doy_eff"].map(clim_lookup)
#         return rep[["time", "rep"]]

#     def pygmt_timeseries(self, ts_dict,
#                         comp_name        : str              = "test",
#                         primary_key      : str              = "FIA",  # "FIA"
#                         smooth           : str|int | None   = None,   # e.g., 15, "15D"
#                         clim_smooth      : int | None       = None,   # 15
#                         climatology      : bool             = False,
#                         ylabel           : str              = None,   # "@[Fast Ice Area (1\\times10^3 km^2)@[",
#                         ylim             : tuple            = [0,1000],
#                         yaxis_pri        : int              = None,
#                         ytick_pri        : int              = 100,
#                         ytick_sec        : int              = 50,
#                         projection       : str              = None,
#                         fig_width        : str              = None,
#                         fig_height       : str              = None,
#                         xaxis_pri        : str              = None,
#                         xaxis_sec        : str              = None,
#                         frame_bndy       : str              = "WS",
#                         legend_pos       : str              = None,
#                         legend_box       : str              = "+gwhite+p0.5p",
#                         fmt_dt_pri       : str              = None,
#                         fmt_dt_sec       : str              = None,
#                         fmt_dt_map       : str              = None,
#                         fnt_type         : str              = "Helvetica",
#                         fnt_wght_lab     : str              = "20p",
#                         fnt_wght_ax      : str              = "18p",
#                         line_pen         : str              = "1p",
#                         grid_wght_pri    : str              = ".25p",
#                         grid_wght_sec    : str              = ".1p",
#                         P_png            : str              = None,
#                         time_coord       : str              = "time",
#                         time_coord_alt   : str              = "date",
#                         keys2plot        : list             = None,
#                         repeat_keys      : list[str] | None = None,
#                         repeat_policy    : str              = "inside_others", # "inside_others" | "outside_others" | "fill_gaps" | "always"
#                         repeat_ref_keys  : list[str] | None = None,            # which keys define the "others" window; default: all except current
#                         clip_x_axis      : bool             = False,
#                         zero_line        : bool             = False,
#                         zero_line_level  : float            = 0.0,
#                         zero_line_pen    : str              = "2p,black",
#                         save_fig         : bool             = None,
#                         show_fig         : bool             = None):
#         """
#         Plot time series of a primary variable (e.g., FIA) for a set of simulations or observations.

#         Parameters
#         ----------
#         ts_dict : dict
#             Dictionary of xarray DataArrays keyed by simulation or dataset name.
#         comp_name : str
#             Name for the comparison (used in figure title and filename).
#         primary_key : str
#             Key used to extract the variable from each dataset (except 'AF2020').
#         climatology : bool
#             If True, plot daily climatology with fill and mean lines.
#         ylim : tuple or None
#             Y-axis limits. If None, inferred from data with 5% padding.
#         ylabel : str
#             Label for the Y-axis.
#         ytick_inc : int
#             Interval between Y-axis ticks.
#         xaxis_pri, xaxis_sec : str
#             GMT frame settings for primary and secondary axes (when climatology is False).
#         P_png : str or None
#             Optional full path to save figure. If None, default filename is constructed.
#         legend_box : str
#             GMT legend box styling.
#         line_pen : str
#             Line thickness for plotted time series.
#         time_coord : str
#             Name of time coordinate in each DataArray.
#         keys2plot : list or None
#             If provided, only datasets with keys in this list are plotted.
#         show_fig : bool or None
#             If True, show figure interactively. Defaults to self.show_fig.
#         repeat_keys : list[str] or None
#             Keys in `ts_dict` to plot as a repeated day-of-year climatology when
#             `climatology=False`.
#         repeat_policy : {"outside_others","fill_gaps","always"}
#             - "outside_others": keep original values where *other* series exist in time,
#             but replace values outside that union with day-of-year climatology of
#             the nominated series (good for fair comparison).
#             - "fill_gaps": use climatology only where the nominated series has no data,
#             keep original where it does.
#             - "always": ignore original values and plot the repeated climatology across
#             the full x-range.
#         repeat_ref_keys : list[str] or None
#             If provided, these keys define the "others" time span for the
#             "outside_others" policy. By default, it's all plotted series except the
#             current one.
#         """
#         show_fig   = show_fig   if show_fig   is not None else self.show_fig
#         save_fig   = save_fig   if save_fig   is not None else self.save_fig
#         fmt_dt_pri = fmt_dt_pri if fmt_dt_pri is not None else "Character"
#         fmt_dt_sec = fmt_dt_sec if fmt_dt_sec is not None else "Abbreviated"
#         fmt_dt_map = fmt_dt_map if fmt_dt_map is not None else "o"
#         # primary_key = primary_key if primary_key is not None else keys2plot[0]
#         # if keys2plot is None:
#         #     raise("either 'primary_key' or 'keys2plot' must be defined")
#         # need to get out the maximum times for plot boundaries
#         tmin, tmax = self.extract_min_max_dates(ts_dict, keys2plot=keys2plot, primary_key=primary_key, time_coord=time_coord)
#         # --- normalize new params ---
#         if isinstance(repeat_keys, (str, bytes)):
#             repeat_keys = [repeat_keys]
#         repeat_keys = set(repeat_keys or [])
#         if isinstance(repeat_ref_keys, (str, bytes)):
#             repeat_ref_keys = [repeat_ref_keys]
#         x_start, x_end = tmin, tmax
#         if (not climatology) and clip_x_axis and repeat_keys:
#             # reference keys = all plotted except the repeated ones, unless user specified
#             if repeat_ref_keys:
#                 ref_keys = set(repeat_ref_keys)
#             else:
#                 ref_keys = {k for k in ts_dict.keys()
#                             if (keys2plot is None or k in (keys2plot or []))
#                             and k not in repeat_keys}

#             other_mins, other_maxs = [], []
#             for k in ref_keys:
#                 da2 = ts_dict[k][primary_key]
#                 tt2 = pd.to_datetime(da2[time_coord].values)
#                 if len(tt2) > 0:
#                     other_mins.append(tt2.min()); other_maxs.append(tt2.max())

#             if other_mins:  # only override if we actually have refs
#                 x_start = pd.to_datetime(min(other_mins)).normalize()
#                 x_end   = pd.to_datetime(max(other_maxs)).normalize()
#         # there are differences in the projection and x-axis for the two types of figures
#         if climatology:
#             fake_year  = 1996
#             fig_width  = fig_width  if fig_width  is not None else "20c"
#             fig_height = fig_height if fig_height is not None else "15c"
#             xaxis_sec  = xaxis_sec  if xaxis_sec  is not None else None
#             xaxis_pri  = xaxis_pri  if xaxis_pri  is not None else "a1Og"
#             region     = [f"{fake_year}-01-01", f"{fake_year}-12-31", ylim[0], ylim[1]]
#         else:
#             fig_width  = fig_width  if fig_width  is not None else "50c"
#             fig_height = fig_height if fig_height is not None else "15c"
#             xaxis_pri  = xaxis_pri  if xaxis_pri  is not None else "a2Of30D"#g30D"
#             xaxis_sec  = xaxis_sec  if xaxis_sec  is not None else "a1Y"
#             region     = [x_start.strftime("%Y-%m-%d"), x_end.strftime("%Y-%m-%d"), ylim[0], ylim[1]]
#         # define the projection and frame of the figure
#         legend_pos = legend_pos if legend_pos is not None else f"JTL+jTL+o0.2c+w{fig_width}"
#         projection = f"X{fig_width}/{fig_height}"
#         yaxis_pri  = yaxis_pri if yaxis_pri is not None else f"a{ytick_pri}f{ytick_sec}+l{ylabel}"
#         if xaxis_sec is not None:
#             frame = [frame_bndy, f"sx{xaxis_sec}", f"px{xaxis_pri}", f"py{yaxis_pri}"]
#         else:
#             frame = [frame_bndy, f"px{xaxis_pri}", f"py{yaxis_pri}"]
#         # make sure the figure has the same configuration by wrapping it in a with condition
#         fig = pygmt.Figure()
#         with pygmt.config(FONT_LABEL                = f"{fnt_wght_lab},{fnt_type}",
#                           FONT                      = f"{fnt_wght_ax},{fnt_type}",
#                           MAP_GRID_PEN_PRIMARY      = grid_wght_pri,
#                           MAP_GRID_PEN_SECONDARY    = grid_wght_sec,
#                           FORMAT_TIME_PRIMARY_MAP   = fmt_dt_pri,
#                           FORMAT_TIME_SECONDARY_MAP = fmt_dt_sec,
#                           FORMAT_DATE_MAP           = fmt_dt_map):
#             fig.basemap(projection=projection, region=region, frame=frame)
#             # loop over each key in the dictionary and if keys2plot is defined only plot those dictionaries
#             cnt=0
#             for i, (dict_key, data) in enumerate(ts_dict.items()):
#                 if data.get("line_clr", None) is not None:
#                     line_color = data['line_clr']
#                 else:
#                     line_color = self.plot_var_dict.get(dict_key, {}).get("line_clr", f"C{i}")
#                 if data.get("leg_lab", None) is not None:
#                     leg_lab = data['leg_lab']
#                 else:
#                     leg_lab = self.plot_var_dict.get(dict_key, {}).get("leg_lab", dict_key )
#                 da = data[primary_key]
#                 self.logger.info(f"pulling out data array for {dict_key} and putting into dataframe")
#                 self.logger.info(f"legend label: {leg_lab}")
#                 self.logger.info(f"line color  : {line_color}")
#                 if dict_key=="AF2020" and primary_key=="FIA":
#                     df = pd.DataFrame({"time": pd.to_datetime(da[time_coord_alt].values), "data": da.values})
#                 else:
#                     df = pd.DataFrame({"time": pd.to_datetime(da[time_coord].values), "data": da.values})
#                 if climatology:
#                     if dict_key=="AF2020" and primary_key=="FIA":
#                         clim = self.compute_doy_climatology(da, time_coord=time_coord_alt)
#                     else:
#                         clim = self.compute_doy_climatology(da)
#                     mean_x = clim['mean'].index
#                     mean_y = clim['mean'].values
#                     if clim_smooth is not None and clim_smooth > 1:
#                         mean_y = (pd.Series(mean_y, index=mean_x)
#                                     .rolling(window=clim_smooth, center=True, min_periods=1)
#                                     .mean()
#                                     .values)
#                     fig.plot(x            = np.concatenate([clim['min'].index, clim['max'].index[::-1]]),
#                              y            = np.concatenate([clim['min'].values, clim['max'].values[::-1]]),
#                              fill         = f"{line_color}@80",
#                              close        = True,
#                              transparency = 80)
#                     fig.plot(x     = mean_x,
#                              y     = mean_y,
#                              pen   = f"{line_pen},{line_color}",
#                              label = leg_lab)
#                 else:       
#                     if dict_key in repeat_keys:
#                         # target x-range (union for figure)
#                         t_start, t_end = tmin.normalize(), tmax.normalize()
#                         full_range = pd.date_range(t_start, t_end, freq="D")
#                         # repeated climatology over the full range
#                         rep_df = self._repeat_doy_over_range(df, full_range)  # cols: time, rep
#                         # original reindexed over full range
#                         base = pd.DataFrame({"time": full_range}).merge(df.rename(columns={"data": "orig"}), on="time", how="left")
#                         mix = base.merge(rep_df, on="time", how="left")
#                         if repeat_policy == "always":
#                             out_df = mix[["time", "rep"]].rename(columns={"rep": "data"})
#                         elif repeat_policy == "fill_gaps":
#                             mix["data"] = mix["orig"].fillna(mix["rep"])
#                             out_df = mix[["time", "data"]]
#                         elif repeat_policy in ("outside_others", "inside_others"):
#                             # Determine the "others" window
#                             if repeat_ref_keys is not None:
#                                 ref_keys = set(repeat_ref_keys)
#                             else:
#                                 ref_keys = set(k for k in ts_dict.keys()
#                                             if k != dict_key and (keys2plot is None or k in keys2plot))
#                             if len(ref_keys) == 0:
#                                 # no refs -> behave like fill_gaps (but still allow inside_others to clip nothing)
#                                 mix["data"] = mix["orig"].fillna(mix["rep"])
#                                 if repeat_policy == "inside_others":
#                                     mix.loc[:, "data"] = np.nan  # nothing to show if no reference window
#                                 out_df = mix[["time", "data"]]
#                             else:
#                                 # compute min/max over "others"
#                                 other_mins, other_maxs = [], []
#                                 for k in ref_keys:
#                                     d2 = ts_dict[k][primary_key]
#                                     tt = pd.to_datetime(d2[time_coord].values)
#                                     if len(tt) > 0:
#                                         other_mins.append(tt.min()); other_maxs.append(tt.max())
#                                 other_min = pd.to_datetime(min(other_mins)).normalize()
#                                 other_max = pd.to_datetime(max(other_maxs)).normalize()
#                                 inside = (mix["time"] >= other_min) & (mix["time"] <= other_max)
#                                 # inside the others’ window: use original where present, else repeated
#                                 mix.loc[inside, "data"] = mix.loc[inside, "orig"].fillna(mix.loc[inside, "rep"])
#                                 if repeat_policy == "inside_others":
#                                     # outside -> hide (clip)
#                                     mix.loc[~inside, "data"] = np.nan
#                                 else:
#                                     # outside -> use repeated climatology
#                                     mix.loc[~inside, "data"] = mix.loc[~inside, "rep"]
#                                 out_df = mix[["time", "data"]]
#                         else:
#                             raise ValueError(f"Unknown repeat_policy: {repeat_policy}")
#                         # optional smoothing after composition
#                         out_df = self._smooth_df_time(out_df, window=smooth)
#                         fig.plot(x=out_df["time"], y=out_df["data"], pen=f"{line_pen},{line_color}", label=leg_lab)
#                     else:
#                         df_sm = self._smooth_df_time(df, window=smooth)
#                         fig.plot(x=df_sm["time"], y=df_sm["data"], pen=f"{line_pen},{line_color}", label=leg_lab)
#             if zero_line:
#                 y0 = float(zero_line_level)
#                 y_min, y_max = ylim
#                 if (y_min <= y0 <= y_max):
#                     x0, x1 = region[0], region[1]
#                     fig.plot(x=[x0, x1], y=[y0, y0], pen=zero_line_pen)
#             fig.legend(position=legend_pos, box=legend_box)
#         if save_fig:
#             F_png = f"{primary_key}_{comp_name}_{'climatology_' if climatology else ''}{tmin.strftime('%Y')}-{tmax.strftime('%Y')}.png"
#             P_png = P_png if P_png is not None else Path(self.D_graph, "timeseries", F_png)
#             fig.savefig(P_png, dpi=300)
#             self.logger.info(f"saved figure to {P_png}")
#         if show_fig:
#             fig.show()
#         pygmt.clib.Session.__exit__

#     def plot_taylor(self, stats_dict, out_path):
#         fig = plt.figure(figsize=(6, 6))
#         ax = fig.add_subplot(111, polar=True)
#         ax.set_theta_zero_location('E')
#         ax.set_theta_direction(-1)

#         for rms in [0.2, 0.5, 1.0, 1.5]:
#             rs = np.linspace(0.5, 2.0, 500)
#             theta = np.arccos(np.clip(1 - (rms**2 - 1)/(2 * rs), -1, 1))
#             ax.plot(theta, rs, '--', color='gray', lw=0.6)
#         ax.plot([0], [1], 'ko', label='Reference')

#         for label, stat in stats_dict.items():
#             angle = np.arccos(stat["corr"])
#             r = stat["std_ratio"]
#             ax.plot(angle, r, 'o', label=label)

#         ax.set_rmax(2)
#         ax.set_rticks([0.5, 1.0, 1.5, 2.0])
#         ax.set_rlabel_position(135)
#         ax.set_title("Taylor Diagram: Sea Ice Speed Comparison", fontsize=12)
#         ax.legend(loc='upper right', bbox_to_anchor=(1.45, 1.1))
#         plt.tight_layout()
#         plt.savefig(out_path)
#         plt.close()