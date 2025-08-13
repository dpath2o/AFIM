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

    def prepare_data_for_pygmt_plot(self, da, bcoords=False, diff_plot=False):
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
        data_dict = {}
        self.load_bgrid(slice_hem=True)
        if bcoords:
            lon2d = self.G_u['lon'].values
            lat2d = self.G_u['lat'].values
        else:
            lon2d = self.G_t['lon'].values
            lat2d = self.G_t['lat'].values
        data2d  = np.asarray(da.data).astype('float32')
        if diff_plot:
            mask = (data2d >= -1) & (data2d <= 1) & np.isfinite(data2d)
        else:
            mask = np.isfinite(data2d)# & (data2d > 0)
        data_dict['data']  = data2d[mask].ravel()
        data_dict['lon']   = lon2d[mask].ravel()
        data_dict['lat']   = lat2d[mask].ravel()
        return data_dict

    def load_GI_lon_lats(self):
        """
        Extract longitude and latitude positions of grounded iceberg (GI) grid cells.

        This method identifies cells where the landmask has been modified to add grounded icebergs
        (i.e., locations where the original landmask had ocean (`kmt_org == 1`) and the modified
        landmask has land (`kmt_mod == 0`)), and returns their corresponding geographic coordinates.

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
        self.load_bgrid(slice_hem=True)
        GI_loc_dict        = {}
        GI_mask            = (self.G_t['kmt_org'] == 1) & (self.G_t['kmt_mod']== 0)
        GI_loc_dict['lon'] = self.G_t['lon'][GI_mask].ravel()
        GI_loc_dict['lat'] = self.G_t['lat'][GI_mask].ravel()
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
                            GI_fill_color     : str   = "red",
                            plot_GI           : bool  = False,
                            min_max_trans_val : int   = 80,
                            yshift_top        : str   = None,
                            yshift_bot        : str   = None,
                            bottom_frame_bndy : str   = "WSne",
                            bottom_yaxis      : str   = None,
                            bottom_xaxis      : str   = None,
                            land_clr          : str   = 'gray',
                            water_clr         : str   = "white",
                            coast_pen         : str   = "1/0.5p,gray30",
                            cbar_pos          : str   = None,
                            lon_coord_name    : str   = None,
                            lat_coord_name    : str   = None,
                            cmap              : str   = None,
                            series            : list  = None,
                            cbar_frame        : str   = None,
                            ANT_IS_pen        : str   = "0.2p,gray",
                            ANT_IS_color      : str   = "lightgray",
                            font_annot_pri    : str   = "24p,Times-Roman",
                            font_annot_sec    : str   = "16p,Times-Roman",
                            font_lab          : str   = "22p,Times-Bold",
                            line_pen         : str   = "2p",
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
                fig.coast(land=land_clr, water=water_clr, shorelines=coast_pen)
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
                fig.plot(data=ANT_IS, pen=ANT_IS_pen, fill=ANT_IS_color)
            fig.colorbar(position=cbar_pos, frame=cbar_frame)
        # -------------------------------------------------------------------------
        # Save / Show
        # -------------------------------------------------------------------------
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
        plot_data_dict = self.prepare_data_for_pygmt_plot(da, diff_plot=diff_plot)
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

    def pygmt_timeseries(self, ts_dict,
                        comp_name    : str   = "test",
                        primary_key  : str   = "FIA",
                        smooth       : str|int|None = None,   # NEW: e.g., 15, "15D"
                        clim_smooth  : int|None = None, 
                        climatology  : bool  = False,
                        ylabel       : str   = "@[Fast Ice Area (1\\times10^3 km^2)@[",
                        ylim         : tuple = [0,1000],
                        yaxis_pri    : int   = None,
                        ytick_pri    : int   = 100,
                        ytick_sec    : int   = 50,
                        projection   : str   = None,
                        fig_width    : str   = None,
                        fig_height   : str   = None,
                        xaxis_pri    : str   = None,
                        xaxis_sec    : str   = None,
                        frame_bndy   : str   = "WS",
                        legend_pos   : str   = None,
                        legend_box   : str   = "+gwhite+p0.5p",
                        fmt_dt_pri   : str   = None,
                        fmt_dt_sec   : str   = None,
                        fmt_dt_map   : str   = None,
                        fnt_type     : str   = "Helvetica",
                        fnt_wght_lab : str   = "20p",
                        fnt_wght_ax  : str   = "18p",
                        line_pen    : str    = "1p",
                        grid_wght_pri: str   = ".25p",
                        grid_wght_sec: str   = ".1p",
                        P_png        : str   = None,
                        time_coord   : str   = "time",
                        keys2plot    : list  = None,
                        save_fig     : bool  = None,
                        show_fig     : bool  = None):
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
        """
        show_fig   = show_fig   if show_fig   is not None else self.show_fig
        save_fig   = save_fig   if save_fig   is not None else self.save_fig
        fmt_dt_pri = fmt_dt_pri if fmt_dt_pri is not None else "Character"
        fmt_dt_sec = fmt_dt_sec if fmt_dt_sec is not None else "Abbreviated"
        fmt_dt_map = fmt_dt_map if fmt_dt_map is not None else "o"
        # need to get out the maximum times for plot boundaries
        tmin, tmax = self.extract_min_max_dates(ts_dict, keys2plot=keys2plot, primary_key=primary_key, time_coord=time_coord)
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
            xaxis_pri  = xaxis_pri  if xaxis_pri  is not None else "a2Of30Dg30D"
            xaxis_sec  = xaxis_sec  if xaxis_sec  is not None else "a1Y"
            region     = [tmin.strftime("%Y-%m-%d"), tmax.strftime("%Y-%m-%d"), ylim[0], ylim[1]]
        # define the projection and frame of the figure
        legend_pos = legend_pos if legend_pos is not None else f"JTL+jTL+o0.2c+w{fig_width}"
        projection = f"X{fig_width}/{fig_height}"
        yaxis_pri  = yaxis_pri if yaxis_pri is not None else f"a{ytick_pri}gf{ytick_sec}+l{ylabel}"
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
                          #MAP_TICK_LENGTH_PRIMARY = "0.1c",
                          FORMAT_TIME_PRIMARY_MAP   = fmt_dt_pri,
                          FORMAT_TIME_SECONDARY_MAP = fmt_dt_sec,
                          FORMAT_DATE_MAP           = fmt_dt_map):
            fig.basemap(projection=projection, region=region, frame=frame)
            # loop over each key in the dictionary and if keys2plot is defined only plot those dictionaries
            cnt=0
            for i, (dict_key, data) in enumerate(ts_dict.items()):
                line_color = self.plot_var_dict.get(dict_key, {}).get("line_clr", f"C{i}")
                leg_lab    = self.plot_var_dict.get(dict_key, {}).get("leg_lab", dict_key )
                da = data[primary_key]
                self.logger.info(f"pulling out data array for {dict_key} and putting into dataframe")
                self.logger.info(f"legend label: {leg_lab}")
                self.logger.info(f"line color  : {line_color}")
                df = pd.DataFrame({"time": pd.to_datetime(da[time_coord].values), "data": da.values})
                if climatology:
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
                             y     = mean_y.values,
                             pen   = f"{line_pen},{line_color}",
                             label = leg_lab)
                else:
                # Repeat AF2020 across the full date range if needed
                    if (dict_key=='AF2020') and primary_key=="FIA":
                        t_start, t_end = tmin.normalize(), tmax.normalize()
                        full_range = pd.date_range(t_start, t_end, freq='D')
                        # Build a repeated time series using day-of-year as lookup
                        df["doy"] = df["time"].dt.dayofyear
                        #print(df)
                        #clim_lookup = df.set_index("doy")["data"]
                        clim_lookup = df.groupby("doy")["data"].mean()
                        repeated_df = pd.DataFrame({"time": full_range})
                        repeated_df["doy"] = repeated_df["time"].dt.dayofyear
                        #print(repeated_df)
                        repeated_df["data"] = repeated_df["doy"].map(clim_lookup)
                        repeated_df = _smooth_df_time(repeated_df, window=smooth)
                        fig.plot(x=repeated_df["time"], y=repeated_df["data"], pen=f"{line_pen},{line_color}", label=leg_lab)
                    else:
                        df_sm = _smooth_df_time(df, window=smooth)
                        fig.plot(x=df_sm["time"], y=df_sm["data"], pen=f"{line_pen},{line_color}", label=leg_lab)                
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