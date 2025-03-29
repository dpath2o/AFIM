import xarray as xr
import pandas as pd
import numpy as np
import json
from pathlib import Path
from datetime import datetime
import pygmt
import os
import time
import sys
import matplotlib.pyplot as plt
sys.path.insert(0, '/home/581/da1339/AFIM/src/AFIM/src')
from grounded_iceberg_processor import GroundedIcebergProcessor
import geopandas as gpd

class SeaIcePlotter:
    """
    Class to plot fast ice and pack ice fields from AFIM simulations.

    This class supports reading spatial and temporal Zarr datasets
    for fast or pack ice and generates static maps or time series
    visualisations using PyGMT.

    Parameters
    ----------
    sim_name       : str
                     Name of the simulation (used to locate output directories).
    dt0_str        : str
                     Start date string (format: YYYY-MM-DD).
    dtN_str        : str, optional
                     End date string. Required for multi-date plots.
    json_path      : str or Path, optional
                     Path to the AFIM project JSON configuration.
    hemisphere     : str, optional
                     Hemisphere to focus on ('north' or 'south'). Defaults to config setting.
    overwrite      : bool, optional
                     If True, overwrite existing figures (default: False).
    show_figs      : bool, optional
                     If True, display figures after creation (default: False).
    prepare        : bool, optional
                     Whether to load plotting data during initialisation (default: True).
    single_figure  : bool, optional
                     Treat input date as single date, even if dtN_str is passed (default: False).
    ice_type       : str, optional
                     Type of ice being plotted ('fast' or 'pack'). Affects file/path resolution.

    Attributes
    ----------
    config         : dict
                     Full project configuration dictionary from JSON.
    P_zarr         : Path
                     Path to the main Zarr dataset being plotted.
    plot_var_dict  : dict
                     Plotting configuration for variables (colorbars, units, ranges).
    plot_df        : pandas.DataFrame
                     Flat dataframe of points to plot (lat/lon/dat triples).
    plot_region    : str
                     Region name if plotting subregions.

    Examples
    --------
    >>> from sea_ice_plotter import SeaIcePlotter
    >>> plotter = SeaIcePlotter("Rothrock", "1994-09-08", ice_type="fast")
    >>> plotter.plot_map(region_name="Weddell", text_str="Fast Ice")
    """

    def __init__(self, sim_name, dt0_str, dtN_str=None, json_path=None, hemisphere=None, overwrite=False,
                 show_figs=False, prepare=True, single_figure=False, ice_type='FI'):
        """
        Initialise the SeaIcePlotter class and resolve data paths and plotting config.

        Sets internal plotting flags and loads project-level JSON config.
        Determines correct Zarr path(s) based on time and ice type.

        Parameters
        ----------
        sim_name      : str
                       Simulation name (used for directory and file lookup).
        dt0_str       : str
                       Date string of first or only timepoint (YYYY-MM-DD).
        dtN_str       : str, optional
                       Optional end date string (YYYY-MM-DD).
        json_path     : str, optional
                       Path to AFIM JSON config. If None, defaults to hardcoded path.
        hemisphere    : str, optional
                       Which hemisphere to plot. Defaults to config value.
        overwrite     : bool
                       Whether to overwrite existing figures.
        show_figs     : bool
                       Whether to show figures interactively (for notebooks).
        prepare       : bool
                       Whether to load Zarr data immediately.
        single_figure : bool
                       Whether to treat as a single frame (even if dtN_str exists).
        ice_type      : str
                       "fast" or "pack" to determine data subdirectories.
        """
        self.sim_name    = sim_name
        self.ice_type    = ice_type.upper()  # 'FI' or 'PI'
        if json_path is None:
            json_path = "/home/581/da1339/AFIM/src/AFIM/JSONs/afim_cice_analysis.json"
        with open(json_path, 'r') as f:
            self.config = json.load(f)

        self.ow_fig         = overwrite
        self.show_figs      = show_figs
        self.hemisphere     = hemisphere or self.config.get('hemisphere', 'south')
        self.plot_var_dict  = self.config["plot_var_dict"]
        self.pygmt_dict     = self.config["pygmt_dict"]
        self.sim_dict       = self.config["sim_dict"]
        self.CICE_dict      = self.config["CICE_dict"]
        self.GI_dict        = self.config["GI_dict"]
        self.default_var    = self.config.get(f"default_plot_var_{self.ice_type}", f"{self.ice_type}C")

        self.D_FI    = Path(self.config['D_dict']['AFIM_out'], sim_name, self.ice_type)
        self.D_graph = Path(self.config['D_dict']['graph'], 'AFIM', sim_name, self.ice_type)
        self.D_ts    = Path(self.config['D_dict']['graph'], 'timeseries')

        self.gi_processor   = GroundedIcebergProcessor(self.config, sim_name)
        self.use_gi         = self.gi_processor.use_gi
        if self.use_gi:
            self.gi_processor.load_AFIM_GI()
            self.GI_total_area  = self.gi_processor.total_area
            self.KMT_path       = self.gi_processor.KMT_path
        else:
            self.GI_total_area = 0
            self.KMT_path      = self.CICE_dict['P_KMT']

        self.define_hemisphere(self.hemisphere)

        if dtN_str is None or single_figure:
            P_zarr, self.plt_dt_str = self._find_nearest_zarr_file(dt0_str)
            print(f"Using nearest Zarr file: {P_zarr.name}")
            self.P_zarr  = P_zarr
            self.dt0_str = self.plt_dt_str
            self.dtN_str = None
        else:
            zarr_pairs = self._find_nearest_zarr_file(dt0_str, end_date_str=dtN_str)
            self.P_zarr_list     = [p for p, _ in zarr_pairs]
            self.plt_dt_str_list = [d for _, d in zarr_pairs]
            self.dt0_str = self.plt_dt_str_list[0]
            self.dtN_str = self.plt_dt_str_list[-1]

        if prepare and dtN_str is None:
            self.prepare_data_for_plotting()

    def define_hemisphere(self, hemisphere):
        if hemisphere.lower().startswith("n"):
            self.hemisphere_projection = "S0.0/90.0/50/15c"
            self.hemisphere_map_extent = [-180, 180, 55, 90]
            self.hemisphere_map_text_location = [-120, 56]
        else:
            self.hemisphere_projection = "S0.0/-90.0/50/15c"
            self.hemisphere_map_extent = [-180, 180, -90, -55]
            self.hemisphere_map_text_location = [0, -89]

    def _find_nearest_zarr_file(self, start_date_str, end_date_str=None):
        """
        Find nearest matching Zarr file(s) to a given start (and optional end) date.

        This method searches the appropriate ice_type subdirectory for Zarr files
        matching the required date(s). If a range is provided, it returns all Zarr
        files within that range. If only one date is provided, it returns the single
        closest file.

        Parameters
        ----------
        start_date_str : str
                         The desired starting date (YYYY-MM-DD).
        end_date_str   : str, optional
                         Optional end date. If provided, will return all Zarrs between.

        Returns
        -------
        Path or list of (Path, str)
            The closest file path (and its matched date string) or list of such pairs.
        """
        if self.ice_type.upper()=='FI':
            F_prefix = "fast_ice"
        elif self.ice_type.upper()=='PI':
            F_prefix = "pack_ice"
        def extract_date(f):
            return datetime.strptime(f.stem.replace(f"{F_prefix}_ice_", ""), "%Y-%m-%d")
        all_files = sorted(self.D_FI.glob(f"{F_prefix}_ice_*.zarr"))
        if not all_files:
            raise FileNotFoundError(f"No Zarr files found in {self.D_FI}")
        file_dates = [(f, extract_date(f)) for f in all_files]
        if end_date_str is None:
            target_date = datetime.strptime(start_date_str, "%Y-%m-%d")
            closest_file, matched_date = min(file_dates, key=lambda x: abs(x[1] - target_date))
            return closest_file, matched_date.strftime("%Y-%m-%d")
        start_dt = datetime.strptime(start_date_str, "%Y-%m-%d")
        end_dt   = datetime.strptime(end_date_str, "%Y-%m-%d")
        in_range = [(f, d.strftime("%Y-%m-%d")) for f, d in file_dates if start_dt <= d <= end_dt]
        if not in_range:
            raise FileNotFoundError(f"No Zarr files found between {start_date_str} and {end_date_str} in {self.D_FI}")
        return in_range

    def prepare_data_for_plotting(self, P_zarr=None, var_name=None):
        """
        Load and prepare Zarr data for plotting, handling both single and multi-date windows.

        Reads the Zarr file(s), extracts the plotting variable, and flattens the
        non-zero/non-NaN points into a pandas DataFrame. This prepares the data
        for visualisation with PyGMT scatter plotting.

        For time-averaged variables (e.g., multi-day averages), this method loops over
        the available Zarr files and computes the mean field across time. For single-date
        variables (e.g., FIC), only one Zarr file is read.

        Parameters
        ----------
        P_zarr   : Path, optional
                   Explicit path to the Zarr file (used for single-date cases).
        var_name : str, optional
                   Variable name to load. Defaults to `self.default_var`.

        Sets
        ----
        self.plot_df      : pandas.DataFrame
                            Columns are 'lat', 'lon', and 'dat' (flattened field).
        self.var_name     : str
                            Name of the loaded variable.
        self.plt_dt_str   : str
                            Date string for plot title or filename.
        """
        if var_name is None:
            var_name = self.default_var

        time_avg_vars = [f"{self.ice_type}P"] + [f"{self.ice_type}{s}_SD" for s in ["HI", "STH", "SH", "DIV", "AG", "AD", "AT", "VD", "VT", "STR"]]

        if self.dtN_str is not None and var_name in time_avg_vars:
            print(f"⏳ Aggregating spatial field over time for variable: {var_name}")
            sum_da = None
            count = 0
            for P_zarr in self.P_zarr_list:
                ds = xr.open_zarr(P_zarr)
                if var_name not in ds:
                    print(f"⚠️ Variable '{var_name}' not found in {P_zarr.name}. Skipping.")
                    continue
                da = ds[var_name].values
                if np.sum(da) > 0:
                    if count == 0:
                        sum_da = da
                    else:
                        sum_da += da
                    count += 1
            if count == 0:
                raise ValueError(f"❌ No valid data found for variable '{var_name}' over time range.")
            avg_da = sum_da / count
            lat = ds[self.config["CICE_dict"][f"{self.ice_type}_lat_coord"]].values
            lon = ds[self.config["CICE_dict"][f"{self.ice_type}_lon_coord"]].values
            dat = avg_da
            mask = (~np.isnan(dat)) & (dat != 0)
            self.var_name = var_name
            self.plot_df = pd.DataFrame({
                "lat": lat.ravel()[mask.ravel()],
                "lon": lon.ravel()[mask.ravel()],
                "dat": dat.ravel()[mask.ravel()]
            })
            dt0 = datetime.strptime(self.dt0_str, "%Y-%m-%d")
            dtN = datetime.strptime(self.dtN_str, "%Y-%m-%d")
            dt_mid = dt0 + (dtN - dt0) / 2
            self.plt_dt_str = dt_mid.strftime("%Y-%m-%d")
            return

        if P_zarr is None:
            P_zarr = self.P_zarr
        ds = xr.open_zarr(P_zarr)
        if var_name not in ds:
            raise KeyError(f"Variable '{var_name}' not found in {P_zarr}")
        da = ds[var_name].squeeze()
        lat = ds[self.config["CICE_dict"][f"{self.ice_type}_lat_coord"]].values
        lon = ds[self.config["CICE_dict"][f"{self.ice_type}_lon_coord"]].values
        dat = da.values
        mask = ~np.isnan(dat)
        self.var_name = var_name
        self.plot_df = pd.DataFrame({
            "lat": lat.ravel()[mask.ravel()],
            "lon": lon.ravel()[mask.ravel()],
            "dat": dat.ravel()[mask.ravel()]
        })

    def plot_path(self, region=None, sim_name=None, var_name=None, D_graph=None, hemisphere=None):
        """
        Generate figure output path and filename for the current plotting configuration.

        This method sets `self.P_png` to a resolved PNG file path that depends on
        simulation name, hemisphere, plotting region, and selected variable.

        Parameters
        ----------
        region     : str, optional
                     Region name, e.g., "Weddell", used for subdirectory.
        sim_name   : str, optional
                     Simulation name override. Defaults to `self.sim_name`.
        var_name   : str, optional
                     Variable name override. Defaults to `self.var_name`.
        D_graph    : Path, optional
                     Directory to write to. Defaults to internal plot path.
        hemisphere : str, optional
                     Hemisphere label for folder structure.

        Sets
        ----
        self.P_png : str
                     Full file path to the output PNG.
        """
        sim_name   = sim_name or self.sim_name
        var_name   = var_name or self.var_name
        D_graph    = Path(D_graph) if D_graph else self.D_graph
        hemisphere = hemisphere or self.hemisphere
        if region:
            D_png = D_graph / "regional" / region / var_name
            filename = f"{self.dt0_str}_to_{self.dtN_str}_{sim_name}_{region}_{var_name}.png" if self.dtN_str else f"{self.dt0_str}_{sim_name}_{region}_{var_name}.png"
        else:
            D_png = D_graph / sim_name / hemisphere / var_name
            filename = f"{self.dt0_str}_to_{self.dtN_str}_{sim_name}_{var_name}.png" if self.dtN_str else f"{self.dt0_str}_{sim_name}_{var_name}.png"
        D_png.mkdir(parents=True, exist_ok=True)
        self.P_png = str(D_png / filename)

    def plot_regions(self, reg_dict=None, include=None, exclude=None, water_color='white'):
        """
        Generate a set of regional maps based on AFIM regional dictionary.

        This method loops through region definitions from config and calls
        `plot_map()` with their geographic bounds.

        Parameters
        ----------
        reg_dict    : dict, optional
                      Dictionary of AFIM-defined regions (defaults to config).
        include     : list, optional
                      If provided, only plot these region names.
        exclude     : list, optional
                      If provided, skip these region names.
        water_color : str, optional
                      Water fill color for PyGMT coast (default: 'white').
        """
        reg_dict = reg_dict or self.config['AF_regions']
        for reg_name, reg_cfg in reg_dict.items():
            if include and reg_name not in include:
                continue
            if exclude and reg_name in exclude:
                continue
            print(f"Plotting region: {reg_name}")
            self.plot_map(
                projection  = f"S{reg_cfg['MC']}/-90/17.5C",
                region      = reg_cfg['plt_ext'],
                region_name = reg_name,
                sq_size     = 0.2,
                water_color = water_color)

    def plot_map(self,
                 projection        = None,
                 region            = None,
                 region_name       = None,
                 cmap              = None,
                 cmap_reverse      = None,
                 series            = None,
                 sq_size           = 0.03,
                 dt_str            = None,
                 cbar_pos          = "JBC+w10c/0.5c+mc+h",
                 cbar_str          = None,
                 units             = None,
                 ice_shelves       = False,
                 water_color       = 'black',
                 land_color        = 'seashell',
                 shoreline_str     = '.2p,white',
                 plain_frame       = False,
                 frame             = None,
                 text_loc          = None,
                 text_str          = None,
                 text_font         = "12p,Courier-Bold,black",
                 text_justify      = "LM",
                 text_fill         = "white",
                 text_pen          = "0.5p,black",
                 P_png             = None,
                 ow_fig            = None,
                 show_figs         = None):
        """
        Plot sea ice field using PyGMT for either regional or hemispheric view.

        This method uses `self.plot_df` to create a PyGMT scatter map.
        It supports configurable options for color map, scale, annotations,
        regional focus, and overlays such as grounded icebergs or shelves.

        Parameters
        ----------
        projection     : str, optional
                         GMT-style projection string.
        region         : list or str, optional
                         Geographic extent as [lon_min, lon_max, lat_min, lat_max].
        region_name    : str, optional
                         Label for regional subdirectory.
        cmap           : str, optional
                         Colormap to use.
        cmap_reverse   : bool, optional
                         Whether to reverse the colormap.
        series         : tuple, optional
                         Min/max data values for colorbar.
        sq_size        : float, optional
                         Size of PyGMT square symbol.
        dt_str         : str, optional
                         Date label override for the plot.
        cbar_pos       : str, optional
                         PyGMT colorbar position string.
        cbar_str       : str, optional
                         Colorbar label.
        units          : str, optional
                         Units string for axis.
        ice_shelves    : bool, optional
                         Whether to overlay Antarctic ice shelves.
        water_color    : str, optional
                         Color to fill water bodies.
        land_color     : str, optional
                         Color to fill land.
        shoreline_str  : str, optional
                         Style string for shorelines.
        plain_frame    : bool, optional
                         If True, removes axis labels and title.
        frame          : str or list, optional
                         Frame formatting.
        text_loc       : list or tuple, optional
                         Text label location (lon, lat).
        text_str       : str, optional
                         Label to overlay on map.
        text_font      : str, optional
                         Font specification for label.
        text_justify   : str, optional
                         Justification of text.
        text_fill      : str, optional
                         Fill behind text.
        text_pen       : str, optional
                         Outline around text.
        P_png          : str, optional
                         Override output file path.
        ow_fig         : bool, optional
                         Force overwrite for this plot.
        show_figs      : bool, optional
                         Force show for this plot.
        """
        meta         = self.plot_var_dict[self.var_name]
        cmap         = cmap         if cmap         is not None else meta.get("cmap", "viridis")
        series       = series       if series       is not None else meta.get("series", [0, 1])
        cmap_reverse = cmap_reverse if cmap_reverse is not None else meta.get("reverse", False)
        units        = units        if units        is not None else meta.get("units", "")
        cbar_str     = cbar_str     if cbar_str     is not None else meta.get("name", self.var_name)
        ow_fig       = ow_fig       if ow_fig       is not None else self.ow_fig
        show_figs    = show_figs    if show_figs    is not None else self.show_figs
        self.create_cbar_frame(series, cbar_str, units=units)
        if P_png is None:
            self.plot_path(region=region_name)
        if os.path.exists(self.P_png) and not ow_fig:
            print(f"figure already exists and not overwriting ... skipping {self.P_png}")
            return
        if (projection is None) or (region is None):
            print("*** PLOTTING HEMISPHERE ***")
            projection = self.hemisphere_projection
            region     = self.hemisphere_map_extent
            text_loc   = self.hemisphere_map_text_location
            print(f"\tprojection: {projection}")
            print(f"\tregion: {region}")
        if text_str is not None:
            print(f"\ttext location: {text_loc}")
        if plain_frame:
            frame = ["af"]
        elif frame is not None:
            frame = frame
        else:
            frame = ["af", f"+t{self.var_name} {self.plt_dt_str}"]
        t0  = time.time()
        fig = pygmt.Figure()
        pygmt.config(FONT_TITLE="20p,Courier-Bold", FONT_ANNOT_PRIMARY="14p,Helvetica")
        fig.basemap(projection=projection, region=region, frame=frame)
        fig.coast(land=land_color, water=water_color)
        pygmt.makecpt(cmap=cmap, reverse=cmap_reverse, series=series)
        if not self.plot_df.empty:
            fig.plot(x=self.plot_df.lon, y=self.plot_df.lat, cmap=True, fill=self.plot_df.dat, style=f"s{sq_size}c")
        else:
            print(f"**WARNING: NO DATA IN DATAFRAME**\n\tdataframe for timestep: {self.plt_dt_str} does not contain any data; writing out map without data")
        fig.colorbar(position=cbar_pos, frame=self.cbar_frame)
        if text_str is not None:
            fig.text(x=text_loc[0], y=text_loc[1], text=text_str,
                     font=text_font, justify=text_justify, fill=text_fill, pen=text_pen, no_clip=True)
        if ice_shelves:
            if not hasattr(self, 'antarctic_ice_shelves'):
                self.load_antarctic_ice_shelves()
            fig.plot(data=self.antarctic_ice_shelves, pen="0.2p,gray", fill="lightgray")
        if self.use_gi:
            fig.plot(x=self.gi_processor.GI_lon_cells, y=self.gi_processor.GI_lat_cells, fill="purple", style=f"c.1c")
            fig.plot(x=self.gi_processor.GI_lons, y=self.gi_processor.GI_lats, fill="black", style=f"s.1c")
        fig.coast(shorelines=shoreline_str)
        if self.show_figs:
            fig.show()
        else:
            print(f"\tsaving figure to: {self.P_png}")
            fig.savefig(self.P_png)
        print(f"\ttime taken {time.time()-t0} seconds")

    def plot_regions_grounded_icebergs(self, use_thinned=True, sq_size=0.15, cmap="viridis"):
        regions = self.config.get("AF_regions", {})
        for reg_name, reg_cfg in regions.items():
            print(f"Plotting region: {reg_name}")
            self.plot_map_grounded_icebergs(
                region=reg_cfg['plt_ext'],
                region_name=reg_name,
                sq_size=sq_size,
                cmap=cmap,
                projection=f"S{reg_cfg['MC']}/-90/17.5c",
                save_path=os.path.join(
                    self.config['D_dict']['graph'],
                    "grounded_icebergs",
                    "regional",
                    reg_name,
                    f"GI_{'thinned' if use_thinned else 'raw'}_{self.sim_name}.png"
                ),
                use_thinned=use_thinned
            )

    def plot_map_grounded_icebergs(self, region, region_name, projection, sq_size, cmap, save_path, use_thinned=True):
        """
        Create a PyGMT map of grounded iceberg data for a single Antarctic region.

        Parameters:
        -----------
        region : list
            GMT-style region bounding box [lon_min, lon_max, lat_min, lat_max].
        region_name : str
            Name of the region (used in filename or title).
        projection : str
            PyGMT polar projection string (e.g., 'S0/-90/17.5c').
        sq_size : float
            Size of square markers for grounded icebergs.
        cmap : str
            Colormap to use for iceberg count shading.
        save_path : str
            Full path to save output PNG.
        use_thinned : bool
            Whether to use thinned counts or raw grounded iceberg counts.
        """
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        fig = pygmt.Figure()
        fig.basemap(region=region, projection=projection, frame=["af", f"+tGrounded Icebergs ({'Thinned' if use_thinned else 'Raw'})"])
        fig.coast(land="lightgray", water="white", shorelines="0.5p,black")

        da = self.GI_thin_da if use_thinned else self.GI_cnts
        lon = self.G_t['lon'].values
        lat = self.G_t['lat'].values
        dat = da.values
        mask = (dat > 0) & ~np.isnan(dat)

        pygmt.makecpt(cmap=cmap, series=[np.min(dat[mask]), np.max(dat[mask])])

        fig.plot(
            x=lon[mask],
            y=lat[mask],
            style=f"s{sq_size}c",
            fill=dat[mask],
            cmap=True,
            pen="thin"
        )

        fig.colorbar(position="JBC+w10c/0.5c+h", frame=["x+lGrounded Iceberg Count"])

        # Add ice shelf overlay if loaded
        if self.ice_shelves_loaded:
            fig.plot(data=self.antarctic_ice_shelves, pen="0.2p,gray", fill="gainsboro")

        # Regional stats
        count_vals = dat[mask]
        stats = [
            f"Total Cells: {count_vals.size}",
            f"Mean Count: {np.mean(count_vals):.2f}",
            f"Max Count: {np.max(count_vals)}"
        ]
        if use_thinned and hasattr(self, 'GI_thin_fact'):
            stats.append(f"Thinning: {self.GI_thin_fact:.2%}")

        fig.text(
            x=region[0] + 1,
            y=region[2] + 1,
            text="\n".join(stats),
            font="10p,Helvetica,black",
            justify="LM",
            fill="white",
            pen="0.5p,black"
        )

        fig.savefig(save_path)
        print(f"Saved: {save_path}")

    def plot_sectors_grounded_icebergs(self, kmt_orig_path, kmt_mod_path, n_sectors=9, sector_width=40, lat_limit=-60, start_lon=0, marker_size="0.1c", save=False, save_dir=None, overlay_ice_shelves=False):
        """
        Plot sector-based PyGMT maps of newly inserted grounded iceberg cells based on changes in KMT.

        Parameters:
        -----------
        kmt_orig_path : str
            Path to original KMT file before grounded iceberg masking.
        kmt_mod_path : str
            Path to modified KMT file after grounded iceberg masking.
        n_sectors : int
            Number of longitude sectors to divide the Southern Ocean.
        sector_width : int
            Width of each longitude sector in degrees.
        lat_limit : float
            Latitude cutoff to restrict plotted points (typically -60 for Antarctic shelf).
        start_lon : float
            Starting longitude for sector division (0 for Greenwich).
        marker_size : str
            Marker size string for PyGMT (e.g., "0.1c").
        """
        ds_orig = xr.open_dataset(kmt_orig_path)
        ds_mod = xr.open_dataset(kmt_mod_path)

        kmt_orig_arr = ds_orig['kmt'].values
        kmt_mod_arr = ds_mod['kmt'].values
        lat_arr = ds_mod['lat'].values
        lon_arr = ds_mod['lon'].values
        lon_arr_wrapped = (lon_arr + 360) % 360

        diff_mask = (kmt_orig_arr == 1) & (kmt_mod_arr == 0)
        print("New grounded iceberg cells:", np.sum(diff_mask))

        nj, ni = np.where(diff_mask)
        lon_points = lon_arr_wrapped[nj, ni]
        lat_points = lat_arr[nj, ni]

        for i in range(n_sectors):
            fig = pygmt.Figure()
            lon_min = (start_lon + i * sector_width) % 360
            lon_max = (lon_min + sector_width) % 360
            if lon_min < lon_max:
                meridian_center = (lon_min + lon_max) / 2
            else:
                meridian_center = ((lon_min + lon_max + 360) / 2) % 360

            if lon_min < lon_max:
                sector_mask = (lon_points >= lon_min) & (lon_points < lon_max) & (lat_points <= lat_limit)
            else:
                sector_mask = ((lon_points >= lon_min) | (lon_points < lon_max)) & (lat_points <= lat_limit)

            sector_lons = lon_points[sector_mask]
            sector_lats = lat_points[sector_mask]

            if len(sector_lats) > 0:
                lat_min = max(-90, np.min(sector_lats) - 1.0)
                lat_max = min(-60, np.max(sector_lats) + 1.5)
            else:
                lat_min = -90
                lat_max = lat_limit
            if lon_max == 0:
                lon_max = 359
            region = [lon_min, lon_max, lat_min, lat_max]

            projection = f"S{meridian_center}/-90/17.5c"
            fig.basemap(region=region, projection=projection, frame=["af", f"+tSector {i+1}: {int(lon_min)} to {int(lon_max)}°"])
            fig.coast(shorelines=True, resolution="i", land="lightgray", water="white")

            if len(sector_lons) > 0:
                fig.plot(x=sector_lons, y=sector_lats, style=f"c{marker_size}", fill="red", pen="black")

                        # Annotate stats
            stat_text = [
                f"# GI cells: {len(sector_lons)}"
            ]
            if len(sector_lats) > 0:
                stat_text.append(f"Lat range: {np.min(sector_lats):.1f}° to {np.max(sector_lats):.1f}°")
                stat_text.append(f"Lon range: {int(lon_min)}° to {int(lon_max)}°")
            fig.text(
                x=region[0] + 1,
                y=region[2] + 1,
                text="\n".join(stat_text),
                font="10p,Helvetica,black",
                justify="LM",
                fill="white",
                pen="0.5p,black"
            )

            # Optional overlay of ice shelves
            if overlay_ice_shelves and self.ice_shelves_loaded:
                fig.plot(data=self.antarctic_ice_shelves, pen="0.2p,gray", fill="gainsboro")

            if save:
                if save_dir is None:
                    save_dir = os.path.join("grounded_iceberg_sectors", self.sim_name)
                os.makedirs(save_dir, exist_ok=True)
                save_path = os.path.join(save_dir, f"sector_{i+1:02d}_grounded_icebergs.png")
                fig.savefig(save_path, dpi=300)
                print(f"Saved: {save_path}")
            else:
                fig.show()
            # Optional: fig.savefig(f"sector_{i+1:02d}_grounded_icebergs.png", dpi=300)

    def plot_timeseries(self, var_names=None, label=None, save_png=True):
        """
        Plot fast ice time series from a multi-day Zarr dataset.

        This method extracts one or more 1D variables (e.g., FIA, FIV, etc.)
        from the Zarr time series dataset and plots them over time.

        Parameters
        ----------
        var_names : list of str, optional
                     List of variable names to plot. Defaults to ['FIA'].
        label     : str, optional
                     Legend label for this simulation (used if comparing multiple).
        save_png  : bool, optional
                     If True, save the figure as a PNG. If False, display interactively.

        Notes
        -----
        Requires that self.P_zarr_list is defined (i.e., for multi-date mode).
        """
        if self.dtN_str is None or not hasattr(self, 'P_zarr_list'):
            raise ValueError("Time series plotting requires a defined date range (dtN_str).")
        if var_names is None:
            var_names = [f"{self.ice_type}A"]
        all_time = []
        data_dict = {var: [] for var in var_names}
        for P_zarr in self.P_zarr_list:
            ds = xr.open_zarr(P_zarr)
            if "time" not in ds:
                raise KeyError(f"'time' coordinate not found in {P_zarr}")
            times = pd.to_datetime(ds["time"].values)
            all_time.extend(times)
            for var in var_names:
                if var in ds:
                    values = ds[var].values
                    if values.ndim != 1 or len(values) != len(times):
                        print(f"⚠️ Shape mismatch for '{var}' in {P_zarr.name}: {values.shape}")
                    data_dict[var].extend(values)
                else:
                    print(f"⚠️ Variable '{var}' missing in {P_zarr.name}, filling with NaN.")
                    data_dict[var].extend([np.nan] * len(times))
        df = pd.DataFrame(data=data_dict, index=pd.to_datetime(all_time))
        df = df.sort_index()
        plt.figure(figsize=(10, 6))
        for var in var_names:
            plot_label = f"{label} ({var})" if label and len(var_names) > 1 else (label or var)
            plt.plot(df.index, df[var], label=plot_label)
        plt.xlabel("Date")
        plt.ylabel("Value")
        title_vars = ", ".join(var_names)
        plt.title(f"{self.sim_name} {self.ice_type} Time Series ({title_vars})")
        plt.grid(True, linestyle="--", alpha=0.5)
        plt.legend()
        if self.show_figs and not save_png:
            plt.show()
        else:
            var_str = "_".join(var_names)
            F_save = f"{self.sim_name}_{var_str}_{self.dt0_str}_to_{self.dtN_str}.png"
            P_save = Path(self.D_ts, F_save)
            plt.savefig(P_save)
            print(f"✅ Saved comparison plot to: {P_save}")

    def plot_timeseries_compare(self, sim_names, var_name=None, label_dict=None, save_png=True,
                                comparison_name='baseline-comparison'):
        """
        Plot the specified fast ice variable from multiple simulations for comparison.

        Parameters
        ----------
        sim_names      : list of str
                         List of simulation names to compare.
        var_name       : str
                         Name of the variable to plot (default is "FIA").
        label_dict     : dict, optional
                         Dictionary mapping sim_name to label to use in plot legend.
        save_png       : bool
                         Whether to save the figure as PNG (default True).
        comparison_name: str
                         Identifier for the saved filename.

        Notes
        -----
        This method instantiates temporary SeaIcePlotter instances
        for each simulation and collects their time series data.
        """
        if self.dtN_str is None:
            raise ValueError("Comparison requires a defined dtN_str at SeaIcePlotter init.")
        if not isinstance(sim_names, list):
            raise TypeError("sim_names must be a list of simulation names.")
        var_name = var_name or f"{self.ice_type}A"
        plt.figure(figsize=(10, 6))
        for sim in sim_names:
            other = SeaIcePlotter(sim, self.dt0_str, self.dtN_str, ice_type=self.ice_type)
            if not hasattr(other, 'P_zarr_list'):
                raise AttributeError(f"No Zarr paths loaded for simulation: {sim}")
            time_list = []
            value_list = []
            for P_zarr in other.P_zarr_list:
                ds = xr.open_zarr(P_zarr)
                if "time" not in ds or var_name not in ds:
                    print(f"⚠️ Missing 'time' or '{var_name}' in {P_zarr}")
                    continue
                times = pd.to_datetime(ds["time"].values)
                values = ds[var_name].values
                if values.ndim != 1 or len(values) != len(times):
                    print(f"⚠️ Shape mismatch in {P_zarr.name}: {values.shape}")
                    continue
                time_list.extend(times)
                value_list.extend(values)
            if len(time_list) == 0:
                print(f"⚠️ No valid data loaded for simulation: {sim}")
                continue
            df = pd.DataFrame({var_name: value_list}, index=pd.to_datetime(time_list))
            df = df.sort_index()
            label = label_dict.get(sim, sim) if label_dict else sim
            plt.plot(df.index, df[var_name], label=label)
        plt.xlabel("Date")
        plt.ylabel(var_name)
        plt.title(f"{self.ice_type} Time Series Comparison: {var_name}")
        plt.grid(True, linestyle="--", alpha=0.5)
        plt.legend()
        if self.show_figs and not save_png:
            plt.show()
        else:
            F_save = f"{comparison_name}_{var_name}_{self.dt0_str}_to_{self.dtN_str}.png"
            P_save = Path(self.D_ts, F_save)
            plt.savefig(P_save)
            print(f"✅ Saved comparison plot to: {P_save}")

    def create_animation(self, var_name=None, region=None, fps=4, codec="libx264"):
        """
        Create an MP4 animation from saved PNG frames.

        Parameters
        ----------
        region_name : str, optional
                      Region label used in filename paths.
        var_name    : str, optional
                      Variable name to animate. Defaults to self.var_name.
        fps         : int
                      Frames per second for the output video.
        suffix      : str
                      File extension for input PNGs (default: .png).

        Notes
        -----
        This method uses `ffmpeg` under the hood. It must be installed and accessible
        in the system path.
        Output file is written to self.D_graph with a standard MP4 filename.
        """
        import imageio.vffmpeg as vffmpeg
        import imageio

        if self.dtN_str is None:
            raise ValueError("Animation can only be created when dtN_str is defined (i.e. for a date range).")

        var_name = var_name or self.default_var
        png_paths = []

        for P_zarr, dt_str in zip(self.P_zarr_list, self.plt_dt_str_list):
            self.prepare_data_for_plotting(P_zarr=P_zarr, var_name=var_name)
            self.plot_map(region_name=region)
            png_paths.append(self.P_png)

        output_path = Path(self.D_graph, f"animation_{self.sim_name}_{var_name}_{self.dt0_str}_to_{self.dtN_str}.mp4")
        with imageio.get_writer(output_path, fps=fps, codec=codec) as writer:
            for png in png_paths:
                img = imageio.imread(png)
                writer.append_data(img)
        print(f"🎞️ MP4 animation saved to: {output_path}")

    def plot_facet_regions(self, var_name=None, region_names=None, cols=3, figsize=(15, 10)):
        """
        Create a multi-panel (faceted) plot of multiple regions.

        Parameters
        ----------
        regions     : list of str
                      List of region names from config['AF_regions'].
        layout      : str
                      Facet layout, e.g., "3x2".
        label_dict  : dict, optional
                      Custom label for each region (title).

        Notes
        -----
        This method generates a grid of subplots using PyGMT subplot
        functionality. It relies on region extents defined in config['AF_regions'].
        """        
        import matplotlib.pyplot as plt
        import matplotlib.cm as cm

        var_name = var_name or self.default_var
        region_names = region_names or list(self.config['AF_regions'].keys())
        n_regions = len(region_names)
        rows = int(np.ceil(n_regions / cols))

        fig, axes = plt.subplots(rows, cols, figsize=figsize, constrained_layout=True)
        axes = axes.flatten()

        self.prepare_data_for_plotting(var_name=var_name)

        for idx, region in enumerate(region_names):
            reg_cfg = self.config['AF_regions'][region]
            ax = axes[idx]

            self.plot_map(
                projection  = f"S{reg_cfg['MC']}/-90/17.5C",
                region      = reg_cfg['plt_ext'],
                region_name = region,
                show_figs   = False,
                ow_fig      = True,
                sq_size     = 0.15,
                P_png       = None
            )
            img = plt.imread(self.P_png)
            ax.imshow(img)
            ax.set_title(region)
            ax.axis('off')

        for j in range(n_regions, rows * cols):
            axes[j].axis('off')

        plt.suptitle(f"{self.sim_name} {var_name} Regional Overview ({self.dt0_str})", fontsize=16)
        out_path = Path(self.D_graph, f"facet_{self.sim_name}_{var_name}_{self.dt0_str}.png")
        fig.savefig(out_path)
        print(f"🗺️ Faceted regional plot saved to: {out_path}")

