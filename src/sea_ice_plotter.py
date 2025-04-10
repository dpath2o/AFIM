import os, sys, time, json, pygmt
import xarray as xr
import pandas as pd
import numpy as np
import geopandas as gpd
from pathlib import Path
from datetime import datetime
import matplotlib.pyplot as plt
sys.path.insert(0, '/home/581/da1339/AFIM/src/AFIM/src')
from grounded_iceberg_processor import GroundedIcebergProcessor

class SeaIcePlotter:
    """
    A plotting class for visualizing sea ice model output (e.g., from CICE or AFIM) in various formats
    including regional maps, hemisphere-wide views, and time series plots.

    The class supports output from different ice types (fast ice, pack ice, sea ice, grounded icebergs)
    and integrates configuration from a JSON settings file to control plotting behavior, file paths,
    and variable metadata.

    Key features:
    - Supports plotting from in-memory xarray datasets or from Zarr files
    - Integrates with grounded iceberg data if available
    - Regional and hemisphere plotting via PyGMT
    - Loads plotting styles from `pygmt_dict` and variable metadata from `plot_var_dict`
    - Fully configurable through JSON to adapt to different simulations and data products
    """
    def __init__(self,
                 P_JSON     = None,
                 sim_name   = None,
                 ice_type   = None, #FI, SI, PI, GI                   -> controls dictionary and for GI controls sub-directory
                 plot_type  = None, #regional, hemisphere, timeseries -> controls PNG directory
                 var_name   = None,
                 dt0_str    = None,
                 dtN_str    = None,
                 overwrite  = None,
                 save_fig   = None,
                 show_fig   = None,
                 hemisphere = None):
        """
        Initialize a SeaIcePlotter instance using configuration provided in a JSON file.

        Parameters
        ----------
        P_JSON : str or Path, optional
            Path to the configuration JSON file. If None, defaults to a pre-defined internal path.
        sim_name : str, optional
            Name of the simulation. Controls file paths and Zarr subdirectories.
        ice_type : str, optional
            Type of ice data ('FI', 'PI', 'SI', 'GI'). Controls dictionary and data sub-paths.
        plot_type : str, optional
            Type of plot to generate: 'regional', 'hemisphere', or 'timeseries'.
        var_name : str, optional
            Name of the variable to plot. Must exist in `plot_var_dict` in the config JSON.
        dt0_str : str, optional
            Start date (inclusive) for plots, in 'YYYY-MM-DD' format. Defaults to value in config.
        dtN_str : str, optional
            End date (inclusive) for plots, in 'YYYY-MM-DD' format. Defaults to value in config.
        overwrite : bool, optional
            If True, overwrite existing figures. Defaults to False.
        save_fig : bool, optional
            If True, save the generated figures to PNG. Defaults to False.
        show_fig : bool, optional
            If True, display figures after plotting. Defaults to True.
        hemisphere : str, optional
            Geographic hemisphere ('north' or 'south'). Determines projection and map extent.

        Notes
        -----
        - Loads multiple dictionaries from JSON, including:
            * D_dict       → file paths
            * CICE_dict    → coordinate names
            * pygmt_dict   → default PyGMT styling
            * plot_var_dict→ plotting metadata per variable
            * AF_regions   → predefined regional extents
        - Initializes a GroundedIcebergProcessor to overlay grounded iceberg locations if enabled.
        """
        if P_JSON is None:
            P_JSON = "/home/581/da1339/AFIM/src/AFIM/src/JSONs/afim_cice_analysis.json"
        with open(P_JSON, 'r') as f:
            self.config = json.load(f)
        self.sim_name    = sim_name         if sim_name   is not None else 'baseline'
        self.ice_type    = ice_type.upper() if ice_type   is not None else 'FI'
        self.plot_type   = plot_type        if plot_type  is not None else 'regional'
        self.var_name    = var_name         if var_name   is not None else 'FIC'
        self.dt0_str     = dt0_str          if dt0_str    is not None else self.config.get('dt0_str', '1993-01-01')
        self.dtN_str     = dtN_str          if dtN_str    is not None else self.config.get('dtN_str', '1999-12-31')
        self.ow_fig      = overwrite        if overwrite  is not None else False
        self.save_fig    = save_fig         if save_fig   is not None else False
        self.show_fig    = show_fig         if show_fig   is not None else True
        self.hemisphere  = hemisphere       if hemisphere is not None else self.config.get('hemisphere', 'south')
        self.define_hemisphere(self.hemisphere)
        self.ice_type_dict             = self.config[f"{self.ice_type}_var_dict"] if self.ice_type!='GI' else None
        self.spatial_distribution_vars = [f"{self.ice_type}P"] + [f"{self.ice_type}{s}_SD" for s in ["HI", "STH", "SH", "DIV", "AG", "AD", "AT", "VD", "VT", "STR"]]
        self.D_zarr                    = Path(self.config['D_dict']['AFIM_out'], sim_name, self.ice_type)
        self.D_graph                   = Path(self.config['D_dict']['graph'], 'AFIM')
        self.CICE_dict                 = self.config.get("CICE_dict", {})
        self.pygmt_dict                = self.config.get("pygmt_dict", {})
        self.plot_var_dict             = self.config.get("plot_var_dict", {})
        self.reg_dict                  = self.config.get('AF_regions', {})
        self.sim_dict                  = self.config.get("sim_dict", {})
        self.GI_dict                   = self.config.get("GI_dict", {})
        self.SIC_scale                 = self.config.get("SIC_scale", 1e12)
        self.gi_processor              = GroundedIcebergProcessor(sim_name=sim_name)
        self.gi_processor.load_AFIM_GI()
        self.use_gi        = self.gi_processor.use_gi
        self.GI_total_area = self.gi_processor.total_area
        self.KMT_mod       = self.gi_processor.KMT_mod
        # if dtN_str is None:
        #     P_zarr, self.plt_dt_str = self._find_nearest_zarr_file(start_date_str=dt0_str)
        #     print(f"Using nearest Zarr file: {P_zarr.name}")
        #     self.P_zarr  = P_zarr
        #     self.dt0_str = self.plt_dt_str
        #     self.dtN_str = None
        # else:
        #     zarr_pairs = self._find_nearest_zarr_file(start_date_str=dt0_str, end_date_str=dtN_str)
        #     if zarr_pairs is not None:
        #         self.P_zarr_list     = [p for p, _ in zarr_pairs]
        #         self.plt_dt_str_list = [d for _, d in zarr_pairs]
        #         self.dt0_str = self.plt_dt_str_list[0]
        #         self.dtN_str = self.plt_dt_str_list[-1]
        #     else:
        #         self.P_zarr_list     = []
        #         self.plt_dt_str_list = []
        #         self.dt0_str         = dt0_str
        #         self.dtN_str         = dtN_str if dtN_str is not None else self.config.get("dtN_str","1999-12-31")

    def define_hemisphere(self, hemisphere):
        """
        THIS FUNCTION SHOULD BE *EXACTLY* THE SAME AS THE ONE IN SEAICEPROCESSOR CLASS
        """
        if hemisphere.lower() in ['north', 'northern', 'nh', 'n', 'no']:
            self.hemisphere_geographic_extent = [0, 360, 0, 90]
            self.hemisphere_map_extent        = [-180,180,55,90]
            self.hemisphere_projection        = "S0.0/90.0/50/15C"
            self.hemisphere_map_text_location = [-120,56]
            self.hemisphere_abbreviation      = 'NH'
            self.hemisphere_nj_slice          = slice(540,1080)
            self.hemisphere                   = 'north'
        elif hemisphere.lower() in ['south', 'southern', 'sh', 's', 'so']:
            self.hemisphere_geographic_extent = [0, 360, -90, 0]
            self.hemisphere_map_extent        = [-180,180,-90,-55]
            self.hemisphere_projection        = "S0.0/-90.0/50/15C"
            self.hemisphere_map_text_location = [0,-90]
            self.hemisphere_abbreviation      = 'SH'
            self.hemisphere_nj_slice          = slice(0,540)
            self.hemisphere                   = 'south'
        else:
            raise ValueError(f"Invalid hemisphere '{hemisphere}'. Valid options are: "
                             "['north', 'south', 'northern', 'southern', 'sh', 'nh', 'SH', 'NH']")

    def load_ice_shelves(self):
        """
        Load Antarctic ice shelf polygons from a shapefile into the processor.

        This sets the `self.antarctic_ice_shelves` attribute for optional overlay
        during plotting. Geometries are cleaned, reprojected, and buffered.
        """
        coast_path              = self.config['pygmt_dict']['P_coast_shape']
        gdf                     = gpd.read_file(coast_path)
        shelves                 = gdf[gdf['POLY_TYPE'] == 'S']
        shelves                 = shelves[~shelves.geometry.is_empty & shelves.geometry.notnull()]
        shelves                 = shelves.to_crs("EPSG:4326")
        shelves.geometry        = shelves.geometry.buffer(0)
        self.ice_shelves        = shelves
        self.ice_shelves_loaded = True

    def _find_nearest_zarr_file(self, start_date_str=None, end_date_str=None):
        start_date_str = start_date_str if start_date_str is not None else self.dt0_str
        end_date_str   = end_date_str if end_date_str is not None else self.dtN_str
        if self.ice_type=='FI':
            F_prefix = "fast_ice"
        elif self.ice_type=='PI':
            F_prefix = "pack_ice"
        elif self.ice_type=='SI':
            F_prefix = "sea_ice"
        def extract_date(f):
            return datetime.strptime(f.stem.replace(f"{F_prefix}_", ""), "%Y-%m-%d")
        all_files = sorted(self.D_zarr.glob(f"{F_prefix}_*.zarr"))
        if not all_files:
            print(f"No Zarr files found in {self.D_zarr}")
            print("returning empty list")
            return
        file_dates = [(f, extract_date(f)) for f in all_files]
        if end_date_str is None:
            target_date = datetime.strptime(start_date_str, "%Y-%m-%d")
            closest_file, matched_date = min(file_dates, key=lambda x: abs(x[1] - target_date))
            return closest_file, matched_date.strftime("%Y-%m-%d")
        start_dt = datetime.strptime(start_date_str, "%Y-%m-%d")
        end_dt   = datetime.strptime(end_date_str, "%Y-%m-%d")
        in_range = [(f, d.strftime("%Y-%m-%d")) for f, d in file_dates if start_dt <= d <= end_dt]
        if not in_range:
            print(f"No Zarr files found between {start_date_str} and {end_date_str} in {self.D_zarr}")
            print("returning empty list")
            return
        else:
            return in_range

    def create_cbar_frame(self, series, label, units=None, max_ann_steps=10):
        """
        Given a data range and label, create GMT-style cbar string.
        """
        vmin, vmax = series
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
        ann_str = f"{ann_step:.3f}".rstrip("0").rstrip(".")
        tick_str = f"{tick_step:.3f}".rstrip("0").rstrip(".")
        if units is not None:
            self.cbar_frame = [f"a{ann_str}f{tick_str}+l{label}",f"y+l {units}"]
        else:
            self.cbar_frame = f"a{ann_str}f{tick_str}+l{label}"

    def define_png_plot_path(self, qualifiers=None):
        """
        """
        if self.plot_type=='regional':
            D = Path(self.D_graph,self.sim_name,self.ice_type,self.plot_type,self.region_name,self.var_name)
            if self.ice_type!='GI':
                F = f"{self.dt0_str}_{self.sim_name}_{self.region_name}_{self.var_name}"
            else:
                F = f"GI_{self.sim_name}_{self.region_name}"
        elif self.plot_type=='timeseries':
            D = Path(self.D_graph,self.sim_name,self.ice_type,self.plot_type,self.var_name)
            F = f"{self.dt0_str}_{self.dtN_str}_{self.sim_name}_{self.var_name}"
        elif self.plot_type=='hemisphere':
            D = Path(self.D_graph,self.sim_name,self.ice_type,self.hemisphere,self.var_name)
            F = f"{self.dt0_str}_{self.sim_name}_{self.hemisphere}_{self.var_name}"
        if qualifiers is not None:
            for q in qualifiers:
                F += f"{F}_{q_}"
        F = f"{F}.png" #add suffix/filetype
        D.mkdir(parents=True, exist_ok=True)
        self.P_png = str(D/F)
        print(f"plot path: {self.P_png}")

    def get_meridian_center_from_geographic_extent(self, geographic_extent):
        """
        Determine the meridian center for PyGMT stereographic projection
        based on a region list [min_lon, max_lon, min_lat, max_lat].

        Works for longitudes in [-180, 180] or [0, 360], and handles dateline wrap.

        Returns:
        float: Central meridian in [-180, 180] for use in 'S<lon>/lat/width' projections.
        """
        lon_min, lon_max = geographic_extent[0], geographic_extent[1]
        lon_min_360 = lon_min % 360
        lon_max_360 = lon_max % 360
        if (lon_max_360 - lon_min_360) % 360 > 180:
            # Geographic_Extent crosses dateline
            center = ((lon_min_360 + lon_max_360 + 360) / 2) % 360
        else:
            center = (lon_min_360 + lon_max_360) / 2
        # Convert back to [-180, 180]
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

    def dataarray_2D_to_dataframe(self, da=None, lat=None, lon=None):
        """
        np.shape(da) == np.shape(lat) == np.shape(lon)
        """
        mask         = ~np.isnan(da)
        self.plot_df = pd.DataFrame({"lat": lat.ravel()[mask.ravel()],
                                     "lon": lon.ravel()[mask.ravel()],
                                     "dat": da.ravel()[mask.ravel()]})

    def aggergate_2D_spatial_distribution(self, ds=None, P_zarr_list=None):
        """
        if ds is none then search for zarr files with _find_nearest_zarr_file(self, start_date_str=None, end_date_str=None), where if dt0_str is None or dtN_str is None then use self.dt0_str and self.dtN_str
        """
        if self.var_name is None or self.var_name not in self.spatial_distribution_vars:
            raise ValueError(f"❌ method requires '{self.var_name}' to be defined\nand for it to be a spatial *known* spatial distribution variable: {self.spatial_distribution_vars:}")
        print(f"⏳ Aggregating spatial field over time for variable: {self.var_name}")
        if ds is not None:
            da = ds[self.var_name].sum(dim=self.time_coord_name).values / len(ds[self.time_coord_name].values)
        else: 
            if P_zarr_list is None:
                zarr_pairs = self._find_nearest_zarr_file(start_date_str=self.dt0_str, end_date_str=self.dtN_str)
                if zarr_pairs is not None:
                    P_zarr_list = [p for p, _ in zarr_pairs]
            ds = xr.open_mfdataset(P_zarr_list, engine="zarr", consolidated=True)
            da = ds[self.var_name].sum(dim=self.time_coord_name).values / len(ds[self.time_coord_name].values)
        t_mid_index        = int(len(ds[self.time_coord_name].values)/2)
        self.plot_date_str = pd.to_datetime(ds[self.time_coord_name].isel({self.time_coord_name : t_mid_index}).values[0]).strftime("%Y-%m-%d") 
        self.frame         = ["af", f"+t{self.var_name} mid-date {self.plot_date_str} for {int(t_mid_index*self.config.get('roll_win'))}-days either side"]
        lat                = ds[self.lat_coord_name].values
        lon                = ds[self.lon_coord_name].values
        mask               = da>0
        lat_mask           = lat[mask]
        lon_mask           = lon[mask]
        da_mask            = da[mask]
        self.dataarray_2D_to_dataframe( da=da_mask, lat=lat_mask, lon=lon_mask)

    def regional_or_hemisphere(self):
        if self.plot_type=='hemisphere':
            print("*** PLOTTING HEMISPHERE ***")
            self.projection    = self.hemisphere_projection
            self.region_extent = self.hemisphere_map_extent
            self.text_loc      = self.hemisphere_map_text_location
            self.create_figure()
        elif self.plot_type=='regional':
            print("*** REGIONAL PLOTTING ***")
            if self.region_name is not None:
                if self.region_name in self.region_dictionary.keys():
                    self.region_extent = self.region_dictionary[self.region_name]['plt_ext']
                    self.get_meridian_center_from_geographic_extent(self.region_extent)
                    self.projection = f"S{self.plot_meridian_center}/-90/17.5C"
                    self.create_figure()
                else:
                    raise ValueError(f"❌ {self.region_name} is not in {self.region_dictionary.keys()}")
            else:
                for self.region_name, region_config in self.region_dictionary.items():
                    self.region_extent = region_config['plt_ext']
                    self.get_meridian_center_from_geographic_extent(self.region_extent)
                    self.projection = f"S{self.plot_meridian_center}/-90/17.5C"
                    self.create_figure()
                self.region_name = None

    def classify_strain(self,val,min_strain=-2,max_strain=2):
        if val < min_strain:
            return 0  # converging
        elif val > max_strain:
            return 2  # diverging
        else:
            return 1  # no strain rate

    def create_figure(self):
        if self.P_png is None:
            self.define_png_plot_path(qualifiers = self.png_name_extras)
        if os.path.exists(self.P_png) and not self.ow_fig:
            print(f"figure already exists and not overwriting ... skipping {self.P_png}")
            return
        if self.units is not None and self.cmap!="categorical":
            self.create_cbar_frame(self.series, self.cbar_str, units=self.units)
        if self.text_str is not None:
            print(f"\ttext location: {self.text_loc}")
        if self.plain_frame:
            frame = ["af"]
        elif self.frame is not None:
            frame = self.frame
        else:
            frame = ["af", f"+t{self.var_name} {self.plot_date_str}"]
        t0  = time.time()
        fig = pygmt.Figure()
        pygmt.config(FONT_TITLE="20p,Courier-Bold", FONT_ANNOT_PRIMARY="14p,Helvetica")
        print(f"plotting figure with projection {self.projection} and frame {frame}")
        fig.basemap(projection=self.projection, region=self.region_extent, frame=frame)
        fig.coast(land=self.land_color, water=self.water_color)
        if not self.plot_df.empty:
            if self.cmap=="categorical":
                cat_labs = self.plot_var_dict[self.var_name]["cat_labs"]
                if self.var_name in ["FIDIV","SIDIV","PIDIV","FIDIV_SD","SIDIV_SD","PIDIV_SD","divu"]:
                    self.plot_df["cat_index"] = self.plot_df["dat"].apply(lambda x: self.classify_strain(x, min_strain=self.series[0], max_strain=self.series[1]))
                pygmt.makecpt(cmap        = "gmt/categorical",
                              series      = [0, len(cat_labs)-1, 1],
                              color_model = "+c" + ",".join(cat_labs) )
                fig.plot(x=self.plot_df.lon, y=self.plot_df.lat, fill=self.plot_df.cat_index, cmap=True, style=f"s{self.sq_size_var}c" )
                fig.colorbar()
            else:
                pygmt.makecpt(cmap=self.cmap, reverse=self.cmap_reverse, series=self.series)
                fig.plot(x=self.plot_df.lon, y=self.plot_df.lat, cmap=True, fill=self.plot_df.dat, style=f"s{self.sq_size_var}c")
                fig.colorbar(position=self.cbar_pos, frame=self.cbar_frame)
        else: # assume just plotting latitudes and longitudes of something with basic fill color
            print(f"**WARNING: NO DATA IN DATAFRAME**\n\ttimestep: {self.plot_date_str}")
            #fig.plot(x=self.plot_df.lon, y=self.plot_df.lat, fill=self.fill_color, style=f"s{self.sq_size_var}c")
        if self.text_str is not None:
            fig.text(x=self.text_loc[0], y=self.text_loc[1], text=self.text_str, font=self.text_font, justify=self.text_justify, fill=self.text_fill, pen=self.text_pen, no_clip=True)
        if self.plot_ice_shelves:
            if not hasattr(self, 'ice_shelves'):
                self.load_ice_shelves()
            fig.plot(data=self.ice_shelves, pen="0.2p,gray", fill="lightgray")
        if self.use_gi and self.var_name is not None:
            fig.plot(x=self.gi_processor.GI_lon_cells, y=self.gi_processor.GI_lat_cells, fill=self.fill_color, style=f"c{self.sq_size_GI}c")
        fig.coast(shorelines=self.shoreline_str)
        if self.show_fig:
            fig.show()
        if self.save_fig:
            fig.savefig(self.P_png)
            print(f"figure saved: {self.P_png}")
        print(f"time taken {time.time()-t0:0.2f} seconds\n")

    def plot_map(self,
                 ds            = None,
                 P_zarr        = None,
                 P_zarr_list   = None,
                 var_name      = None,
                 dt0_str       = None,
                 dtN_str       = None,
                 time_index    = None,
                 single_figure = False,
                 region_name   = None,
                 plain_frame   = False,
                 frame         = None,
                 text_str      = None,
                 P_png         = None,
                 save_fig      = None,
                 ow_fig        = None,
                 show_fig      = None,
                 **kwargs):
        """
        This method plots 2D gridded data using PyGMT.

        ---- INPUT OPTIONS ----
        There are four ways data can be passed:
        1. ds (xarray.Dataset)      → Directly pass a dataset (e.g., from SeaIceProcessor or raw CICE).
        2. P_zarr (Path)            → A single Zarr file to plot a single timestep.
        3. P_zarr_list (list[Path])→ Multiple Zarr files to loop over and plot multiple figures.
        4. None                     → Method will auto-detect files using `dt0_str` and `dtN_str` date bounds.

        For options (1) and (4), time-based loops can be triggered using `dt0_str` and `dtN_str`.

        ---- COORDINATE METADATA ----
        You can override coordinate names using:
        - lon_coord_name (default from config["CICE_dict"]["proc_lon_coord"])
        - lat_coord_name (default from config["CICE_dict"]["proc_lat_coord"])
        - time_coord_name (default from config["CICE_dict"]["proc_time_dim"])

        ---- REGION CONTROLS ----
        `plot_type` determines whether to generate regional or hemisphere-wide plots.
        For regional plotting:
        - `region_dictionary` defines regions (defaults to config['AF_regions'])
        - `region_name` allows selecting a specific one

        ---- VARIABLE METADATA (from plot_var_dict) ----
        These are automatically pulled unless overridden:
        - cmap: colormap name
        - series: [min, max] range for colorbar
        - cmap_reverse: whether to reverse colormap
        - units: units string for colorbar
        - cbar_str: label for colorbar

        ---- PyGMT DISPLAY SETTINGS (from pygmt_dict) ----
        These can be overridden via kwargs, otherwise they fall back to config["pygmt_dict"]:
        - cbar_pos, sq_size_var, ice_shelves, water_color, land_color, fill_color, shoreline_str
        - text_font, text_justify, text_fill, text_pen

        ---- OUTPUT CONTROL ----
        - save_fig: if True, saves the PNG
        - ow_fig: if True, overwrites existing file
        - show_fig: if True, displays figure using `fig.show()`
        - P_png: manually specify output path
        - png_name_extras: list of strings appended to PNG name
        """
        # Set primary variable
        var_name = var_name if var_name is not None else self.var_name
        self.var_name = var_name
        if self.var_name not in self.plot_var_dict:
            raise ValueError(f"❌ Unknown variable '{self.var_name}'. Must be defined in `plot_var_dict`:\n{self.plot_var_dict.keys()}")
        var_dict   = self.plot_var_dict[var_name]
        pygmt_dict = self.pygmt_dict
        self.region_name   = region_name
        self.single_figure = single_figure
        self.plain_frame   = plain_frame
        self.frame         = frame
        self.P_png         = P_png
        self.plot_date_str = None
        self.text_str      = text_str
        # Non-kwarg fallback settings
        self.time_index = time_index if time_index is not None else 0
        self.dt0_str    = dt0_str    if dt0_str    is not None else getattr(self, "dt0_str")
        self.dtN_str    = dtN_str    if dtN_str    is not None else getattr(self, "dtN_str")
        self.ow_fig     = ow_fig     if ow_fig     is not None else getattr(self, "ow_fig")
        self.save_fig   = save_fig   if save_fig   is not None else getattr(self, "save_fig")
        self.show_fig   = show_fig   if show_fig   is not None else getattr(self, "show_fig")
        # General kwargs
        self.png_name_extras = kwargs.get("png_name_extras", None)
        # Class-level fallback values
        self.ice_type          = kwargs.get("ice_type"           , getattr(self, "ice_type"))
        self.plot_type         = kwargs.get("plot_type"          , getattr(self, "plot_type"))
        self.sim_name          = kwargs.get("sim_name"           , getattr(self, "sim_name"))
        self.D_graph           = kwargs.get("graphical_directory", getattr(self, "D_graph" ))
        self.hemisphere        = kwargs.get("hemisphere"         , getattr(self, "hemisphere"))
        self.region_dictionary = kwargs.get("region_dictionary"  , getattr(self, "reg_dict"))
        self.lon_coord_name    = kwargs.get("lon_coord_name"     , self.CICE_dict["proc_lon_coord"])
        self.lat_coord_name    = kwargs.get("lat_coord_name"     , self.CICE_dict["proc_lat_coord"])
        self.time_coord_name   = kwargs.get("time_coord_name"    , self.CICE_dict["proc_time_dim"])
        # Variable-dependent (from plot_var_dict)
        self.cmap         = kwargs.get("cmap"        , var_dict.get("cmap", "viridis"))
        self.series       = kwargs.get("series"      , var_dict.get("series", [0, 1]))
        self.cmap_reverse = kwargs.get("cmap_reverse", var_dict.get("reverse", False))
        self.units        = kwargs.get("units"       , var_dict.get("units", ""))
        self.cbar_str     = kwargs.get("cbar_str"    , var_dict.get("name", var_name))
        # PyGMT-style settings (from pygmt_dict)
        self.cbar_pos         = kwargs.get("cbar_pos"        , pygmt_dict.get("cbar_pos", "JBC+w10c/0.5c+mc+h"))
        self.sq_size_var      = kwargs.get("sq_size_var"     , pygmt_dict.get("sq_size_var", 0.15))
        self.sq_size_GI       = kwargs.get("sq_size_GI"      , pygmt_dict.get("sq_size_GI", 0.2))
        self.plot_ice_shelves = kwargs.get("plot_ice_shelves", pygmt_dict.get("ice_shelves", True))
        self.water_color      = kwargs.get("water_color"     , pygmt_dict.get("water_color", "white"))
        self.land_color       = kwargs.get("land_color"      , pygmt_dict.get("land_color", "seashell"))
        self.fill_color       = kwargs.get("fill_color"      , pygmt_dict.get("fill_color", "red"))
        self.shoreline_str    = kwargs.get("shoreline_str"   , pygmt_dict.get("shoreline_str", ".2p,white"))
        self.text_font        = kwargs.get("text_font"       , pygmt_dict.get("text_font", "12p,Courier-Bold,black"))
        self.text_justify     = kwargs.get("text_justify"    , pygmt_dict.get("text_justify", "LM"))
        self.text_fill        = kwargs.get("text_fill"       , pygmt_dict.get("text_fill", "white"))
        self.text_pen         = kwargs.get("text_pen"        , pygmt_dict.get("text_pen", "0.5p,black"))
        # 
        dt0_obj = pd.to_datetime(dt0_str)
        dtN_obj = pd.to_datetime(dtN_str)
        # OPTION 1:
        if ds is not None:
            print("using plot_map with option 1: dataset provided")
            if self.var_name not in ds:
                raise KeyError(f"Variable '{var_name}' not found in dataset\n{ds}")
            if self.var_name in self.spatial_distribution_vars:
                self.aggergate_2D_spatial_distribution(ds=ds, P_zarr_list=P_zarr_list)
                self.regional_or_hemisphere()
                return
            if self.time_index is not None:
                print(f"time index provided as index {self.time_index}, plot_map will use this on time dimension '{self.time_coord_name}' to create a single figure")
                da  = ds[self.var_name].isel({self.time_coord_name : self.time_index}).squeeze().values
                lat = ds[self.lat_coord_name].values
                lon = ds[self.lon_coord_name].values
                self.dataarray_2D_to_dataframe( da=da, lat=lat, lon=lon)
                self.plot_date_str = self.plot_date_str if self.plot_date_str is not None else pd.to_datetime(ds[self.time_coord_name].values[0]).strftime("%Y-%m-%d")
                self.regional_or_hemisphere()
                return
            elif self.single_figure: # find the nearest date ... single figure
                print(f"single figure switched on, plot_map will use dt0_str '{self.dt0_str}' to select from the time coordinate '{self.time_coord_name}' to create a single figure")
                da  = ds[self.var_name].sel({self.time_coord_name : dt0_str}, method='nearest').squeeze().values
                lat = ds[self.lat_coord_name].values
                lon = ds[self.lon_coord_name].values
                self.dataarray_2D_to_dataframe( da=da, lat=lat, lon=lon)
                self.plot_date_str = self.plot_date_str if self.plot_date_str is not None else pd.to_datetime(ds[self.time_coord_name].values[0]).strftime("%Y-%m-%d")
                self.regional_or_hemisphere()
                return
            else: # loop over dataaset time coordinate
                print(f"looping over time coordinate '{self.time_coord_name}' and using dt0_str '{self.dt0_str}' and dtN_str '{selt.dtN_str}' to determine which figures to create")
                for i in range(len(ds[self.time_coord_name].values)):
                    lat = ds[self.lat_coord_name].values
                    lon = ds[self.lon_coord_name].values
                    t   = ds[self.time_coord_name].isel({self.time_coord_name: i}).values
                    da  = ds[self.var_name].isel({self.time_coord_name: i}).squeeze().values
                    if dt0_obj <= pd.to_datetime(t) <= dtN_obj:
                        self.dataarray_2D_to_dataframe( da=da, lat=lat, lon=lon)
                        self.plot_date_str = pd.to_datetime(t).strftime("%Y-%m-%d")
                        print(f"plotting for date {plot_date_str}")
                        self.regional_or_hemisphere()
                    else:
                        print(f"{t} not in period")
                return
        # OPTION 2:
        elif P_zarr is not None:
            print("using plot_map with option 2: single zarr file path given")
            print(f"loading zarr file: {P_zarr}")
            ds = xr.open_dataset(P_zarr, engine='zarr', consolidated=True)
            if self.var_name not in ds:
                raise KeyError(f"Variable '{self.var_name}' not found in {ds}")
            da  = ds[self.var_name].squeeze().values
            lat = ds[self.lat_coord_name].values
            lon = ds[self.lon_coord_name].values
            self.dataarray_2D_to_dataframe( da=da, lat=lat, lon=lon)
            self.plot_date_str = self.plot_date_str if self.plot_date_str is not None else pd.to_datetime(ds[self.time_coord_name].values[0]).strftime("%Y-%m-%d")
            self.regional_or_hemisphere()
            return
        # OPTION 3:
        elif P_zarr_list is not None:
            print("using plot_map with option 3: list of zarr file paths given")
            for P_zarr in P_zarr_list:
                print(f"loading zarr file: {P_zarr}")
                ds = xr.open_dataset(P_zarr, engine='zarr', consolidated=True)
                if self.var_name not in ds:
                    raise KeyError(f"Variable '{self.var_name}' not found in {ds}")
                da  = ds[self.var_name].squeeze().values
                lat = ds[self.lat_coord_name].values
                lon = ds[self.lon_coord_name].values
                self.dataarray_2D_to_dataframe( da=da, lat=lat, lon=lon)
                self.plot_date_str = self.plot_date_str if self.plot_date_str is not None else pd.to_datetime(ds[self.time_coord_name].values[0]).strftime("%Y-%m-%d")
                self.regional_or_hemisphere()
            return
        # OPTION 4:
        else:
            print("using plot_map with option 4: dt0_str and dtN_str to load list of zarr files")
            zarr_pairs  = self._find_nearest_zarr_file(start_date_str=dt0_str, end_date_str=dtN_str)
            P_zarr_list = [p for p, _ in zarr_pairs]
            if self.var_name in self.spatial_distribution_vars:
                self.aggergate_2D_spatial_distribution(ds=ds, P_zarr_list=P_zarr_list)
                self.regional_or_hemisphere()
            else:
                for P_zarr in P_zarr_list:
                    print(f"loading zarr file: {P_zarr}")
                    ds = xr.open_dataset(P_zarr, engine='zarr', consolidated=True)
                    if self.var_name not in ds:
                        raise KeyError(f"Variable '{self.var_name}' not found in {ds}")
                    da  = ds[self.var_name].squeeze().values
                    lat = ds[self.lat_coord_name].values
                    lon = ds[self.lon_coord_name].values
                    self.dataarray_2D_to_dataframe( da=da, lat=lat, lon=lon)
                    self.plot_date_str = pd.to_datetime(ds[self.time_coord_name].values[0]).strftime("%Y-%m-%d")
                    self.regional_or_hemisphere()
            return

    def plot_grounded_icebergs(self):
        """
        Option 1: plot 'raw' grounded iceberg locations per regional dictionary; no simulation name required -> saves to grounded_iceberg director
        Option 2: plot grounded iceberg grid cell counts per regional dictionary; simulation name required -> saves to simulation directory
        Option 3: plot grounded iceberg grid cells; simulation name required -> saves to simulation directory
        """
        if P_KMT1 is None:
            GI_proc = GroundedIcebergProcessor()
            GI_proc.load_grid_and_landmask()
            kmt1 = GI_proc.KMT_org
        else: # assume KMT file is CICE original format (may not be but that's OK)
            GI_proc = GrounededIcebergProcessor()
            GI_proc.load_grid_and_landmask()
            kmt1 = xr.DataArray(data   = xr.open_dataset(P_KMT1).kmt.values,
                                dims   = self.CICE_dict.spatial_dims,
                                coords = {'lat': (self.CICE_dict.spatial_dims, GI_proc.G_t['lat'].values),
                                          'lon': (self.CICE_dict.spatial_dims, GI_proc.G_t['lon'].values)},
                                name   = 'kmt' )
        if P_KMT2 is None:
            GI_proc = GroundedIcebergProcessor(sim_name=self.sim_name)
            GI_proc.load_AFIM_GI()
            kmt2 = GI_proc.KMT_mod
        else: # assume KMT file is CICE original format (may not be but that's OK)
            GI_proc = GrounedIcebergProcessor()
            GI_proc.load_grid_and_landmask()
            kmt1 = xr.DataArray(data   = xr.open_dataset(P_KMT2).kmt.values,
                                dims   = self.CICE_dict.spatial_dims,
                                coords = {'lat': (self.CICE_dict.spatial_dims, GI_proc.G_t['lat'].values),
                                          'lon': (self.CICE_dict.spatial_dims, GI_proc.G_t['lon'].values)},
                                name   = 'kmt' )
        kmt_diffs = kmt1 - kmt2
        mask      = ((kmt_diffs == -1) | (kmt_diffs == 1))
        lon       = xr.DataArray(kmt_diffs['lon'].values, dims=kmt_diffs.dims)
        lat       = xr.DataArray(kmt_diffs['lat'].values, dims=kmt_diffs.dims)
        plot_lon  = lon.where(mask).values.flatten()
        plot_lat  = lat.where(mask).values.flatten()
        plot_dat  = kmt_diffs.where(mask).values.flatten()
        self.plot_df = pd.DataFrame({"lat": plot_lat,
                                     "lon": plot_lon,
                                     "dat": plot_dat})
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
            fig.plot(data=self.ice_shelves, pen="0.2p,gray", fill="gainsboro")

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
                fig.plot(data=self.ice_shelves, pen="0.2p,gray", fill="gainsboro")

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

    def plot_timeseries_groupby_month(self, ds=None, sim_name=None, var_name=None, dt0_str=None, dtN_str=None):
        sim_name = 'baseline'
        ds = xr.open_mfdataset(f"/g/data/gv90/da1339/afim_output/{sim_name}/FI/fast_ice_199*.zarr", engine='zarr')
        time = pd.to_datetime(ds['t_dim'].values)
        fia = ds['FIA']
        fia_obs = ds['FIA_OBS'].sel(sector='circumpolar')
        # --- Model: Monthly cycle by year ---
        fia_df = pd.DataFrame({'FIA': fia.values}, index=time)
        fia_df['Year'] = fia_df.index.year
        fia_df['Month'] = fia_df.index.month
        monthly_fia = fia_df.groupby(['Year', 'Month']).mean().reset_index()
        monthly_cycle = monthly_fia.pivot(index='Month', columns='Year', values='FIA')
        # --- Obs: Monthly climatology (mean over available time) ---
        # If FIA_OBS is already time-averaged, we assume it has a 12-month length
        obs_df = fia_obs.to_dataframe().dropna().reset_index()
        # Try parsing month info from the time (or fallback if no datetime index)
        if 't_dim' in obs_df.columns and np.issubdtype(obs_df['t_dim'].dtype, np.datetime64):
            obs_df['Month'] = obs_df['t_dim'].dt.month
        elif 'Month' not in obs_df.columns:
            obs_df['Month'] = range(1, len(obs_df)+1)
            monthly_obs = obs_df.groupby('Month')['FIA_OBS'].mean()
        # --- Plot ---
        plt.figure(figsize=(12, 6))
        plt.style.use('ggplot')
        # Plot each year from model
        for year in monthly_cycle.columns:
            plt.plot(monthly_cycle.index, monthly_cycle[year], label=f"{year}", linewidth=2)
        # Plot observed climatology
        plt.plot(monthly_obs.index, monthly_obs.values, label='Observed FIA (AF2020db)', 
                 color='black', linestyle='--', linewidth=3)
        plt.title("Monthly Fast Ice Area Cycle (1993–1999)", fontsize=16)
        plt.xlabel("Month", fontsize=14)
        plt.ylabel(f"Fast Ice Area ({fia.attrs.get('units', 'unknown')})", fontsize=14)
        plt.xticks(ticks=range(1,13), labels=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'])
        plt.grid(True, linestyle='--', alpha=0.5)
        plt.legend(loc='upper right', fontsize=10, title='Year')
        plt.tight_layout()
        plt.show()
        plt.savefig(f"/g/data/gv90/da1339/GRAPHICAL/timeseries/{sim_name}_FIA_grouped_by_month.png")

    def plot_timeseries(self, var_names=None, label=None, save_png=True):
        """
        Plot sea ice time series for a single simulation.

        For pack ice (PI), it plots two subplots:
        - Top: PIA and SIA
        - Bottom: PIE and SIE

        If Zarr files are missing, will auto-generate using PackIceProcessor.
        """
        if self.dtN_str is None:
            raise ValueError("Time series plotting requires a defined dtN_str.")
        if var_names is None:
            if self.ice_type == "PI":
                var_names = ["PIA", "SIA", "PIE", "SIE"]
            else:
                var_names = [f"{self.ice_type}A"]
        all_time = []
        data_dict = {var: [] for var in var_names}
        if not getattr(self, 'P_zarr_list', []):
            print(f"⏳ No Zarr time series found, generating for vars: {var_names}")
            ds = self._process_pack_ice(self.sim_name, var_names)
            if "time" not in ds:
                raise KeyError("Processed dataset missing 'time' coordinate")
            times = pd.to_datetime(ds["time"].values)
            all_time.extend(times)
            for var in var_names:
                if var not in ds:
                    print(f"⚠️ Variable '{var}' missing from processed dataset.")
                    data_dict[var] = [np.nan] * len(times)
                else:
                    data_dict[var] = ds[var].values.tolist()
        else:
            for P_zarr in self.P_zarr_list:
                ds = xr.open_zarr(P_zarr)
                if "time" not in ds:
                    raise KeyError(f"'time' coordinate not found in {P_zarr}")
                times = pd.to_datetime(ds["time"].values)
                all_time.extend(times)
                for var in var_names:
                    values = ds.get(var, np.full(len(times), np.nan))
                    data_dict[var].extend(values)

        df = pd.DataFrame(data=data_dict, index=pd.to_datetime(all_time)).sort_index()

        if self.ice_type == "PI":
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
            fig.subplots_adjust(hspace=0.3)

            # Top: PIA and SIA
            for var in ["PIA", "SIA"]:
                if var in df:
                    plot_label = f"{label} ({var})" if label else var
                    ax1.plot(df.index, df[var], label=plot_label)
            ax1.set_ylabel("Area (10³ km²)")
            ax1.set_title(f"{self.sim_name} Pack Ice Area (PIA/SIA)")
            ax1.grid(True, linestyle="--", alpha=0.5)
            ax1.legend()

            # Bottom: PIE and SIE
            for var in ["PIE", "SIE"]:
                if var in df:
                    plot_label = f"{label} ({var})" if label else var
                    ax2.plot(df.index, df[var], label=plot_label)
            ax2.set_xlabel("Date")
            ax2.set_ylabel("Extent (10³ km²)")
            ax2.set_title(f"{self.sim_name} Pack Ice Extent (PIE/SIE)")
            ax2.grid(True, linestyle="--", alpha=0.5)
            ax2.legend()

            if self.show_fig and not save_png:
                plt.show()
            else:
                F_save = f"{self.sim_name}_PIA_PIE_{self.dt0_str}_to_{self.dtN_str}.png"
                P_save = Path(self.D_ts, F_save)
                plt.savefig(P_save)
                print(f"✅ Saved dual subplot to: {P_save}")

        else:
            # Fallback: single timeseries plot for non-PI types
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

            if self.show_fig and not save_png:
                plt.show()
            else:
                var_str = "_".join(var_names)
                F_save = f"{self.sim_name}_{var_str}_{self.dt0_str}_to_{self.dtN_str}.png"
                P_save = Path(self.D_ts, F_save)
                plt.savefig(P_save)
                print(f"✅ Saved plot to: {P_save}")

    def plot_timeseries_compare(self, sim_names, var_name=None, label_dict=None,
                                save_png=True, comparison_name="comparison"):
        """
        Compare a variable across multiple simulations. If var_name == 'PIA',
        will also include NSIDC SIA by default if configured.
        """
        if self.dtN_str is None:
            raise ValueError("Comparison requires a defined dtN_str at SeaIcePlotter init.")

        var_name = var_name or f"{self.ice_type}A"
        plt.figure(figsize=(10, 6))

        for sim in sim_names:
            t0 = time.time()
            print(f"🔄 Processing simulation {sim} for variable {var_name}")
            sim_plotter = SeaIcePlotter(sim, self.dt0_str, self.dtN_str, ice_type=self.ice_type)
            if not sim_plotter.check_zarr_exists(var_name):
                print(f"ℹ️ Zarr for '{var_name}' not found, generating with PackIceProcessor...")
                from AFIM.src.pack_ice_processor import PackIceProcessor
                proc = PackIceProcessor(sim_name=sim, dt0_str=self.dt0_str, dtN_str=self.dtN_str)
                proc.process_window(var_list=[var_name], save_zarr=True)

            time_list, value_list = [], []
            for P_zarr in sim_plotter.P_zarr_list:
                ds = xr.open_zarr(P_zarr)
                if "time" in ds and var_name in ds:
                    time_list.extend(pd.to_datetime(ds["time"].values))
                    value_list.extend(ds[var_name].values)
            if not time_list:
                print(f"⚠️ No data loaded for simulation: {sim}")
                continue
            df = pd.DataFrame({var_name: value_list}, index=pd.to_datetime(time_list)).sort_index()
            label = label_dict.get(sim, sim) if label_dict else sim
            plt.plot(df.index, df[var_name], label=label)
            print(f"✅ Added {sim}.{var_name} to plot in {time.time()-t0:.2f} seconds")

        # NSIDC overlay
        if var_name == "PIA" and self.config["NSIDC_dict"].get("plot_SIA", True):
            self._plot_nsidc_overlay()

        plt.xlabel("Date")
        plt.ylabel(var_name)
        plt.title(f"{self.ice_type} Time Series Comparison: {var_name}")
        plt.grid(True, linestyle="--", alpha=0.5)
        plt.legend()
        if self.show_fig and not save_png:
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

