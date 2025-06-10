import os, sys, time, json, imageio, shutil, pygmt, imageio, shutil
import xarray as xr
import pandas as pd
import numpy as np
import geopandas as gpd
#import geopandas as gpd
from pathlib import Path
from datetime import datetime
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
sys.path.insert(0, '/home/581/da1339/AFIM/src/AFIM/src')
from sea_ice_science import SeaIceScience
from PIL import Image

class SeaIcePlotter(SeaIceScience):
    def __init__(self, sim_name, **kwargs):
        super().__init__(sim_name, **kwargs)
    # """
    # A plotting class for visualizing sea ice model output (e.g., from CICE or AFIM) in various formats
    # including regional maps, hemisphere-wide views, and time series plots.

    # The class supports output from different ice types (fast ice, pack ice, sea ice, grounded icebergs)
    # and integrates configuration from a JSON settings file to control plotting behavior, file paths,
    # and variable metadata.

    # Key features:
    # - Supports plotting from in-memory xarray datasets or from Zarr files
    # - Integrates with grounded iceberg data if available
    # - Regional and hemisphere plotting via PyGMT
    # - Loads plotting styles from `pygmt_dict` and variable metadata from `plot_var_dict`
    # - Fully configurable through JSON to adapt to different simulations and data products
    # """
    # def __init__(self,
    #              P_json     = None,
    #              sim_name   = None,
    #              dt0_str    = None,
    #              dtN_str    = None,
    #              overwrite  = False,
    #              save_fig   = False,
    #              show_fig   = False,
    #              hemisphere = None):
    #     """
    #     Initialize a SeaIcePlotter instance using configuration provided in a JSON file.

    #     Parameters
    #     ----------
    #     P_JSON : str or Path, optional
    #         Path to the configuration JSON file. If None, defaults to a pre-defined internal path.
    #     sim_name : str, optional
    #         Name of the simulation. Controls file paths and Zarr subdirectories.
    #     ice_type : str, optional
    #         Type of ice data ('FI', 'PI', 'SI', 'GI'). Controls dictionary and data sub-paths.
    #     plot_type : str, optional
    #         Type of plot to generate: 'regional', 'hemisphere', or 'timeseries'.
    #     var_name : str, optional
    #         Name of the variable to plot. Must exist in `plot_var_dict` in the config JSON.
    #     dt0_str : str, optional
    #         Start date (inclusive) for plots, in 'YYYY-MM-DD' format. Defaults to value in config.
    #     dtN_str : str, optional
    #         End date (inclusive) for plots, in 'YYYY-MM-DD' format. Defaults to value in config.
    #     overwrite : bool, optional
    #         If True, overwrite existing figures. Defaults to False.
    #     save_fig : bool, optional
    #         If True, save the generated figures to PNG. Defaults to False.
    #     show_fig : bool, optional
    #         If True, display figures after plotting. Defaults to True.
    #     hemisphere : str, optional
    #         Geographic hemisphere ('north' or 'south'). Determines projection and map extent.

    #     Notes
    #     -----
    #     - Loads multiple dictionaries from JSON, including:
    #         * D_dict       → file paths
    #         * CICE_dict    → coordinate names
    #         * pygmt_dict   → default PyGMT styling
    #         * plot_var_dict→ plotting metadata per variable
    #         * AF_regions   → predefined regional extents
    #     - Initializes a GroundedIcebergProcessor to overlay grounded iceberg locations if enabled.
    #     """
    #     if P_json is None:
    #         P_json = "/home/581/da1339/AFIM/src/AFIM/src/JSONs/afim_cice_analysis.json"
    #     with open(P_json, 'r') as f:
    #         self.config = json.load(f)
    #     self.sim_name      = sim_name
    #     self.dt0_str       = dt0_str          if dt0_str    is not None else self.config.get('dt0_str', '1993-01-01')
    #     self.dtN_str       = dtN_str          if dtN_str    is not None else self.config.get('dtN_str', '1999-12-31')
    #     self.ow_fig        = overwrite        if overwrite  is not None else False
    #     self.save_fig      = save_fig         if save_fig   is not None else False
    #     self.show_fig      = show_fig         if show_fig   is not None else True
    #     self.hemisphere    = hemisphere       if hemisphere is not None else self.config.get('hemisphere', 'south')
    #     self.CICE_dict     = self.config.get("CICE_dict", {})
    #     self.pygmt_dict    = self.config.get("pygmt_dict", {})
    #     self.plot_var_dict = self.config.get("plot_var_dict", {})
    #     self.reg_dict      = self.config.get('AF_regions', {})
    #     self.sim_dict      = self.config.get("sim_dict", {})
    #     self.GI_dict       = self.config.get("GI_dict", {})
    #     self.SIC_scale     = self.config.get("SIC_scale", 1e12)
    #     self.D_graph       = Path(self.config['D_dict']['graph'], 'AFIM')
    #     self.define_hemisphere(self.hemisphere)
    #     self.P_KMT_org           = Path(self.GI_dict["D_GI_thin"],self.GI_dict['KMT_org_fmt'])
    #     self.GI_version          = self.sim_config.get('GI_version')
    #     self.GI_version_str      = f"{self.GI_version:.2f}".replace('.', 'p')
    #     if self.GI_thin_fact>0:
    #         self.use_gi    = True
    #         self.P_KMT_mod = os.path.join(self.GI_dict['D_GI_thin'],
    #                                       self.GI_dict['KMT_mod_fmt'].format(GI_thin_fact = self.GI_thin_str,
    #                                                                          version      = self.GI_version_str))

    # def define_hemisphere(self, hemisphere):
    #     """
    #     THIS FUNCTION SHOULD BE *EXACTLY* THE SAME AS THE ONE IN SEAICEPROCESSOR CLASS
    #     """
    #     if hemisphere.lower() in ['north', 'northern', 'nh', 'n', 'no']:
    #         self.hemisphere_geographic_extent = [0, 360, 0, 90]
    #         self.hemisphere_map_extent        = [-180,180,55,90]
    #         self.hemisphere_projection        = "S0.0/90.0/50/15C"
    #         self.hemisphere_map_text_location = [-120,56]
    #         self.hemisphere_abbreviation      = 'NH'
    #         self.hemisphere_nj_slice          = slice(540,1080)
    #         self.hemisphere                   = 'north'
    #     elif hemisphere.lower() in ['south', 'southern', 'sh', 's', 'so']:
    #         self.hemisphere_geographic_extent = [0, 360, -90, 0]
    #         self.hemisphere_map_extent        = [-180,180,-90,-55]
    #         self.hemisphere_projection        = "S0.0/-90.0/50/15C"
    #         self.hemisphere_map_text_location = [0,-90]
    #         self.hemisphere_abbreviation      = 'SH'
    #         self.hemisphere_nj_slice          = slice(0,540)
    #         self.hemisphere                   = 'south'
    #     else:
    #         raise ValueError(f"Invalid hemisphere '{hemisphere}'. Valid options are: "
    #                          "['north', 'south', 'northern', 'southern', 'sh', 'nh', 'SH', 'NH']")

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
            return [f"a{ann_str}f{tick_str}+l{label}",f"y+l {units}"]
        else:
            return f"a{ann_str}f{tick_str}+l{label}"

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

    def create_figure(self):
        if self.P_png is None:
            self.define_png_plot_path(qualifiers = self.png_name_extras)
        if os.path.exists(self.P_png) and not self.ow_fig:
            print(f"figure already exists and not overwriting ... skipping {self.P_png}")
            self.P_png = None
            return
        if self.units is not None and self.cmap!="categorical":
            cbar_frame_fore = self.create_cbar_frame(self.series, self.cbar_str, units=self.units)
        if self.units_back is not None:
            cbar_frame_back = self.create_cbar_frame(self.series_back, self.cbar_str_back, units=self.units_back)
        if self.text_str is not None:
            print(f"\ttext location: {self.text_loc}")
        if self.plain_frame:
            frame = ["af"]
        elif self.frame is not None:
            frame = self.frame
        else:
            if self.title_extras is not None:
                frame = ["af", f"+t{self.var_name} {self.plot_date_str}, {self.title_extras}"]
            else:
                frame = ["af", f"+t{self.var_name} {self.plot_date_str}"]
        t0  = time.time()
        fig = pygmt.Figure()
        pygmt.config(FONT_TITLE="16p,Courier-Bold", FONT_ANNOT_PRIMARY="14p,Helvetica")
        #pygmt.config(NAN_COLOR="none")
        print(f"plotting figure with projection {self.projection} and frame {frame}")
        fig.basemap(projection=self.projection, region=self.region_extent, frame=frame)
        fig.coast(land=self.land_color, water=self.water_color)
        if self.plot_df_background is not None:
            print("plotting background graphic")
            pygmt.makecpt(cmap=self.cmap_back, reverse=self.cmap_reverse_back, series=self.series_back)
            fig.plot(x=self.plot_df_background.lon, y=self.plot_df_background.lat, cmap=True, fill=self.plot_df_background.dat, style=f"s{self.sq_size_var_back}c")
            fig.colorbar(position=self.cbar_pos_back, frame=cbar_frame_back)
        if not self.plot_df.empty:
            if self.var_name == "ispd":
                df      = self.plot_df.dropna(subset=["dat"])
                df_slow = df[(df['dat'] > 0) & (df['dat'] <= self.ispd_thresh)]
                df_fast = df[df['dat'] > self.ispd_thresh]
                if not df_slow.empty:
                    slow_min = max(df_slow['dat'].min(), 1e-10)
                    pygmt.makecpt(cmap="gray", series=[slow_min, self.ispd_thresh])
                    fig.plot(x=df_slow.lon, y=df_slow.lat, cmap=True, fill=df_slow.dat, style=f"s{self.sq_size_var}c")
                    cbar_frame_vert = [f"a{self.ispd_thresh/5}f{self.ispd_thresh/10}+lslow ice speed",f"y+l m/s"]
                    fig.colorbar(position="JMR+o1.2c/0c+w8c/0.5c+mc", frame=cbar_frame_vert)#f"af+l'Slow ice speed (≤ {self.ispd_thresh:.0e} m/s)'")
                if not df_fast.empty:
                    pygmt.makecpt(cmap=self.cmap, reverse=self.cmap_reverse, series=[self.ispd_thresh,self.series[1]])
                    fig.plot(x=df_fast.lon, y=df_fast.lat, cmap=True, fill=df_fast.dat, style=f"s{self.sq_size_var}c")
                    fig.colorbar(position=self.cbar_pos, frame=cbar_frame_fore)
            else:
                pygmt.makecpt(cmap=self.cmap, reverse=self.cmap_reverse, series=self.series)
                fig.plot(x=self.plot_df.lon, y=self.plot_df.lat, cmap=True, fill=self.plot_df.dat, style=f"s{self.sq_size_var}c")
                fig.colorbar(position=self.cbar_pos, frame=cbar_frame_fore)
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
            fig.plot(x=self.GI_proc.GI_lon_cells, y=self.GI_proc.GI_lat_cells, fill=self.fill_color, style=f"c{self.sq_size_GI}c")
        fig.coast(shorelines=self.shoreline_str)
        if self.show_fig:
            fig.show()
        if self.save_fig:
            fig.savefig(self.P_png)
            print(f"figure saved: {self.P_png}")
            self.P_png = None
        print(f"time taken {time.time()-t0:0.2f} seconds\n")
        pygmt.clib.Session.__exit__

    def plot_cice_field(self, x1, y1, z1,
                        projection   = "S75/-90/30c",
                        region       = None,
                        title        = None,
                        cmap         = "cmocean/haline",
                        cmap_series  = [0.9, 1],
                        reverse_cmap = False,
                        GI_coords    = None,
                        GI_color     = "red",
                        GI_style     = "c0.3c",
                        var_style    = "s0.5c",
                        cbar_label   = "sea ice concentration",
                        units        = "1/100",
                        plot_cbar    = True,
                        **kwargs):
        pygmt_dict       = self.pygmt_dict
        cbar_pos         = kwargs.get("cbar_pos"        , pygmt_dict.get("cbar_pos", "JBC+w10c/0.5c+mc+h"))
        plot_ice_shelves = kwargs.get("plot_ice_shelves", pygmt_dict.get("ice_shelves", True))
        water_color      = kwargs.get("water_color"     , pygmt_dict.get("water_color", "white"))
        land_color       = kwargs.get("land_color"      , pygmt_dict.get("land_color", "seashell"))
        shoreline_str    = kwargs.get("shoreline_str"   , pygmt_dict.get("shoreline_str", ".2p,white"))
        if title is not None:
            frame = ["af", f"+t{title}"]
        else:
            frame = ["af"]
        fig = pygmt.Figure()
        pygmt.config(FORMAT_GEO_MAP="ddd.x", FONT_TITLE="16p,Courier-Bold", FONT_ANNOT_PRIMARY="14p,Helvetica")
        pygmt.makecpt(cmap=cmap, reverse=reverse_cmap, series=cmap_series)
        fig.basemap(projection=projection, region=region, frame=frame)
        fig.coast(land=land_color, water=water_color)
        fig.plot(x=x1, y=y1, fill=z1, cmap=True, style=var_style)
        if GI_coords:
            fig.plot(x=GI_coords[0], y=GI_coords[1], fill=GI_color, style=GI_style)
        if plot_ice_shelves:
            if not hasattr(self, 'ice_shelves'):
                self.load_ice_shelves()
            fig.plot(data=self.ice_shelves, pen="0.2p,gray", fill="lightgray")
        fig.coast(shorelines=shoreline_str)
        if plot_cbar:
            fig.colorbar(position=cbar_pos, frame=[f"x+l{cbar_label}", f"y+l{units}"])
        return fig

    def plot_persistence_map(self, DA,
                             lon_coord_name = "TLON",
                             lat_coord_name = "TLAT",
                             tit_str        = None,
                             ispd_str       = None,
                             ice_type       = None,
                             sim_name       = None,
                             regional       = False,
                             plot_GI        = False,
                             plot_cbar      = True,
                             dt_range_str   = None,
                             overwrite_png  = False,
                             show_fig       = False):
        if plot_GI:
            from grounded_iceberg_processor import GroundedIcebergProcessor
            GI_proc   = GroundedIcebergProcessor(sim_name=self.sim_name)
            G_t       = GI_proc.align_modified_landmask()
            GI_coords = (G_t['GI_lon'], G_t['GI_lat'])
        else:
            GI_coords = None
        cmap = self.pygmt_dict.get("FIP_CPT")
        cbar_lab = "fast ice persistence"
        cmap_ser = [0.01,1,.01]
        da   = DA.values.ravel()
        lon  = DA[lon_coord_name].values.ravel()
        lat  = DA[lat_coord_name].values.ravel()
        mask = ~np.isnan(da)
        lat_plt, lon_plt, da_plt = lat[mask], lon[mask], da[mask]
        if regional:
            for reg_name, reg_cfg in self.reg_dict.items():
                if ispd_str is not None and ice_type is not None:
                    D_plt = Path(self.D_graph, sim_name, reg_name, "FIP", f"ispd_thresh_{ispd_str}", ice_type)
                else:
                    D_plt = Path(self.D_graph, sim_name, reg_name, "FIP")
                D_plt.mkdir(parents=True, exist_ok=True)
                P_plt = Path(D_plt,f"FIP_{dt_range_str}.png")
                if P_plt.exists() and not overwrite_png:
                    print(f"figure {P_plt} exists and not overwriting")
                    continue
                else:
                    MC  = self.get_meridian_center_from_geographic_extent(reg_cfg['plt_ext'])
                    fig = self.plot_cice_field(lon_plt, lat_plt, da_plt, 
                                               projection  = f"S{MC}/-90/30c",
                                               region      = reg_cfg['plt_ext'],
                                               title       = tit_str,
                                               cmap        = cmap,
                                               cmap_series = cmap_ser,
                                               reverse_cmap = False,
                                               var_style   = "s0.25c",
                                               cbar_label  = cbar_lab)
                    print(f"📸 Saving figure: {P_plt}")
                    fig.savefig(P_plt)
                    if show_fig:
                        fig.show()
        else:
            print("*** PLOTTING HEMISPHERE ***")
            if ispd_str is not None and ice_type is not None:
                D_plt = Path(self.D_graph, sim_name, "SO", "FIP", f"ispd_thresh_{ispd_str}", ice_type)
            else:
                D_plt = Path(self.D_graph, sim_name, "SO", "FIP")
            D_plt.mkdir(parents=True, exist_ok=True)
            P_plt = Path(D_plt,f"FIP_{dt_range_str}.png")
            if P_plt.exists() and not overwrite_png:
                print(f"figure {P_plt} exists and not overwriting")
            else:
                fig = self.plot_cice_field(lon_plt, lat_plt, da_plt,
                                           projection  = self.hemisphere_projection,
                                           region      = self.hemisphere_map_extent,
                                           title       = tit_str,
                                           cmap        = cmap,
                                           cmap_series = cmap_ser,
                                           reverse_cmap = False,
                                           var_style   = "s0.1c",
                                           cbar_label = cbar_lab)
                print(f"📸 Saving figure: {P_plt}")
                fig.savefig(P_plt)
                if show_fig:
                    fig.show()

    def animate_over_time(self, ds, var_name, D_mp4=None, region_name=None, fps=2, time_coordinate_name=None, clobber_temporaries=True, **kwargs):
        """
        Generate an animation (GIF) over the 'time' dimension using plot_map.
        Parameters
        ----------
        """
        def round_up_to_multiple(x, multiple=16):
            return int(np.ceil(x / multiple) * multiple)
        P_mp4      = D_mp4 if D_mp4 is not None else Path(self.D_graph,"animations","AFIM",self.sim_name,f"{var_name}_{region_name}.mp4")
        time_coord = time_coordinate_name if time_coordinate_name is not None else self.CICE_dict["time_dim"]
        kwargs['time_coord_name'] = time_coord
        D_frames   = Path(self.config['D_dict']['tmp'],"animation_frames")
        D_frames.mkdir(parents=True, exist_ok=True)
        if clobber_temporaries:
            for p in D_frames.glob("*.png"):
                p.unlink()
        n_times  = ds[var_name].sizes[time_coord]
        P_frames = []
        for i in range(n_times):
            plot_date_str = pd.to_datetime(ds[time_coord].isel({time_coord: i}).values).strftime("%Y-%m-%d")
            P_png         = str(Path(D_frames,f"frame_{i:03d}_{plot_date_str}.png"))
            print(f"Rendering frame: {i+1}/{n_times}: {plot_date_str}")
            self.plot_map(ds=ds, time_index=i, var_name=var_name, region_name=region_name, plot_date_str=plot_date_str, P_png=P_png, **kwargs)
            P_frames.append(P_png)
        P_frames = [p for p in P_frames if p is not None and Path(p).exists()]
        images   = [imageio.v3.imread(png) for png in P_frames]
        P_mp4.parent.mkdir(parents=True, exist_ok=True)
        print(f"🌀 Encoding MP4 from {len(P_frames)} frames ...")
        # Load first image to get target size
        target_image  = imageio.v3.imread(P_frames[0])
        target_width  = round_up_to_multiple(target_image.shape[1])
        target_height = round_up_to_multiple(target_image.shape[0])
        target_size   = (target_width, target_height)
        #target_size  = (target_image.shape[1], target_image.shape[0])  # (width, height)
        images = []
        for png in P_frames:
            img = imageio.v3.imread(png)
            if img.shape[:2] != target_size[::-1]:
                print(f"⚠️ Resizing image {png} from {img.shape[:2]} to {target_size[::-1]}")
                img = np.array(Image.fromarray(img).resize(target_size, Image.BILINEAR))
            images.append(img)
        imageio.mimsave(P_mp4, images, fps=fps, codec='libx264', macro_block_size=None)
        #imageio.mimsave(P_mp4, images, fps=fps, codec='libx264', quality=8)
        print(f"Animation saved to: {P_mp4}")
        if clobber_temporaries:
            for p in D_frames.glob("*.png"):
                p.unlink()

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
        fig.plot(x     = lon[mask],
                 y     = lat[mask],
                 style = f"s{sq_size}c",
                 fill  = dat[mask],
                 cmap  = True,
                 pen   = "thin" )
        fig.colorbar(position="JBC+w10c/0.5c+h", frame=["x+lGrounded Iceberg Count"])
        # Add ice shelf overlay if loaded
        if self.ice_shelves_loaded:
            fig.plot(data=self.ice_shelves, pen="0.2p,gray", fill="gainsboro")
        # Regional stats
        count_vals = dat[mask]
        stats = [f"Total Cells : {count_vals.size}",
                 f"Mean Count  : {np.mean(count_vals):.2f}",
                 f"Max Count   : {np.max(count_vals)}" ]
        if use_thinned and hasattr(self, 'GI_thin_fact'):
            stats.append(f"Thinning: {self.GI_thin_fact:.2%}")
        fig.text(x       = region[0] + 1,
                 y       = region[2] + 1,
                 text    = "\n".join(stats),
                 font    = "10p,Helvetica,black",
                 justify = "LM",
                 fill    = "white",
                 pen     = "0.5p,black" )
        fig.savefig(save_path)

        print(f"Saved: {save_path}")

    def _plot_fip_east_west(self, fig, DA, sim_name, plot_GI, lon_coord_name="TLON", lat_coord_name="TLAT"):
        ANT_IS = gpd.read_file(self.pygmt_dict['P_coast_shape'])
        ANT_IS = ANT_IS[(ANT_IS['POLY_TYPE'] == 'S') & ANT_IS.geometry.notnull()].to_crs("EPSG:4326")
        ANT_IS.geometry = ANT_IS.geometry.buffer(0)

        DA_mask = (DA.data > 0) & np.isfinite(DA.data)
        dat2d = DA.data
        lon2d = DA[lon_coord_name].data
        lat2d = DA[lat_coord_name].data
        plt_dat = dat2d[DA_mask].ravel()
        plt_lon = lon2d[DA_mask].ravel()
        plt_lat = lat2d[DA_mask].ravel()

        if plot_GI:
            G_t = self.define_cice_grid(grid_type='t').isel(nj=slice(0,540))
            kmt_mod = xr.open_dataset(self.P_KMT_mod).isel(nj=slice(0,540)).kmt.data
            kmt_org = xr.open_dataset(self.P_KMT_org).isel(nj=slice(0,540)).kmt.data
            GI_mask = (kmt_org == 1) & (kmt_mod == 0)
            plt_GI_lon = lon2d[GI_mask].ravel()
            plt_GI_lat = lat2d[GI_mask].ravel()

        pygmt.makecpt(cmap=self.pygmt_dict.get("FIP_CPT"), series=[0.01, 1.0, 0.01], continuous=True)

        for region, proj in [([-20, 160, -73, -62], "M70/-70/30"),
                             ([160, 340, -79, -61], "M70/-70/30")]:
            fig.shift_origin(yshift="-10.5c")
            fig.basemap(region=region, projection=proj, frame=["af"])
            fig.plot(x=plt_lon, y=plt_lat, fill=plt_dat, style="s0.2c", cmap=True)
            if plot_GI:
                fig.plot(x=plt_GI_lon, y=plt_GI_lat, fill="red", style="c0.05c")
            fig.plot(data=ANT_IS, pen="0.2p,gray", fill="lightgray")

        fig.colorbar(position="JBC+w10c/0.5c+mc+h", frame=["x+lFast Ice Persistence", "y+l1/100"])

    def plot_FIA_FIP_faceted(self, FIA_dict, FIP_DA,
                              sim_name       = None,
                              dt_range_str   = None,
                              P_png          = None,
                              plot_GI        = False,
                              ylim           = (0, 1000),
                              roll_days      = 0,
                              lon_coord_name = None,
                              lat_coord_name = None,
                              overwrite_png  = False,
                              show_fig       = False):
        """
        Composite plot showing monthly climatology of FIA and east/west persistence map.
        """
        sim_name = sim_name or self.sim_name

        # ----- Create and save x-axis annotation file for months -----
        xticks = [1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 360]
        month_labels = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
                        "Sep", "Oct", "Nov", "Dec", "Dec"]
        xannot_path = Path("xannots.txt")
        xannot_lines = [f"{tick}\tig\t{label}" for tick, label in zip(xticks, month_labels)]
        xannot_path.write_text("\n".join(xannot_lines) + "\n")

        # ----- Plot FIA climatology -----
        fig = pygmt.Figure()
        pygmt.config(FORMAT_GEO_MAP="ddd.x", MAP_GRID_PEN="0p,white")
        fig.basemap(
            region=[1, 360, 100, ylim[1]],
            projection="X30c/15c",
            **{"frame": [
                f"sxc{xannot_path}",
                f"ya100f+lFast Ice Area (1000 km\u00b2)",
                "WSne"]}
        )

        plot_style = {"FI_BT":       {"label": "daily ice-speed mask", "color": "black"},
                      "FI_BT_bool":  {"label": "Binary-days (6 of 7) ice-speed mask", "color": "orange"},
                      "FI_BT_roll":  {"label": "15-day-avg ice-speed mask", "color": "darkgreen"},}

        for i, (name, da) in enumerate(FIA_dict.items()):
            if isinstance(da, xr.Dataset):
                da = da.to_array().squeeze()

            if name == "AF2020db_cli":
                doy = da["doy"].values
                values = da.values
                fig.plot(x=doy, y=values, pen="1.5p,blue,--", label="AF2020 Climatology (2000–2018)")
                continue

            time = pd.to_datetime(da["time"].values)
            area = da.rolling(time=roll_days, center=True, min_periods=1).mean().values if roll_days else da.values
            df = pd.DataFrame({"area": area}, index=time)
            df["doy"] = df.index.dayofyear
            df["year"] = df.index.year
            df = df[df["year"] > df["year"].min()]
            grouped = df.groupby("doy")["area"]
            area_min = grouped.min()
            area_max = grouped.max()
            area_mean = grouped.mean()

            style = plot_style.get(name, {"label": name, "color": f"tab{i}"})
            fig.plot(x=np.concatenate([area_min.index, area_max.index[::-1]]),
                     y=np.concatenate([area_min.values, area_max.values[::-1]]),
                     fill=f"{style['color']}@70", close=True, transparency=60)
            fig.plot(x=area_mean.index, y=area_mean.values,
                     pen=f"2p,{style['color']}", label=f"{style['label']} climatology")

        fig.legend(position="JTL+jTL+o0.2c", box="+gwhite+p.5p")

        # ----- Append east/west FIP maps -----
        fig.shift_origin(yshift="-6c")
        self._plot_fip_east_west(fig, FIP_DA, sim_name, plot_GI, lon_coord_name=lon_coord_name, lat_coord_name=lat_coord_name)

        if P_png:
            P_png.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(P_png)
            print(f"📏 Saved figure to {P_png}")
        if show_fig:
            fig.show()

    def plot_sectors_grounded_icebergs(self, 
                                       KMT1=None, 
                                       KMT2=None,
                                       n_sectors           = 9,
                                       sector_width        = 40,
                                       lat_limit           = -60,
                                       start_lon           = 0,
                                       marker_size         = "0.1c",
                                       save                = False,
                                       save_dir            = None,
                                       filename_extras     = None,
                                       overlay_ice_shelves = False):
        """
        """
        KMT1 = KMT1 if KMT1 is not None else xr.open_dataset(self.P_KMT_org)
        KMT2 = KMT2 if KMT2 is not None else xr.open_dataset(self.P_KMT_mod)
        da1  = KMT1['kmt'].values
        da2  = KMT2['kmt'].values
        lat1 = KMT1['lat'].values
        lon1 = KMT1['lon'].values
        lon1 = (lon1 + 360) % 360
        lat2 = KMT2['lat'].values
        lon2 = KMT2['lon'].values
        lon2 = (lon2 + 360) % 360
        diff_mask = (da1 == 1) & (da2 == 0)
        print("New grounded iceberg cells:", np.sum(diff_mask))
        nj, ni = np.where(diff_mask)
        lon1_pts = lon1[nj, ni]
        lat1_pts = lat1[nj, ni]
        for i in range(n_sectors):
            fig = pygmt.Figure()
            lon_min = (start_lon + i * sector_width) % 360
            lon_max = (lon_min + sector_width) % 360
            if lon_min < lon_max:
                meridian_center = (lon_min + lon_max) / 2
            else:
                meridian_center = ((lon_min + lon_max + 360) / 2) % 360
            if lon_min < lon_max:
                sector_mask = (lon1_pts >= lon_min) & (lon1_pts < lon_max) & (lat1_pts <= lat_limit)
            else:
                sector_mask = ((lon1_pts >= lon_min) | (lon1_pts < lon_max)) & (lat1_pts <= lat_limit)
            sector_lons = lon1_pts[sector_mask]
            sector_lats = lat1_pts[sector_mask]
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
            stat_text = [f"# GI cells: {len(sector_lons)}"]
            if len(sector_lats) > 0:
                stat_text.append(f"Lat range: {np.min(sector_lats):.1f}° to {np.max(sector_lats):.1f}°")
                stat_text.append(f"Lon range: {int(lon_min)}° to {int(lon_max)}°")
            fig.text(x       = region[0] + 1,
                     y       = region[2] + 1,
                     text    = "\n".join(stat_text),
                     font    = "10p,Helvetica,black",
                     justify = "LM",
                     fill    = "white",
                     pen     = "0.5p,black")
            if overlay_ice_shelves and self.ice_shelves_loaded:
                fig.plot(data=self.ice_shelves, pen="0.2p,gray", fill="gainsboro")
            if self.save_fig:
                D_save = Path(self.D_graph,save_dir)
                os.makedirs(D_save, exist_ok=True)
                P_save = os.path.join(D_save, f"sector_{i+1:02d}_grounded_icebergs_{filename_extras}.png")
                fig.savefig(P_save, dpi=300)
                print(f"Saved: {P_save}")
            if self.show_fig:
                fig.show()
            # Optional: fig.savefig(f"sector_{i+1:02d}_grounded_icebergs.png", dpi=300)

    def plot_ice_area(self, area_dict, 
                      ice_type="FI",
                      roll_days=0,
                      tit_str=None,
                      P_png=None,
                      ylim=(0,1000),
                      figsize=(20,12),
                      obs_clim=None,
                      keys_to_plot=None):
        """
        Plot time series of ice area for one or more simulations, optionally with observational climatology.

        Parameters
        ----------
        area_dict : dict
        Dictionary of {sim_name: DataArray} containing 1D time series of area.
        ice_type : str
        One of ["FI", "PI", "SO"]. Used for titles and default y-label.
        roll_days : int
        Rolling mean window size in days.
        tit_str : str
        Title string.
        P_png : Path or str
        Output path for PNG file.
        ylim : tuple
        Y-axis limits.
        figsize : tuple
        Figure size.
        obs_clim : xarray.DataArray, optional
        1D array of 365-day climatology (e.g., AF2020) to repeat and overlay for comparison.
        """
        # Normalize to list if needed
        if keys_to_plot is not None:
            if isinstance(keys_to_plot, str):
                keys_to_plot = [keys_to_plot]
            area_dict = {k: v for k, v in area_dict.items() if k in keys_to_plot}
        df = pd.DataFrame()
        time_array = None
        series_list = []
        for name, da in area_dict.items():
            if isinstance(da, xr.Dataset):
                da = da.to_array().squeeze()
            if roll_days >= 3:
                da_rolled = da.rolling(time=roll_days, center=True, min_periods=1).mean()
            else:
                da_rolled = da
            ser = pd.Series(data=da_rolled.compute().values,
                            index=pd.to_datetime(da["time"].values),
                            name=name)
            series_list.append(ser)
        df = pd.concat(series_list, axis=1)
        df.index.name = "time"
        plt.figure(figsize=figsize, constrained_layout=True)
        # Plot observational climatology as black background line
        if obs_clim is not None:
            doy = df.index.dayofyear
            obs_vals = np.interp(doy, obs_clim["doy"].values, obs_clim.values)
            plt.plot(df.index, obs_vals, label="AF2020 obs", color="black", linewidth=2.0, linestyle="--", alpha=0.7)
        # Plot simulation curves
        for label in df.columns:
            color = self.plot_var_dict.get(label, {}).get("line_clr", None)
            plt.plot(df.index, df[label], label=label, color=color)
        # Year separators
        sep30s = pd.date_range(start=df.index.min(), end=df.index.max(), freq='YE-SEP')
        for dt in sep30s:
            plt.axvline(dt, color='gray', linestyle='--', linewidth=0.8)
        # Vertical dashed lines at year-day 61 (typically March 1st)
        years = range(df.index.year.min(), df.index.year.max() + 1)
        for year in years:
            try:
                day61 = pd.Timestamp(f"{year}-03-01")
            except ValueError:
                # In rare cases, fallback to next day if something goes wrong
                day61 = pd.Timestamp(f"{year}-03-02")
            if df.index.min() <= day61 <= df.index.max():
                plt.axvline(day61, color='gray', linestyle='--', linewidth=0.8, alpha=0.6)
        # Yearly min/max markers and annotations
        years = range(df.index.year.min(), df.index.year.max() + 1)
        for year in years:
            yr_mask = (df.index.year == year)
            for label in df.columns:
                series = df[label][yr_mask].dropna()
                if series.empty:
                    continue
                dt_max, val_max = series.idxmax(), series.max()
                plt.plot(dt_max, val_max, marker='^', color='black')
                plt.text(dt_max, val_max + 10, f'{val_max:.0f}', ha='center', fontsize=8)
                if year not in [df.index.year.min(), df.index.year.max()]:
                    dt_min, val_min = series.idxmin(), series.min()
                    plt.plot(dt_min, val_min, marker='v', color='black')
                    plt.text(dt_min, val_min - 25, f'{val_min:.0f}', ha='center', fontsize=8)
        # --- Labels and output
        plt.title(tit_str)
        if ice_type=="FI":
            plt.ylabel(f"{ice_type} Area (1000-km²)")
        else:
            plt.ylabel(f"{ice_type} Area (1e6-km²)")
        plt.xlabel("Time")
        plt.grid(True)
        plt.ylim(ylim)
        plt.xlim(pd.Timestamp(f"{years.start}-01-01"), pd.Timestamp(f"{years.stop}-01-01"))
        plt.legend()
        plt.tight_layout()
        if self.save_fig:
            plt.savefig(P_png, dpi=100)
            print(f"💾 Saved plot to {P_png}")
        if self.show_fig:
           plt.show()

    def plot_monthly_ice_area_by_year(self, area_dict,
                                      ice_type         = "FI",
                                      roll_days        = 0,
                                      ylim             = (0, 1000),
                                      figsize          = (18, 10),
                                      tit_str          = None,
                                      P_png            = None,
                                      tick_fontsize    = 12,
                                      label_fontsize   = 14,
                                      title_fontsize   = 18,
                                      legend_fontsize  = 12,
                                      plot_annotations = False):
        """
        Plot daily ice area by month (Jan–Dec) for each simulation year (excluding the first).
        For AF2020db_cli: plot as dashed line with max/min markers.
        For each model type (FI_BT, FI_BT_bool, FI_BT_roll):
            - Plot climatology as solid line
            - Shade min/max bounds of each day of year
        If ice_type == "SI", also attempt to plot NSIDC climatology.

        INPUTS:
        area_dict : dict; dictionary of {sim_name: DataArray} containing 1D time series of ice area.
        ice_type  : str; type of ice ("FI", "PI", or "SO" or "SI") for labeling.
        roll_days : int; optional smoothing window (rolling mean in days).
        ylim      : tuple; y-axis limits.
        figsize   : tuple; figure size.
        tit_str   : str; plot title.
        P_png     : Path or str, optional; path to save PNG output.
        """
        plt.figure(figsize=figsize, constrained_layout=True)
        cmap         = plt.get_cmap("tab10")
        month_labels = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
        xticks       = [1,32,62,92,122,152,182,212,242,272,302,332]
        plot_style   = {"FI_BT"      : {"label" : "daily ice-speed mask",
                                        "color" : "black"},
                        "FI_BT_bool" : {"label" : "Binary-days (6 of 7) ice-speed mask",
                                        "color" : "orange"},
                        "FI_BT_roll" : {"label" : "15-day-avg ice-speed mask",
                                        "color" : "green"} }

        # For collecting max/min summary for legend
        sim_annots       = []
        sim_annot_colors = []
        for i, (name, da) in enumerate(area_dict.items()):
            if isinstance(da, xr.Dataset):
                da = da.to_array().squeeze()
            if name == "AF2020db_cli" and ice_type == "FI":
                doy = da["doy"].values
                values = da.values
                plt.plot(doy, values, label="AF2020 Climatology (2000–2018)", linestyle="--", color="blue", linewidth=1.5)
                if plot_annotations:
                    max_doy = doy[np.argmax(values)]
                    min_doy = doy[np.argmin(values)]
                    plt.plot(max_doy, values[np.argmax(values)], "b^")
                    plt.plot(min_doy, values[np.argmin(values)], "bv")
                    sim_annots.append(f"AF2020 Max: {values[np.argmax(values)]:.1f} @ DOY {max_doy}")
                    sim_annots.append(f"AF2020 Min: {values[np.argmin(values)]:.1f} @ DOY {min_doy}")
                    sim_annot_colors.extend(["blue", "blue"])
                continue
            if name == "NSIDC" and ice_type == "SI":
                time = pd.to_datetime(da["time"].values)
                area = da.rolling(time=roll_days, center=True, min_periods=1).mean().values if roll_days >= 1 else da.values
                df = pd.DataFrame({"area": area}, index=time)
                df["doy"] = df.index.dayofyear
                grouped = df.groupby("doy")["area"]
                clim_mean = grouped.mean()
                plt.plot(clim_mean.index, clim_mean.values, label="NSIDC Climatology", linestyle="--", color="black", linewidth=1.5)
                if plot_annotations:
                    max_doy = clim_mean.idxmax()
                    min_doy = clim_mean.idxmin()
                    plt.plot(max_doy, clim_mean[max_doy], "k^")
                    plt.plot(min_doy, clim_mean[min_doy], "kv")
                    sim_annots.append(f"NSIDC Max: {clim_mean[max_doy]:.1f} @ DOY {max_doy}")
                    sim_annots.append(f"NSIDC Min: {clim_mean[min_doy]:.1f} @ DOY {min_doy}")
                    sim_annot_colors.extend(["black", "black"])
                continue
            # simulations
            time       = pd.to_datetime(da["time"].values)
            area       = da.rolling(time=roll_days, center=True, min_periods=1).mean().values if roll_days >= 1 else da.values
            df         = pd.DataFrame({"area": area}, index=time)
            df["doy"]  = df.index.dayofyear
            df["year"] = df.index.year
            # Drop the first year
            first_year = df["year"].min()
            df         = df[df["year"] > first_year]
            # Compute min/max/shaded envelope and climatology
            grouped    = df.groupby("doy")["area"]
            area_min   = grouped.min()
            area_max   = grouped.max()
            area_mean  = grouped.mean()
            style      = plot_style.get(name, {"label": name, "color": cmap(i)})
            label_name = style["label"]
            line_color = style["color"]
            plt.fill_between(area_min.index, area_min, area_max, alpha=0.2, color=line_color, label=f"{label_name} min/max range")
            plt.plot(area_mean.index, area_mean, color=line_color, linewidth=2, label=f"{label_name} climatology")
            if plot_annotations:
                # Markers and legend text for simulation max/min
                max_doy = area_mean.idxmax()
                min_doy = area_mean.idxmin()
                plt.plot(max_doy, area_mean[max_doy], marker="^", color=line_color)
                plt.plot(min_doy, area_mean[min_doy], marker="v", color=line_color)
                sim_annots.append(f"{name} Max: {area_mean[max_doy]:.1f} @ DOY {max_doy}")
                sim_annots.append(f"{name} Min: {area_mean[min_doy]:.1f} @ DOY {min_doy}")
                sim_annot_colors.extend([line_color, line_color])
        plt.xticks(xticks, labels=month_labels, fontsize=tick_fontsize)
        plt.yticks(fontsize=tick_fontsize)
        plt.xlim(1, 366)
        plt.ylim(ylim)
        plt.ylabel(f"Fast Ice Area (1000 km²)" if ice_type == "FI" else f"Sea Ice Area (1e6-km²)", fontsize=label_fontsize)
        plt.xlabel("Month", fontsize=label_fontsize)
        if tit_str:
            plt.title(tit_str, fontsize=title_fontsize)
        plt.grid(axis="y", linestyle="--", alpha=0.5)
        plt.legend(loc="upper left", fontsize=legend_fontsize)
        # Summary annotation box in lower right
        if sim_annots and plot_annotations:
            ax = plt.gca()
            x0, y0 = 0.78, 0.02
            for i, (text, color) in enumerate(zip(sim_annots, sim_annot_colors)):
                ax.text(x0, y0 + i * 0.03, text,
                        transform           = ax.transAxes,
                        fontsize            = 11,
                        verticalalignment   = 'bottom',
                        horizontalalignment = 'center',
                        color               = color)
        if P_png and self.save_fig:
            plt.savefig(P_png, dpi=100)
            print(f"📏 Saved plot to {P_png}")
        if self.show_fig:
            plt.show()

    def plot_timeseries_groupby_month(self, ds, var_name, sim_name=None, dt0_str=None, dtN_str=None):
        self.sim_name = sim_name if sim_name is not None else getattr(self, "sim_name")
        self.dt0_str  = dt0_str  if dt0_str  is not None else getattr(self, 'dt0_str')
        self.dtN_str  = dtN_str  if dtN_str  is not None else self.config.get('dtN_str', '1999-12-31')
        time = pd.to_datetime(ds['t_dim'].values)
        fia = ds['FIA']
        fia_obs = ds['FIA_OBS'].sel(sector='circumpolar')
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
        plt.figure(figsize=(12, 6))
        plt.style.use('ggplot')
        for year in monthly_cycle.columns:
            plt.plot(monthly_cycle.index, monthly_cycle[year], label=f"{year}", linewidth=2)
        plt.plot(monthly_obs.index, monthly_obs.values, label='Observed FIA (AF2020db)', color='black', linestyle='--', linewidth=3)
        plt.title("Monthly Fast Ice Area Cycle (1993–1999)", fontsize=16)
        plt.xlabel("Month", fontsize=14)
        plt.ylabel(f"Fast Ice Area ({fia.attrs.get('units', 'unknown')})", fontsize=14)
        plt.xticks(ticks=range(1,13), labels=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'])
        plt.grid(True, linestyle='--', alpha=0.5)
        plt.legend(loc='upper right', fontsize=10, title='Year')
        plt.tight_layout()
        plt.show()
        plt.savefig(f"/g/data/gv90/da1339/GRAPHICAL/timeseries/{sim_name}_FIA_grouped_by_month.png")

    def plot_facet_regions(self, 
                           da_dict         = None,
                           var_name        = 'FIP',
                           tit_str         = None,
                           filename_extras = None,
                           report_stats    = False,
                           time_coord_name = 'time',
                           lat_coord_name  = 'lat',
                           lon_coord_name  = 'lon',
                           figure_size     = ("60c","90c"),
                           panel_margins   = [".1c", ".01c"],
                           cbar_position   = "JBC+w10c/1c+h",
                           **kwargs):
        """
        """
        self.load_ice_shelves()
        var_dict     = self.plot_var_dict[var_name]
        cmap         = kwargs.get("cmap"        , var_dict.get("cmap", "viridis"))
        series       = kwargs.get("series"      , var_dict.get("series", [0, 1]))
        cmap_reverse = kwargs.get("cmap_reverse", var_dict.get("reverse", False))
        units        = kwargs.get("units"       , var_dict.get("units", ""))
        cbar_str     = kwargs.get("cbar_str"    , var_dict.get("name", var_name))
        land_color   = kwargs.get("land_color"  , self.pygmt_dict.get("land_color", "seashell"))
        water_color  = kwargs.get("water_color" , self.pygmt_dict.get("water_color", "white"))
        sq_size_var  = kwargs.get("sq_size_var" , self.pygmt_dict.get("sq_size_var", 0.15))
        sq_size_GI   = kwargs.get("sq_size_GI"  , self.pygmt_dict.get("sq_size_GI", 0.05))
        cbar_frame   = self.create_cbar_frame(series, cbar_str, units=units)
        self.define_png_plot_path( var_name=var_name, qualifiers = filename_extras)
        if os.path.exists(self.P_png) and not self.ow_fig:
            print(f"figure already exists and not overwriting ... skipping {self.P_png}")
            self.P_png = None
            return
        t0    = time.time()
        fig   = pygmt.Figure()
        pygmt.config(FORMAT_GEO_MAP="ddd.x", FONT_TITLE="16p,Courier-Bold", FONT_ANNOT_PRIMARY="14p,Helvetica")
        pygmt.makecpt(cmap=cmap, reverse=cmap_reverse, series=series)
        #n_panel_rows = len(self.reg_dict.keys())+1
        n_panel_cols = len(da_dict.keys())+1
        sim_names    = list(da_dict.keys())
        region_items = list(self.reg_dict.items())
        region_name  = "EIO"
        region_config = self.reg_dict[region_name]
        with fig.subplot( nrows=1, ncols=n_panel_cols, figsize=figure_size, title=tit_str, margins=panel_margins, sharex="b", sharey="l"):
            for i, sim_name in enumerate(sim_names):
                print(f"working on simulation: {sim_name}")
                da_full  = da_dict[sim_name]
                lat_full = da_full[lat_coord_name].values
                lon_full = da_full[lon_coord_name].values
                da_np    = da_full.values if hasattr(da_full, "values") else da_full
                mask_na  = ~np.isnan(da_np)
                lat_plt  = lat_full.ravel()[mask_na.ravel()]
                lon_plt  = lon_full.ravel()[mask_na.ravel()]
                da_plt   = da_np.ravel()[mask_na.ravel()]
                GI_proc  = GroundedIcebergProcessor(sim_name=sim_name)
                GI_proc.load_AFIM_GI()
                use_gi   = GI_proc.use_gi
                #for j, (region_name, region_config) in enumerate(region_items):
                print(f"working on region: {region_name}")
                region_extent              = region_config['plt_ext']   
                region_MC                  = self.get_meridian_center_from_geographic_extent(region_extent)
                projection                 = f"S{region_MC}/-90/?"#20c
                    #with fig.set_panel(panel=[j, i]):
                    #    if j==0:
                frame = ["af",f"+t{sim_name}"]
                    #    else:
                    #        frame = ["af"]
                fig.basemap(projection=projection, region=region_extent, frame=frame, panel=i)
                fig.coast(land=land_color, water=water_color, panel=i)
                fig.plot(x=lon_plt, y=lat_plt, cmap=True, fill=da_plt, style=f"s{sq_size_var}c", panel=i)
                fig.plot(data=self.ice_shelves, pen="0.2p,gray", fill="lightgray", panel=i)
                if use_gi:
                    fig.plot(x=GI_proc.GI_lon_cells, y=GI_proc.GI_lat_cells, fill='red', style=f"c{sq_size_GI}c", panel=i)
                        # if report_stats:
                        #   fig.text(x=self.text_loc[0], y=self.text_loc[1], text=self.text_str, font=self.text_font, justify=self.text_justify, fill=self.text_fill, pen=self.text_pen, no_clip=True)
        fig.colorbar(position=cbar_position, frame=cbar_frame)
        if self.show_fig:
            fig.show()
        if self.save_fig:
            fig.savefig(self.P_png)
            print(f"figure saved: {self.P_png}")
            self.P_png = None
        print(f"time taken {time.time()-t0:0.2f} seconds\n")
        pygmt.clib.Session.__exit__

