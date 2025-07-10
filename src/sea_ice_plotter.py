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
        Load Antarctic ice shelf polygons from a shapefile into the processor.

        This sets the `self.antarctic_ice_shelves` attribute for optional overlay
        during plotting. Geometries are cleaned, reprojected, and buffered.
        """
        gdf                     = gpd.read_file(self.config['pygmt_dict']['P_coast_shape'])
        shelves                 = gdf[gdf['POLY_TYPE'] == 'S']
        shelves                 = shelves[~shelves.geometry.is_empty & shelves.geometry.notnull()]
        shelves                 = shelves.to_crs("EPSG:4326")
        shelves.geometry        = shelves.geometry.buffer(0)
        return shelves.geometry

    def create_IBCSO_bath(self):
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
        return xr.open_dataset(self.config['pygmt_dict']["P_IBCSO_bath"]).bath

    def prepare_data_for_pygmt_plot(self, da, lon_coord_name=None, lat_coord_name=None, diff_plot=False):
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
        GI_loc_dict        = {}
        kmt_mod            = xr.open_dataset(self.P_KMT_mod).isel(nj=self.hemisphere_dict['nj_slice']).kmt.data
        kmt_org            = xr.open_dataset(self.P_KMT_org).isel(nj=self.hemisphere_dict['nj_slice']).kmt.data
        GI_mask            = (kmt_org == 1) & (kmt_mod == 0)
        GI_loc_dict['lon'] = data_dict['lon2d'][GI_mask].ravel()
        GI_loc_dict['lat'] = data_dict['lat2d'][GI_mask].ravel()
        return GI_loc_dict

    def create_cbar_frame(self, series, label, units=None, extend_cbar=False, max_ann_steps=10):
        """
        Given a data range and label, create GMT-style colorbar annotation string.
        
        Parameters
        ----------
        series : list or tuple of [vmin, vmax]
        label  : str
            Label for the colorbar
        units  : str, optional
            Units to append on secondary axis
        extend_cbar : bool, optional
            Whether to append extension arrows to the colorbar
        max_ann_steps : int
            Target number of annotated steps on the colorbar

        Returns
        -------
        str or list[str]
            GMT-style colorbar annotation frame(s)
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
                              FIA_ylim       = (100, 1000),
                              roll_days      = 0,
                              lon_coord_name = None,
                              lat_coord_name = None,
                              cmap           = None,
                              series         = None,
                              cbar_frame     = None,
                              overwrite_fig  = None,
                              show_fig       = None):
        """
        Composite plot showing monthly climatology of FIA and east/west persistence map.
        """
        sim_name   = sim_name      if sim_name      is not None else self.sim_name
        show_fig   = show_fig      if show_fig      is not None else self.show_fig
        ow_fig     = overwrite_fig if overwrite_fig is not None else self.ow_fig
        cmap       = cmap          if cmap          is not None else self.pygmt_dict.get("FIP_CPT")
        series     = series        if series        is not None else [0.01, 1.0, 0.01]
        cbar_frame = cbar_frame    if cbar_frame    is not None else ["x+lFast Ice Persistence", "y+l1/100"]
        # load fast ice observation climatology
        AF2020_CSV = self.load_AF2020_FIA_summary( start=self.dt0_str, end=self.dtN_str )
        obs_FIA    = AF2020_CSV['FIA_clim_repeat'].sel(region='circumpolar')
        # ----- Create and save x-axis annotation file for months -----
        xticks       = [1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 360]
        month_labels = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec", "Dec"]
        xannot_path  = Path("xannots.txt")
        xannot_lines = [f"{tick}\tig\t{label}" for tick, label in zip(xticks, month_labels)]
        xannot_path.write_text("\n".join(xannot_lines) + "\n")
        # FIA figure
        default_line_colors = ["orange", "green", "blue", "red", "magenta", "cyan"]
        fig = pygmt.Figure()
        pygmt.config(FORMAT_GEO_MAP="ddd.x", MAP_GRID_PEN="0p,white")
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
                da = da.load()  # explicitly load once (Dask-safe)
                # Convert time to day-of-year
                df = pd.DataFrame({
                    "time": pd.to_datetime(da["time"].values),
                    "area": da.values
                })
                df["doy"]  = df["time"].dt.dayofyear
                df["year"] = df["time"].dt.year
                # Drop first year (optional)
                df = df[df["year"] > df["year"].min()]
                grouped = df.groupby("doy")["area"]
                area_min  = grouped.min()
                area_max  = grouped.max()
                area_mean = grouped.mean()
                # Observation climatology
                if i == 0 and not diff_plot:
                    obs_df = pd.DataFrame({
                        "time": pd.to_datetime(obs_FIA["time"].values),
                        "area": obs_FIA.values
                    })
                    obs_df["doy"] = obs_df["time"].dt.dayofyear
                    fig.plot(x=obs_df["doy"], y=obs_df["area"], pen="1.5p,blue,--", label="AF2020 Climatology (2000–2018)")
                # Labels/colors
                if name in self.plot_ice_area_dict:
                    leg_lab    = self.plot_ice_area_dict[name]["label"].format(min_win=self.bool_min_days, win=self.bool_window)
                    line_color = self.plot_ice_area_dict[name]["color"]
                else:
                    leg_lab    = name
                    line_color = default_line_colors[i % len(default_line_colors)]
                # Plot either raw or difference
                if diff_plot:
                    obs_interp = np.interp(area_mean.index, obs_df["doy"], obs_df["area"])
                    area_mean_diff = obs_interp - area_mean.values
                    fig.plot(x=area_mean.index, y=area_mean_diff, pen=f"2p,{line_color}", label=f"{leg_lab} diff")
                else:
                    fig.plot(
                        x=np.concatenate([area_min.index, area_max.index[::-1]]),
                        y=np.concatenate([area_min.values, area_max.values[::-1]]),
                        fill=f"{line_color}@70", close=True, transparency=60
                    )
                    fig.plot(x=area_mean.index, y=area_mean.values, pen=f"2p,{line_color}", label=f"{leg_lab} climatology")
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
        Composite plot showing monthly climatology of FIA and east/west persistence map.
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

    def generate_regional_annotation_stats(self, da, region, lon_coord_name, lat_coord_name):
        """
        Add statistics (mean, std, min/max) and their locations on the map.
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
                        roll_days : int   = None,
                        tit_str   : str   = None,
                        ylim      : tuple = None,
                        ylabel    : str   = "Fast Ice Area (1000-km²)",
                        ytick_inc : int   = 100,
                        xlabel    : str   = "Date",
                        xtick_inc : int   = 1,
                        P_png     : str   = None,
                        fig_width : str   = "15c",
                        fig_height: str   = "5c",
                        legend_pos: str   = "JTL+jTL+o0.2c",
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
        show_fig = self.show_fig if show_fig is None else show_fig
        if keys2plot is not None:
            if isinstance(keys2plot, str):
                keys2plot = [keys2plot]
            ts_dict = {k: v for k, v in ts_dict.items() if k in keys2plot}
        dfs = []
        for name, da in ts_dict.items():
            if isinstance(da, xr.Dataset):
                da = da.to_array().squeeze()
            if time_coord not in da.coords:
                raise ValueError(f"Missing time coordinate '{time_coord}' in variable '{name}'")
            time        = pd.to_datetime(da[time_coord].values)
            values      = da.rolling({time_coord: roll_days}, center=True, min_periods=1).mean().values if roll_days and roll_days >= 3 else da.values
            df          = pd.DataFrame({"time": time, "value": values})
            df["label"] = name
            dfs.append(df)
        df_all     = pd.concat(dfs, ignore_index=True).dropna()
        all_times  = df_all["time"]
        tmin, tmax = all_times.min(), all_times.max()
        if ylim is None:
            ymin = df_all["value"].min()
            ymax = df_all["value"].max()
            pad  = 0.05 * (ymax - ymin)
            ylim = [ymin - pad, ymax + pad]
        region     = [tmin.strftime("%Y-%m-%dT00:00:00"), tmax.strftime("%Y-%m-%dT00:00:00"), ylim[0], ylim[1]]
        projection = f"X{fig_width}/{fig_height}"
        frame      = [f"xaf{xtick_inc}y+l\"{xlabel}\"", f"yaf{ytick_inc}+l\"{ylabel}\"", "WSne"]
        if tit_str:
            frame.append(f"+t\"{tit_str}\"")
        fig = pygmt.Figure()
        pygmt.config(FONT_TITLE              = "18p,Helvetica-Bold",
                     FONT_ANNOT_PRIMARY      = "14p,Helvetica",
                     FORMAT_DATE_MAP         = "o",
                     FORMAT_TIME_PRIMARY_MAP = "full")
        fig.basemap(region=region, projection=projection, frame=frame)
        colors = ["red", "blue", "green", "purple", "orange"]
        for i, (label, _) in enumerate(ts_dict.items()):
            df = df_all[df_all["label"] == label]
            color = colors[i]#self.plot_var_dict.get(label, {}).get("line_clr", f"C{i}")
            legend_label = sim_name#self.plot_var_dict.get(label, {}).get("label", label)
            fig.plot(x=df["time"], y=df["value"], pen=f"{pen_weight},{color}", label=legend_label)
        fig.legend(position=legend_pos, box=legend_box)
        if P_png and self.save_fig:
            fig.savefig(P_png)
            print(f"📏 Saved plot to {P_png}")
        if show_fig:
            fig.show()
        return fig

    def plot_timeseries_matplotlib(self, ts_dict,
                                    roll_days: int = None,
                                    tit_str: str = None,
                                    ylim: tuple = None,
                                    ylabel: str = "Fast Ice Area (1000-km²)",
                                    ytick_inc: int = 100,
                                    xlabel: str = "Date",
                                    xtick_inc: int = 1,
                                    P_png: str = None,
                                    fig_width: str = "15c",
                                    fig_height: str = "5c",
                                    legend_pos: str = "best",  # Matplotlib uses standard positioning
                                    legend_box: str = "upper left",
                                    pen_weight: int = 2,
                                    time_coord: str = "time",
                                    keys2plot: list = None,
                                    show_fig: bool = None):
        """
        Plot one or more time series using Matplotlib.

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
            Interval for x-axis ticks.
        P_png : str, optional
            Path to save PNG output. If None, figure is not saved.
        fig_width : str
            Width of the figure (e.g., "15c").
        fig_height : str
            Height of the figure (e.g., "5c").
        legend_pos : str
            Matplotlib legend position string (e.g., "best").
        legend_box : str
            Matplotlib legend box style string.
        pen_weight : int
            Line thickness for the time series (default 2).
        time_coord : str
            Name of the time coordinate in each xarray object.
        keys2plot : list, optional
            Subset of keys in ts_dict to plot. If None, all keys are used.
        show_fig : bool, optional
            If True, show the figure interactively. Defaults to self.show_fig.

        Returns
        -------
        fig : matplotlib.figure.Figure
            The generated Matplotlib figure.
        """
        show_fig = self.show_fig if show_fig is None else show_fig
        if keys2plot is not None:
            if isinstance(keys2plot, str):
                keys2plot = [keys2plot]
            ts_dict = {k: v for k, v in ts_dict.items() if k in keys2plot}

        # Prepare the data for plotting
        dfs = []
        for name, da in ts_dict.items():
            if isinstance(da, xr.Dataset):
                da = da.to_array().squeeze()
            if time_coord not in da.coords:
                raise ValueError(f"Missing time coordinate '{time_coord}' in variable '{name}'")
            
            # Convert to Pandas DataFrame
            time = pd.to_datetime(da[time_coord].values)
            values = da.values
            if roll_days and roll_days >= 3:
                values = pd.Series(values).rolling(window=roll_days, center=True, min_periods=1).mean().values
            
            df = pd.DataFrame({"time": time, "value": values})
            df["label"] = name
            dfs.append(df)
        
        # Combine all DataFrames into one
        df_all = pd.concat(dfs, ignore_index=True).dropna()

        # Set axis limits
        all_times = df_all["time"]
        tmin, tmax = all_times.min(), all_times.max()
        if ylim is None:
            ymin = df_all["value"].min()
            ymax = df_all["value"].max()
            pad = 0.05 * (ymax - ymin)
            ylim = [ymin - pad, ymax + pad]

        # Set up the plot
        fig, ax = plt.subplots(figsize=(int(fig_width[:-1]), int(fig_height[:-1])))  # Convert to inches
        ax.set_xlim([tmin, tmax])
        ax.set_ylim(ylim)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        
        # Add title if provided
        if tit_str:
            ax.set_title(tit_str)
        
        # Plot the data
        for label, _ in ts_dict.items():
            df = df_all[df_all["label"] == label]
            color = self.plot_var_dict.get(label, {}).get("line_clr", f"C{list(ts_dict.keys()).index(label)}")
            legend_label = self.plot_var_dict.get(label, {}).get("label", label)
            ax.plot(df["time"], df["value"], label=legend_label, lw=pen_weight, color=color)

        # Set x-axis ticks
        ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True, prune="both"))
        ax.xaxis.set_major_formatter(plt.matplotlib.dates.DateFormatter('%Y-%m-%d'))

        # Customize the y-axis ticks
        ax.yaxis.set_major_locator(plt.MultipleLocator(ytick_inc))

        # Add legend
        ax.legend(loc=legend_pos, frameon=True, facecolor="white", edgecolor="black")

        # Show or save the plot
        if P_png and self.save_fig:
            fig.savefig(P_png)
            print(f"📏 Saved plot to {P_png}")
        
        if show_fig:
            plt.show()
        
        return fig


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