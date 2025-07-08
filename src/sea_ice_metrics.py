import gc
import xarray      as xr
import numpy       as np
import pandas      as pd
from pathlib       import Path
from scipy.ndimage import binary_dilation, distance_transform_edt

class SeaIceMetrics:

    def __init__(self, **kwargs):
        return

    def sea_ice_metrics_wrapper(self,
                                sim_name        = None,
                                ice_type        = None,
                                dt0_str         = None,
                                dtN_str         = None,
                                ispd_thresh     = None,
                                D_out           = None,
                                overwrite_zarr  = False,
                                overwrite_png   = False,
                                smooth_FIA_days = 15,
                                drop_first_year = True):
        """
        """
        sim_name           = sim_name if sim_name is not None else self.sim_name
        dt0_str            = dt0_str  if dt0_str  is not None else "1994-01-01"
        dtN_str            = dtN_str  if dtN_str  is not None else "1999-12-31"
        ispd_thresh        = ispd_thresh or self.ispd_thresh
        ispd_thresh_str    = f"{ispd_thresh:.1e}".replace("e-0", "e-")
        ice_type           = ice_type or self.ice_type
        dt_range_str       = f"{dt0_str[:4]}-{dtN_str[:4]}"
        FIA_dict           = {}
        AF2020_CSV         = self.load_AF2020_FIA_summary( start=dt0_str, end=dtN_str )
        FIA_dict["AF2020"] = AF2020_CSV['FIA_clim_repeat'].sel(region='circumpolar')
        ice_types          = [ice_type, f"{ice_type}_roll", f"{ice_type}_bool"]
        for i_type in ice_types:
            P_METS = Path(self.D_metrics, f"{i_type}_mets.zarr")
            P_sum  = Path(self.D_metrics, f"{i_type}_summary.csv")
            if P_METS.exists() and not self.overwrite_zarr_group:
                self.logger.info(f"{P_METS} exists and not overwriting--loading")
                METS             = xr.open_zarr(P_METS)
                FIA_dict[i_type] = METS['FIA']
            else:
                self.logger.info(f"{P_METS} does NOT exists and/or overwriting--computing")
                if hasattr(i_type, "endswith"):
                    roll = i_type.endswith("_roll")
                else:
                    roll = False
                # important to know that "FI" (below) is either daily (``raw``; i.e. classified without any other criteria
                # other than ice speed threshold) or daily-rolling-mean (``rolled``; i.e. classified after ice speed has had
                # a rolling-mean applied to it) ... the load of either is controlled by the True/False "rolling" option 
                FI, CICE_SO = self.load_processed_cice(ispd_thresh = ispd_thresh,
                                                       ice_type    = ice_type,
                                                       dt0_str     = dt0_str,
                                                       dtN_str     = dtN_str,
                                                       zarr_CICE   = True,
                                                       rolling     = roll,
                                                       slice_hem   = True)
                if i_type==f"{ice_type}_bool":
                    bool_mask          = self.boolean_fast_ice(FI['FI_mask'])
                    FI_bool            = CICE_SO.where(bool_mask)
                    FI_bool["FI_mask"] = bool_mask
                    FI                 = FI_bool
                METS = self.compute_sea_ice_metrics( FI, FIA_dict["AF2020"], P_METS )
            FIA_dict[i_type] = METS['FIA']
        P_png = Path(self.D_graph, sim_name, f"FIA_FIP_{sim_name}_{ispd_thresh_str}_{dt_range_str}.png")
        self.plot_FIA_FIP_faceted(FIA_dict, METS['FIP'], P_png=P_png, plot_GI=True if self.use_gi else False)

    def compute_sea_ice_metrics(self, FI_sim, FIA_obs, P_METS=None):
        """
        Compute sea ice metrics from a processed fast ice simulation and compare against observations.

        Parameters
        ----------
        FI_sim : xr.Dataset
            Processed fast ice dataset containing 'aice', 'hi', 'tarea', and 'FI_mask'.
        FIA_obs : xr.DataArray
            Observational fast ice area time series (e.g., AF2020 climatology repeated).
        P_METS : Path or str, optional
            Output path to store Zarr archive of computed metrics.

        Returns
        -------
        xr.Dataset
            Dataset containing computed FIA, FIP, and additional seasonal and spatial metrics.
        """
        METS = {}
        # --- 1D Time Series Metrics ---
        FIA = self.compute_ice_area(FI_sim['aice'], FI_sim['tarea'])
        FIV = self.compute_ice_volume(FI_sim['aice'], FI_sim['hi'], FI_sim['tarea'])
        FIP = self.compute_variable_aggregate(FI_sim['aice'])
        # Compute immediately before use in downstream stat methods
        FIA_comp = FIA.compute()
        del FIA
        gc.collect()
        FIV_comp = FIV.compute()
        del FIV
        gc.collect()
        FIP_comp = FIP.compute()
        del FIP
        gc.collect()
        METS["FIA"] = FIA_comp
        METS["FIV"] = FIV_comp
        METS["FIP"] = FIP_comp
        # --- Seasonal Statistics (1D time series) ---
        try:
            FIA_seasonal = self.compute_seasonal_statistics(FIA_comp)
            FIV_seasonal = self.compute_seasonal_statistics(FIV_comp)
        except Exception as e:
            self.logger.warning(f"⚠️ Seasonal statistics failed: {e}")
            FIA_seasonal = {}
            FIV_seasonal = {}
        # --- Spatial Statistics from FIP (2D map) ---
        try:
            FIP_stats = self.compute_fip_spatial_stats(FIP_comp)
        except Exception as e:
            self.logger.warning(f"⚠️ FIP spatial statistics failed: {e}")
            FIP_stats = {}
        # --- Fast Ice Distance Metrics ---
        try:
            FI_dist = self.compute_fast_ice_distance_extent(FI_sim['FI_mask'])
        except Exception as e:
            self.logger.warning(f"⚠️ Distance extent computation failed: {e}")
            FI_dist = {}
        # --- Skill Statistics vs Observations ---
        try:
            FIA_skill = self.compute_skill_statistics(FIA_comp, FIA_obs)
        except Exception as e:
            self.logger.warning(f"⚠️ Skill statistics failed: {e}")
            FIA_skill = {}
        # --- Build Dataset from METS ---
        DS_METS = xr.Dataset()
        for k, v in METS.items():
            if isinstance(v, xr.DataArray):
                DS_METS[k] = v
            elif np.ndim(v) == 2:
                DS_METS[k] = xr.DataArray(v, dims=('nj', 'ni'))
            elif np.ndim(v) == 1:
                DS_METS[k] = xr.DataArray(v, dims=('time',))
            else:
                DS_METS[k] = xr.DataArray(v, dims=())
        # --- Merge All Metrics into Dataset ---
        summary = {**{f"FIA_{k}": v for k, v in FIA_seasonal.items()},
                   **{f"FIV_{k}": v for k, v in FIV_seasonal.items()},
                   **FIP_stats,
                   **FI_dist,
                   **FIA_skill}
        for k, v in summary.items():
            DS_METS[k] = xr.DataArray(v, dims=())  # All summary stats are scalars
        # --- Save to Zarr (Optional) ---
        if P_METS is not None:
            DS_METS.to_zarr(P_METS, mode="w", consolidated=True)
            self.logger.info(f"📊 Metrics written to {P_METS}")
        return DS_METS

    def compute_ice_area(self, SIC, GC_area, ice_area_scale=None, spatial_dim_names=None, add_grounded_iceberg_area=None, grounded_iceberg_area=None):
        """
        Compute the total sea ice area (IA) by integrating sea ice concentration (SIC) over the grid area, and optionally 
        including the grounded iceberg area (GI_total_area).

        This method computes the sea ice area by multiplying the sea ice concentration (SIC) by the grid cell area (GC_area),
        summing the result over the specified spatial dimensions, and optionally adding the grounded iceberg area. The final
        result is then scaled by an optional ice area scale factor.

        Parameters:
        -----------
        SIC : xarray.DataArray
            Sea ice concentration field with dimensions (nj, ni) representing the concentration of sea ice at each grid point.
            
        GC_area : xarray.DataArray
            Grid cell area with the same dimensions (nj, ni) as SIC, representing the area of each grid cell in square meters.
            
        ice_area_scale : float, optional
            A scale factor for the sea ice area calculation. If not provided, the default scale is used from the `FIC_scale` attribute.
            
        spatial_dim_names : list of str, optional
            The dimension names over which to sum the sea ice area (typically latitude and longitude). If not provided, the default spatial dimensions 
            are used from the `CICE_dict` attribute.
            
        add_grounded_iceberg_area : bool, optional
            A flag indicating whether to include the grounded iceberg area in the ice area calculation. Defaults to the class attribute `use_gi`.
            
        grounded_iceberg_area : float, optional
            The grounded iceberg area in square meters to be added to the sea ice area. If not provided, the grounded iceberg area is computed using 
            `compute_grounded_iceberg_area()`.

        Returns:
        --------
        IA : xarray.DataArray
            The total sea ice area, including the optional grounded iceberg area, with the same spatial dimensions as SIC and GC_area.
            The value is returned after summing over the spatial dimensions and scaling by the provided `ice_area_scale`.

        Notes:
        ------
        - The grounded iceberg area is added to the sea ice area calculation if `add_grounded_iceberg_area` is `True`. If no grounded iceberg area 
        is provided, the method will compute it using the `compute_grounded_iceberg_area()` method.
        - The sea ice area is computed as the sum of SIC * GC_area across the specified spatial dimensions, and then scaled by `ice_area_scale`.
        """
        add_grounded_iceberg_area = add_grounded_iceberg_area if add_grounded_iceberg_area is not None else self.use_gi
        if add_grounded_iceberg_area:
            if grounded_iceberg_area is not None:
                GI_total_area = grounded_iceberg_area
            else:
                GI_total_area = self.compute_grounded_iceberg_area()
        else:
            GI_total_area = 0
        spatial_dim_names = spatial_dim_names if spatial_dim_names is not None else self.CICE_dict['spatial_dims']
        self.logger.info(f"{GI_total_area:0.2f} m^2 total circumpolar grounded iceberg area for {self.sim_name}")
        ice_area_scale = ice_area_scale if ice_area_scale is not None else self.FIC_scale
        self.logger.info(f"🧮 Spatially-integrating the product of sea ice concentrations and grid cell areas")
        IA = ((SIC * GC_area).sum(dim=spatial_dim_names))
        IA = IA + GI_total_area
        IA = IA/ice_area_scale
        return IA

    def compute_ice_extent(self, SIC, GC_area, ice_extent_threshold=0.15, spatial_dim_names=None, 
                        add_grounded_iceberg_area=None, grounded_iceberg_area=None):
        """
        Compute sea ice extent (SIE) by summing the grid cell areas where the sea ice concentration exceeds the threshold.
        Parameters:
            SIC                   : Sea ice concentration (array)
            GC_area               : Grid cell area (array)
            ice_extent_threshold  : Threshold for sea ice extent (default 0.15)
            spatial_dim_names     : List of spatial dimension names for summing (default to class attribute)
            add_grounded_iceberg_area : Whether to include grounded iceberg area in the computation
            grounded_iceberg_area : Grounded iceberg area (optional)
        Returns:
            SIE                   : Sea ice extent
        """
        
        add_grounded_iceberg_area = add_grounded_iceberg_area if add_grounded_iceberg_area is not None else self.use_gi
        if add_grounded_iceberg_area:
            if grounded_iceberg_area is not None:
                GI_total_area = grounded_iceberg_area
            else:
                GI_total_area = self.compute_grounded_iceberg_area()
        else:
            GI_total_area = 0
            
        spatial_dim_names = spatial_dim_names if spatial_dim_names is not None else self.CICE_dict['spatial_dims']
        self.logger.info(f"{GI_total_area:0.2f} m^2 total circumpolar grounded iceberg area for {self.sim_name}")

        # Apply the threshold for sea ice concentration (e.g., SIC >= 0.15 for SIE)
        SIC_thresholded = SIC.where(SIC >= ice_extent_threshold, 0)

        # Compute sea ice extent by summing the grid cell areas where SIC > threshold
        self.logger.info(f"🧮 Spatially-integrating sea ice concentration exceeding threshold ({ice_extent_threshold * 100}%)")
        SIE = (SIC_thresholded * GC_area).sum(dim=spatial_dim_names)

        # Optionally add grounded iceberg area
        SIE = SIE + GI_total_area

        return SIE


    def compute_ice_volume(self, SIC, HI, GC_area, ice_volume_scale=1e12, spatial_dim_names=None):
        """
        Compute total sea ice volume in m³, optionally scaled for output.

        Parameters:
        - SIC: Sea ice concentration [unitless]
        - HI: Sea ice thickness [m]
        - GC_area: Grid cell area [m²]
        - ice_volume_scale: Scale factor (e.g. 1e12 for 1000 km³)

        Returns:
        - Scaled total sea ice volume
        """
        spatial_dim_names = spatial_dim_names if spatial_dim_names is not None else self.CICE_dict['spatial_dims']
        GI_total_area = self.compute_grounded_iceberg_area() if self.use_gi else 0
        self.logger.info(f"{GI_total_area:0.2f} m² total grounded iceberg area for {self.sim_name}")
        self.logger.info("🧮 Computing ice volume: sum(SIC × HI × area)")
        IV = (SIC * HI * GC_area).sum(dim=spatial_dim_names)
        IV += GI_total_area
        IV /= ice_volume_scale
        return IV

    def compute_variable_aggregate(self, da, time_coord_name='time'):
        return da.sum(dim=time_coord_name) / da[time_coord_name].sizes.get(time_coord_name, 1)

    # SEASONAL METRICS
    def detect_onset_via_derivative(self, da, window=15, polyorder=2, min_onset_doy=50):
        """
        Detect onset of ice growth using Savitzky-Golay derivative filter.
        
        Parameters
        ----------
        da : xr.DataArray
            Time series with daily time axis (single year).
        window : int
            Smoothing window for Savitzky-Golay filter (must be odd).
        polyorder : int
            Polynomial order for filter.
        min_onset_doy : int
            DOY to start searching for growth onset.
        
        Returns
        -------
        int or None
            Estimated DOY of onset, or None if not found.
        """
        from scipy.signal import savgol_filter
        if da.time.size < window:
            return None
        y     = da.values
        dy_dt = savgol_filter(y, window_length=window, polyorder=polyorder, deriv=1)
        doy   = da["time"].dt.dayofyear.values
        mask  = doy >= min_onset_doy
        growth_start_candidates = np.where((dy_dt > 0) & mask)[0]
        if growth_start_candidates.size == 0:
            return None
        # Optional: Find first continuous upward segment
        onset_idx = growth_start_candidates[0]
        return int(doy[onset_idx])

    def compute_seasonal_duration_derivative(da, window=15, polyorder=2, min_onset_doy=50, max_retreat_doy=330):
        """
        Compute ice season duration using first derivative of smoothed time series.

        Seasonal Duration = Retreat DOY − Onset DOY
        Where: Onset DOY is the beginning of sustained growth, and Retreat DOY is the end of sustained presence (start of ablation),
        
        Parameters
        ----------
        da : xr.DataArray
            Time series for a single year (daily resolution).
        window : int
            Smoothing window for Savitzky-Golay filter.
        polyorder : int
            Polynomial order for smoothing.
        min_onset_doy : int
            Minimum DOY to search for onset.
        max_retreat_doy : int
            Maximum DOY to search for retreat.

        Returns
        -------
        (onset_doy, retreat_doy, duration_days)
            Tuple of day-of-year values (int) and duration (int), or (None, None, None)
        """
        from scipy.signal import savgol_filter
        if da.time.size < window:
            return None, None, None
        y   = da.values
        doy = da["time"].dt.dayofyear.values
        dy_dt = savgol_filter(y, window_length=window, polyorder=polyorder, deriv=1)
        # --- Onset detection (first sustained positive trend) ---
        onset_candidates = np.where((dy_dt > 0) & (doy >= min_onset_doy))[0]
        if onset_candidates.size == 0:
            return None, None, None
        onset_idx = onset_candidates[0]
        onset_doy = doy[onset_idx]
        # --- Retreat detection (first sustained negative trend after peak) ---
        peak_idx = np.argmax(y)
        retreat_candidates = np.where((dy_dt < 0) & (np.arange(len(y)) > peak_idx) & (doy <= max_retreat_doy))[0]
        if retreat_candidates.size == 0:
            return onset_doy, None, None
        retreat_idx = retreat_candidates[-1]
        retreat_doy = doy[retreat_idx]
        # Handle wraparound (e.g., retreat in next year)
        duration = retreat_doy - onset_doy if retreat_doy > onset_doy else (365 - onset_doy + retreat_doy)
        return onset_doy, retreat_doy, duration

    def compute_retreat_via_inflection(self, da, window=15, polyorder=2, min_retreat_doy=240):
        from scipy.signal import savgol_filter
        y           = da.values
        doy         = da["time"].dt.dayofyear.values
        d2y_dt2     = savgol_filter(y, window_length=window, polyorder=polyorder, deriv=2)
        inflections = np.where((d2y_dt2 < 0) & (doy >= min_retreat_doy))[0]
        if inflections.size == 0:
            return None
        return doy[inflections[-1]]

    def compute_seasonal_statistics(self, da, 
                                    growth_range        = (71 , 273),
                                    retreat_early_range = (273, 330),
                                    retreat_late_range1 = (331, 365),
                                    retreat_late_range2 = (1  , 71)  ):
        """
        Generalized method to compute seasonal growth and retreat statistics for any sea ice metric (e.g., SIA, FIA, SIV, FIV).

        Parameters
        ----------
        da                  : xr.DataArray;  Time series of a sea ice metric with a 'time' dimension.
        var_name            : str, optional; Label for result entry (e.g., 'SIA', 'FIV').
        growth_range        : tuple of int;  Day-of-year range for ice growth period.
        retreat_early_range : tuple of int;  DOY range for early melt (e.g., melt onset to late spring).
        retreat_late_range1 : tuple of int;  DOY range for late melt in previous year.
        retreat_late_range2 : tuple of int;  DOY range for late melt in next year.

        Returns
        -------
        dict; Summary statistics: mean/std growth and retreat rates.
        """
        da                  = da.sel(time=~da["time"].dt.is_leap_year)
        da                  = da.sel(time=da["time"].dt.year > da["time"].dt.year.min())
        maximums            = []
        minimums            = []
        growth_rates        = []
        retreat_rates       = []
        retreat_early_rates = []
        retreat_late_rates  = []
        min_doys            = []
        max_doys            = []
        onset_doys          = []
        duration_days       = []
        for year in np.unique(da.time.dt.year):
            this_year = da.sel(time=da.time.dt.year == year).compute()
            maximums.append(this_year.max().values)
            minimums.append(this_year.min().values)
            time_min  = this_year["time"].values[this_year.argmin("time").values]
            time_max  = this_year["time"].values[this_year.argmax("time").values]
            min_doys.append(pd.to_datetime(time_min).dayofyear)
            max_doys.append(pd.to_datetime(time_max).dayofyear)
            onset_doy   = self.detect_onset_via_derivative(this_year)
            retreat_doy = self.compute_retreat_via_inflection(this_year)
            duration_days.append(retreat_doy - onset_doy if retreat_doy > onset_doy else (365 - onset_doy + retreat_doy))
            if onset_doy:
                onset_doys.append(onset_doy)          
            grow = this_year.sel(time=((this_year.time.dt.dayofyear >= growth_range[0]) & (this_year.time.dt.dayofyear <= growth_range[1])))
            if grow.time.size >= 2:
                days  = (grow.time - grow.time[0]) / np.timedelta64(1, 'D')
                slope = np.polyfit(days, grow.values.flatten(), 1)[0]
                growth_rates.append(slope * 1e6)
            early = this_year.sel(time=((this_year.time.dt.dayofyear >= retreat_early_range[0]) & (this_year.time.dt.dayofyear <= retreat_early_range[1])))
            if early.time.size >= 2:
                days  = (early.time - early.time[0]) / np.timedelta64(1, 'D')
                slope = np.polyfit(days, early.values.flatten(), 1)[0]
                retreat_early_rates.append(slope * 1e6)
            fall      = this_year.sel(time=this_year.time.dt.dayofyear >= retreat_late_range1[0])
            next_year = da.sel(time=da.time.dt.year == (year + 1))
            spring    = next_year.sel(time=next_year.time.dt.dayofyear <= retreat_late_range2[1])
            late      = xr.concat([fall, spring], dim="time")
            if late.time.size >= 2:
                days  = (late.time - late.time[0]) / np.timedelta64(1, 'D')
                slope = np.polyfit(days, late.values.flatten(), 1)[0]
                retreat_late_rates.append(slope * 1e6)
            if early.time.size >= 1 and late.time.size >= 1:
                decay = xr.concat([early, late], dim="time")
                if decay.time.size >= 2:
                    days  = (decay.time - decay.time[0]) / np.timedelta64(1, 'D')
                    slope = np.polyfit(days, decay.values.flatten(), 1)[0]
                    rate  = slope * 1e6
                    if rate > 0:
                        rate *= -1
                    retreat_rates.append(rate)
        return {"Maximum Mean"       : np.mean(maximums),
                "Maximum Std"        : np.std(maximums),
                "Minimum Mean"       : np.mean(minimums),
                "Minimum Std"        : np.std(minimums),
                "Growth Mean"        : np.mean(growth_rates),
                "Growth Std"         : np.std(growth_rates),
                "retreat Mean"       : np.mean(retreat_rates),
                "retreat Std"        : np.std(retreat_rates),
                "retreat Early Mean" : np.mean(retreat_early_rates),
                "retreat Early Std"  : np.std(retreat_early_rates),
                "retreat Late Mean"  : np.mean(retreat_late_rates),
                "retreat Late Std"   : np.std(retreat_late_rates),
                "Duration-days Mean" : np.mean(duration_days),
                "Duration-days std"  : np.std(duration_days),
                "DOY Min Mean"       : np.mean(min_doys),
                "DOY Min Std"        : np.std(min_doys),
                "DOY Max Mean"       : np.mean(max_doys),
                "DOY Max Std"        : np.std(max_doys),
                "DOY Onset Mean"     : np.mean(onset_doys) if onset_doys else None,
                "DOY Onset Std"      : np.std(onset_doys)  if onset_doys else None}

    # STABILITY METRICS
    def compute_fip_spatial_stats(self, FIP):
        valid = FIP.where(~np.isnan(FIP))
        return {'FIP_spatial_mean' : float(valid.mean()),
                'FIP_spatial_std'  : float(valid.std())}

    def compute_cellwise_stability(self, FI_mask):
        total = FI_mask.sizes['time']
        return FI_mask.sum(dim='time') / total

    def compute_stability_index(self, persistent_FIA, total_FIA):
        return persistent_FIA / total_FIA if total_FIA > 0 else np.nan

    def compute_fast_ice_distance_extent(self, FI_mask, grid_dx_km=9.8):
        """
        Compute mean and maximum distance (in km) of fast ice from the coastline.

        Parameters
        ----------
        FI_mask : xarray.DataArray
        Boolean mask of fast ice presence over time (time, nj, ni)
        grid_dx_km : float
        Approximate horizontal grid spacing in kilometers (default 9.8 km; Antarctic coast)

        Returns
        -------
        mean_dist : float
        Mean distance from coast for fast ice presence
        max_dist : float
        Maximum distance from coast for fast ice presence
        """
        if self.use_gi:
            P_kmt = self.P_KMT_mod
        else:
            P_kmt = self.P_KMT_org
        kmt             = xr.open_dataset(P_kmt)['kmt']
        kmt             = kmt.isel(nj=self.hemisphere_dict['nj_slice']).values
        land_mask       = (kmt == 0)
        sea_mask        = ~land_mask
        coast_mask      = sea_mask & binary_dilation(land_mask)
        coast_distances = distance_transform_edt(~coast_mask) * grid_dx_km
        TLAT = FI_mask.TLAT
        TLON = FI_mask.TLON
        if TLAT.ndim == 3:
            TLAT = TLAT.isel(time=0)
        if TLON.ndim == 3:
            TLON = TLON.isel(time=0)
        coast_dist_da   = xr.DataArray(data   = coast_distances,
                                       dims   = ('nj', 'ni'),
                                       coords = {'nj'   : FI_mask.nj,
                                                 'ni'   : FI_mask.ni,
                                                 'TLAT' : (('nj','ni'), TLAT.values),
                                                 'TLON' : (('nj','ni'), TLON.values)})
        fi_mask_time_mean = FI_mask.mean(dim="time") > 0.5
        fast_ice_dists    = coast_dist_da.where(fi_mask_time_mean)
        mean_dist         = float(fast_ice_dists.mean().values)
        max_dist          = float(fast_ice_dists.max().values)
        return {'FI_mean_dist_ext' : mean_dist,
                'FI_max_dist_ext'  : max_dist}

    # INTER-COMPARISONS
    def compute_skill_statistics(self, model, obs, dropna=True):
        """
        Compute standard statistics between two aligned time series.

        Parameters
        ----------
        model : xr.DataArray or np.ndarray
        obs   : xr.DataArray or np.ndarray
        dropna : bool
            Drop NaNs before computing metrics (default: True)

        Returns
        -------
        dict
            Dictionary of Bias, RMSE, MAE, Corr, SD_model, SD_obs

        Notes:
        -------------------------------------
        Metric	                Definition	                        Good For
        Bias	                Mean(model - obs)	                Systematic offset
        RMSE	                √(mean((model - obs)²))	            Total error magnitude
        RMSD	                Same as RMSE 	                    Shape error (variance around mean)
        MAE	                    Mean absolute error	Error           magnitude without squaring
        Correlation	Pearson     correlation coefficient	            Phase & shape agreement
        Standard deviation	    Variability	                        Amplitude bias
        Skill Score	            Normalized error relative to cliimatology Relative improvement over a reference baseline        
        """
        from sklearn.metrics import mean_squared_error, mean_absolute_error
        from scipy.stats import pearsonr        
        t_xsect = np.intersect1d(model.time.values, obs.time.values)
        model   = model.sel(time=t_xsect)
        obs     = obs.sel(time=t_xsect)        
        if isinstance(model, xr.DataArray):
            model = model.values
        if isinstance(obs, xr.DataArray):
            obs = obs.values
        if dropna:
            valid = (~np.isnan(model)) & (~np.isnan(obs))
            model = model[valid]
            obs = obs[valid]
        bias     = np.mean(model - obs)
        rmse     = np.sqrt(mean_squared_error(obs, model))
        mae      = mean_absolute_error(obs, model)
        corr, _  = pearsonr(model, obs)
        sd_model = np.std(model)
        sd_obs   = np.std(obs)
        return {"Bias"    : bias,
                "RMSE"    : rmse,
                "MAE"     : mae,
                "Corr"    : corr,
                "SD_Model": sd_model,
                "SD_Obs"  : sd_obs,}

    def compute_taylor_stats(self, da1, da2, chunk_size=30):
        all_model, all_obs = [], []
        for i in range(0, da1.sizes['time'], chunk_size):
            a = da1.isel(time=slice(i, i+chunk_size)).compute().values.flatten()
            b = da2.isel(time=slice(i, i+chunk_size)).compute().values.flatten()
            valid = np.isfinite(a) & np.isfinite(b)
            if valid.sum() > 0:
                all_model.append(a[valid])
                all_obs.append(b[valid])
        model, obs         = np.concatenate(all_model), np.concatenate(all_obs)
        corr               = np.corrcoef(model, obs)[0, 1]
        std_model, std_obs = np.std(model), np.std(obs)
        rmsd               = np.sqrt(np.mean((model - obs - (np.mean(model) - np.mean(obs)))**2))
        return {"corr": corr, "std_ratio": std_model / std_obs, "rmsd_ratio": rmsd / std_obs}