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
                                dt0_str         = None,
                                dtN_str         = None,
                                ice_type        = None,
                                ispd_thresh     = None,
                                overwrite_zarr  = None,
                                overwrite_png   = None):
        """
        Compute, load, and plot summary fast ice metrics (FIA and FIP) for a given simulation and configuration.

        This method wraps the full processing and plotting pipeline for fast ice area (FIA) and fast ice persistence (FIP)
        by optionally loading pre-computed metrics from disk (Zarr), or re-computing them from processed CICE model output.
        It also loads the observational AF2020 climatology for comparison and generates a faceted figure comparing
        multiple fast ice masking types (raw, rolling, and boolean). Metrics are written to disk if not already present
        or if overwrite is enabled.

        Parameters
        ----------
        sim_name : str, optional
            Name of the model simulation to process. Defaults to `self.sim_name`.
        dt0_str : str, optional
            Start date (inclusive) of the processing window in 'YYYY-MM-DD' format. Defaults to `self.dt0_str`.
        dtN_str : str, optional
            End date (inclusive) of the processing window in 'YYYY-MM-DD' format. Defaults to `self.dtN_str`.
        ice_type : str, optional
            Type of fast ice mask to use (e.g., 'FI_BT'). Defaults to `self.ice_type`.
        ispd_thresh : float, optional
            Ice speed threshold (in m/s) used to define fast ice. Defaults to `self.ispd_thresh`.
        overwrite_zarr : bool, optional
            If True, re-compute and overwrite existing metric Zarr files. Defaults to `self.overwrite_zarr_group`.
        overwrite_png : bool, optional
            If True, overwrite the output PNG plot file. Defaults to `self.overwrite_png`.

        Notes
        -----
        - The method processes and compares three variants of the fast ice mask:
          1. Raw (e.g., 'FI_BT') – direct thresholding of daily speed
          2. Rolling mean (e.g., 'FI_BT_roll') – speed is smoothed before thresholding
          3. Boolean (e.g., 'FI_BT_bool') – based on a persistence filter applied after masking

        - Each mask is used to compute and/or load:
          * FIA (Fast Ice Area) — integrated area over time
          * FIP (Fast Ice Persistence) — fraction of time with fast ice presence

        - Observational climatology is loaded from AF2020 via `load_AF2020_FIA_summary()`.

        - All computed metrics are saved as Zarr and CSV files under `self.D_metrics`.

        - The final plot is saved under `self.D_graph/sim_name/`.

        Output
        ------
        - Saves summary metrics as:
          * Zarr: {ice_type}_mets.zarr
          * CSV:  {ice_type}_summary.csv
        - Generates faceted figure:
          * PNG: FIA_FIP_{sim_name}_{ispd_thresh}_{year-range}.png

        See Also
        --------
        - self.load_AF2020_FIA_summary
        - self.load_processed_cice
        - self.boolean_fast_ice
        - self.compute_sea_ice_metrics
        - self.plot_FIA_FIP_faceted
        """
        sim_name           = sim_name       if sim_name       is not None else self.sim_name
        ice_type           = ice_type       if ice_type       is not None else self.ice_type
        dt0_str            = dt0_str        if dt0_str        is not None else self.dt0_str
        dtN_str            = dtN_str        if dtN_str        is not None else self.dtN_str
        ispd_thresh        = ispd_thresh    if ispd_thresh    is not None else self.ispd_thresh
        overwrite_zarr     = overwrite_zarr if overwrite_zarr is not None else self.overwrite_zarr_group
        overwrite_png      = overwrite_png  if overwrite_png  is not None else self.ow_fig
        ispd_thresh_str    = f"{ispd_thresh:.1e}".replace("e-0", "e-")
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
        self.plot_FIA_FIP_faceted(FIA_dict, METS['FIP'],
                                  P_png         = P_png,
                                  overwrite_fig = overwrite_png,
                                  plot_GI       = True if self.use_gi else False)

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
        sim_config = self.sim_config
        # --- Merge All Metrics into Dataset ---
        summary = {**{f"FIA_{k}": v for k, v in FIA_seasonal.items()},
                   **{f"FIV_{k}": v for k, v in FIV_seasonal.items()},
                   **FIP_stats,
                   **FI_dist,
                   **FIA_skill,
                   **sim_config}
        # Add summary statistics as scalar variables
        for k, v in summary.items():
            if k in sim_config:
                DS_METS.attrs[k] = v  # Store sim_config as dataset-level attributes
            else:
                DS_METS[k] = xr.DataArray(v, dims=())  # Summary metrics remain as DataArrays
        # --- Save to Zarr (Optional) ---
        if P_METS is not None:
            DS_METS.to_zarr(P_METS, mode="w", consolidated=True)
            self.logger.info(f"📊 Metrics written to {P_METS}")
        return DS_METS

    def compute_ice_area(self, SIC, GC_area,
                         ice_area_scale            = None,
                         spatial_dim_names         = None,
                         sic_threshold             = None,
                         add_grounded_iceberg_area = None,
                         grounded_iceberg_area     = None):
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

        sic_threshold : float, optional
            Minimum SIC value to be included in the area calculation. Points with SIC below this threshold are excluded.
            Defaults to `self.icon_thresh`.
            
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
        sic_threshold             = sic_threshold             if sic_threshold             is not None else self.icon_thresh
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
        self.logger.info(f"spatially-integrating the product of sea ice concentrations and grid cell areas")
        mask    = SIC > sic_threshold
        SIC     = SIC.where(mask)
        GC_area = GC_area.where(mask)
        IA      = ((SIC * GC_area).sum(dim=spatial_dim_names))
        IA      = IA + GI_total_area
        IA      = IA/ice_area_scale
        return IA

    def compute_ice_volume(self, SIC, HI, GC_area,
                           ice_volume_scale  = 1e12, 
                           spatial_dim_names = None,
                           sic_threshold     = None):
        """
        Compute total sea ice volume by integrating sea ice concentration, thickness, and grid cell area.

        This method calculates the total sea ice volume as the sum of the product of sea ice concentration (SIC),
        thickness (HI), and grid cell area (GC_area) across the model domain. A SIC threshold can be applied to exclude
        grid cells with low concentrations. Optionally includes grounded iceberg volume and applies a scaling factor 
        for unit conversion.

        Parameters
        ----------
        SIC : xarray.DataArray
            Sea ice concentration (unitless, typically between 0 and 1).
            
        HI : xarray.DataArray
            Sea ice thickness in meters.
            
        GC_area : xarray.DataArray
            Grid cell area in square meters (m²).

        ice_volume_scale : float, optional
            Scale factor for the output volume. Default is `1e12`, converting m³ to 1000 km³.

        spatial_dim_names : list of str, optional
            Names of spatial dimensions to sum over (e.g., ['nj', 'ni']). Defaults to `self.CICE_dict['spatial_dims']`.

        sic_threshold : float, optional
            Minimum SIC value to be included in the volume calculation. Grid cells with SIC ≤ `sic_threshold` are masked.
            Defaults to `self.icon_thresh`.

        Returns
        -------
        float
            Total sea ice volume, optionally including grounded ice and scaled (e.g., in 1000 km³ if default scale is used).

        Notes
        -----
        - SIC, HI, and GC_area must share the same grid and be broadcastable.
        - Grid cells below `sic_threshold` are excluded from the volume calculation.
        - If `self.use_gi` is True, a grounded iceberg area is converted to volume and included.
        """
        spatial_dim_names = spatial_dim_names if spatial_dim_names is not None else self.CICE_dict['spatial_dims']
        sic_threshold     = sic_threshold     if sic_threshold     is not None else self.icon_thresh
        GI_total_area = self.compute_grounded_iceberg_area() if self.use_gi else 0
        self.logger.info(f"{GI_total_area:0.2f} m² total grounded iceberg area for {self.sim_name}")
        mask    = SIC > sic_threshold
        HI      = HI.where(mask)
        SIC     = SIC.where(mask)
        GC_area = GC_area.where(mask)
        self.logger.info("computing ice volume: sum(SIC × HI × area)")
        IV = (SIC * HI * GC_area).sum(dim=spatial_dim_names)
        IV += GI_total_area
        IV /= ice_volume_scale
        return IV

    def compute_ice_thickness(self, HI, SIC, GC_area, spatial_dim_names=None, sic_threshold=None):
        """
        Compute average sea ice thickness weighted by grid cell area and sea ice concentration.

        This method calculates the domain-averaged sea ice thickness as the ratio of the total sea ice volume to 
        the total sea ice area. A SIC threshold is used to exclude low-concentration grid cells from both the numerator 
        and denominator.

        Parameters
        ----------
        HI : xarray.DataArray
            Sea ice thickness in meters.

        SIC : xarray.DataArray
            Sea ice concentration (unitless, typically between 0 and 1).

        GC_area : xarray.DataArray
            Grid cell area in square meters (m²).

        spatial_dim_names : list of str, optional
            Names of spatial dimensions over which to compute the sums. Defaults to `self.CICE_dict['spatial_dims']`.

        sic_threshold : float, optional
            Minimum SIC value to include in the average. Grid cells with SIC ≤ `sic_threshold` are excluded.
            Defaults to `self.icon_thresh`.

        Returns
        -------
        float or xarray.DataArray
            Domain-averaged sea ice thickness, computed as the total volume divided by total area.

        Notes
        -----
        - Computed as: ∑(HI × area) / ∑(SIC × area), with masking by `sic_threshold`.
        - Excludes cells with SIC below threshold from both numerator and denominator.
        - Units are in meters.
        """
        spatial_dim_names = spatial_dim_names if spatial_dim_names is not None else self.CICE_dict['spatial_dims']
        sic_threshold     = sic_threshold     if sic_threshold     is not None else self.icon_thresh
        self.logger.info(f"masking for sea ice concentration above {sic_threshold}")
        mask    = SIC > sic_threshold
        HI      = HI.where(mask)
        SIC     = SIC.where(mask)
        GC_area = GC_area.where(mask)
        self.logger.info(f"computing area and volume")
        sia     = (SIC * GC_area).sum(dim=spatial_dim_names)
        siv     = (HI * GC_area).sum(dim=spatial_dim_names)
        return siv / sia

    def compute_variable_aggregate(self, da, time_coord_name='time'):
        """
        Compute the time-mean aggregate of a variable by summing over time and normalizing by time dimension length.

        This method computes the average of a DataArray along a given time dimension using the total sum divided by 
        the number of time steps. Useful for handling Dask arrays or preserving lazy evaluation.

        Parameters
        ----------
        da : xarray.DataArray
            The input data array to be aggregated over time.
        time_coord_name : str, optional
            Name of the time coordinate/dimension. Defaults to `'time'`.

        Returns
        -------
        xarray.DataArray
            Time-averaged data array.

        Notes
        -----
        - Falls back to division by 1 if the time dimension is missing.
        - This is equivalent to `da.mean(dim='time')` but avoids triggering eager computation.
        """
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
        """
        Estimate the day of sea ice retreat based on the latest inflection point in a time series.

        This method uses the second derivative of a smoothed time series to identify the last significant
        inflection point indicating retreat (change from decelerating to accelerating melt) during the melt season.
        The signal is smoothed using a Savitzky–Golay filter.

        Parameters
        ----------
        da : xarray.DataArray
            A 1D time series of sea ice area (or similar metric) with a `time` coordinate.
        window : int, optional
            Window length (in days) for the Savitzky–Golay filter. Must be an odd integer.
            Default is 15.
        polyorder : int, optional
            Polynomial order used in the Savitzky–Golay filter. Default is 2.
        min_retreat_doy : int, optional
            Minimum day-of-year (DOY) to begin searching for inflection points. Default is 240 (~end of August).

        Returns
        -------
        int or None
            Day-of-year (DOY) of the last negative inflection point (i.e., where curvature becomes negative),
            indicating sea ice retreat. Returns `None` if no such inflection point is found after `min_retreat_doy`.

        Notes
        -----
        - Inflection points are defined where the second derivative (`d²y/dt²`) becomes negative.
        - This method assumes the time series is sufficiently smooth and uniformly spaced in time.
        - The last inflection point after `min_retreat_doy` is returned as the estimated retreat day.
        - No output is returned if no qualifying inflection points are found.
        """
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
        Compute seasonal growth, retreat, and timing statistics for a sea ice metric (e.g., SIA, FIA, SIV, FIV).

        This method loops over all years in the provided time series and calculates:
        - Minimum and maximum values and their day-of-year (DOY)
        - Growth rate during the ice advance season
        - Retreat rate, split into early (spring) and late (autumn–next year) seasons
        - Overall retreat rate (from both early and late phases)
        - Duration of the seasonal ice cycle (onset to retreat)
        - DOY of growth onset (based on derivative)
        - DOY of retreat (based on inflection analysis)

        Parameters
        ----------
        da : xarray.DataArray
            Time series of a sea ice variable (must have a 'time' coordinate). Should be continuous and span multiple years.
        growth_range : tuple of int, optional
            Day-of-year (DOY) range used to compute the seasonal **growth** rate. Default is (71, 273), i.e. mid-March to late-September.
        retreat_early_range : tuple of int, optional
            DOY range for computing **early melt/retreat** rate (e.g., spring melt). Default is (273, 330).
        retreat_late_range1 : tuple of int, optional
            DOY range for the **late melt** phase in the current year (e.g., late December). Default is (331, 365).
        retreat_late_range2 : tuple of int, optional
            DOY range for the **late melt** phase in the following year (e.g., January). Default is (1, 71).

        Returns
        -------
        dict
            Dictionary of summary statistics, including:
            - "Maximum Mean", "Maximum Std"
            - "Minimum Mean", "Minimum Std"
            - "Growth Mean", "Growth Std"
            - "retreat Mean", "retreat Std"
            - "retreat Early Mean", "retreat Early Std"
            - "retreat Late Mean", "retreat Late Std"
            - "Duration-days Mean", "Duration-days std"
            - "DOY Min Mean", "DOY Min Std"
            - "DOY Max Mean", "DOY Max Std"
            - "DOY Onset Mean", "DOY Onset Std" (if detectable)

        Notes
        -----
        - Only non-leap years are included in the analysis to maintain consistent DOY ranges.
        - The first year is excluded to ensure late-season comparisons can access next year's data.
        - Growth and retreat rates are computed as the slope (Δvalue / Δdays) × 1e6 for clearer units (e.g., 10⁶ km²/day).
        - Retreat DOY is detected using the inflection point of a smoothed curve (via `compute_retreat_via_inflection`).
        - Onset DOY is computed from the maximum positive slope in early-season growth (via `detect_onset_via_derivative`).
        - Retreat duration handles wrap-around years (e.g., onset in September, retreat in March).

        See Also
        --------
        - self.compute_retreat_via_inflection : Detects retreat date via curvature analysis
        - self.detect_onset_via_derivative    : Detects growth onset via slope analysis
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
        """
        Compute spatial mean and standard deviation of Fast Ice Persistence (FIP).

        This method calculates the mean and standard deviation of FIP values over a spatial domain,
        ignoring any NaN values (i.e., outside the coastline or invalid grid cells).

        Parameters
        ----------
        FIP : xarray.DataArray
            A 2D (or broadcastable) array of fast ice persistence values ranging from 0 to 1,
            where each cell represents the fraction of time with fast ice presence over a given period.

        Returns
        -------
        dict
            Dictionary containing:
            - 'FIP_spatial_mean' : float
                Mean FIP across valid spatial cells.
            - 'FIP_spatial_std' : float
                Standard deviation of FIP across valid spatial cells.

        Notes
        -----
        - NaNs are excluded from all statistics using `xarray.where()`.
        - Intended to summarize spatial characteristics of fast ice persistence across the domain.
        """
        valid = FIP.where(~np.isnan(FIP))
        return {'FIP_spatial_mean' : float(valid.mean()),
                'FIP_spatial_std'  : float(valid.std())}

    def compute_cellwise_stability(self, FI_mask):
        """
        Compute cellwise fast ice stability as the fraction of time fast ice is present.

        This method calculates, for each grid cell, the proportion of time steps where fast ice is detected,
        using a boolean mask over time.

        Parameters
        ----------
        FI_mask : xarray.DataArray
            A boolean DataArray of shape (time, nj, ni), where True indicates presence of fast ice
            and False/NaN indicates absence.

        Returns
        -------
        xarray.DataArray
            A 2D field (nj, ni) of fractional fast ice stability, where each value represents the
            fraction of time that fast ice was present.

        Notes
        -----
        - Assumes the 'time' dimension exists and is the first dimension.
        - This metric can be interpreted as the temporal persistence of fast ice at each location.
        """
        total = FI_mask.sizes['time']
        return FI_mask.sum(dim='time') / total

    def compute_stability_index(self, persistent_FIA, total_FIA):
        """
        Compute the fast ice stability index as the ratio of persistent to total fast ice area.

        This index quantifies the stability of fast ice by comparing the area of long-lived (persistent)
        fast ice to the total fast ice area over the same period.

        Parameters
        ----------
        persistent_FIA : float
            Area (in consistent units, e.g., km²) of fast ice that is persistent across the full season
            or meets a defined persistence threshold.
        total_FIA : float
            Total fast ice area (in same units) accumulated or detected over the season.

        Returns
        -------
        float
            Stability index, defined as persistent_FIA / total_FIA.

        Notes
        -----
        - Returns NaN or inf if `total_FIA` is zero.
        - Values close to 1 indicate stable fast ice with little seasonal turnover,
          whereas values near 0 indicate transient or unstable ice conditions.
        """
        return persistent_FIA / total_FIA if total_FIA > 0 else np.nan

    def compute_fast_ice_distance_extent(self, FI_mask, grid_dx_km=9.8):
        """
        Compute the mean and maximum distance of fast ice from the coastline.

        This method estimates how far fast ice extends from the coast by computing the
        Euclidean distance of each grid cell to the nearest coastal cell (defined as the
        boundary between land and ocean). It then computes statistics on those distances,
        weighted by the mean presence of fast ice.

        Parameters
        ----------
        FI_mask : xarray.DataArray
            Boolean mask of fast ice presence (shape: time × nj × ni), where True indicates fast ice.
        grid_dx_km : float, optional
            Approximate horizontal grid spacing in kilometers (default is 9.8 km for ~0.25° Antarctic grid).

        Returns
        -------
        dict
            Dictionary containing:
            - 'FI_mean_dist_ext' : float
                Mean distance (km) of fast ice presence from the coastline.
            - 'FI_max_dist_ext' : float
                Maximum distance (km) of fast ice from the coastline.

        Notes
        -----
        - The coastline is computed from the land mask in the model's `kmt` field.
        - A morphological dilation is used to identify coastal boundary cells.
        - Distance is computed using a 2D Euclidean distance transform.
        - Fast ice extent is summarized from the time-mean of the mask (`FI_mask.mean(dim="time") > 0.5`).
        - Results are reported in kilometers.
        - Requires `self.P_KMT_mod` and `self.P_KMT_org` to point to NetCDF files with a `kmt` variable.
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
        Compute statistical skill metrics between model and observational time series.

        This method compares two aligned 1D time series and returns standard evaluation metrics
        including bias, RMSE, MAE, correlation, and standard deviations.

        Parameters
        ----------
        model : xarray.DataArray or np.ndarray
            Time series of modelled values.
        obs : xarray.DataArray or np.ndarray
            Time series of observed values (must be aligned with model).
        dropna : bool, optional
            Whether to drop NaN values prior to computing statistics. Default is True.

        Returns
        -------
        dict
            Dictionary of metrics including:
            - "Bias"     : Mean(model - obs)
            - "RMSE"     : Root Mean Squared Error
            - "MAE"      : Mean Absolute Error
            - "Corr"     : Pearson correlation coefficient
            - "SD_Model" : Standard deviation of model
            - "SD_Obs"   : Standard deviation of observations

        Notes
        -----
        - Time alignment is ensured using `np.intersect1d()` on time coordinates if xarray is used.
        - If `dropna=True`, only valid (non-NaN) pairs are used in the statistics.
        - RMSE captures total error magnitude, while MAE gives average absolute error.
        - Pearson correlation (`Corr`) quantifies phase agreement.
        - The method assumes that input arrays represent the same variable over the same time span.

        See Also
        --------
        - sklearn.metrics.mean_squared_error
        - sklearn.metrics.mean_absolute_error
        - scipy.stats.pearsonr
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
