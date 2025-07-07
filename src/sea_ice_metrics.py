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
                                ispd_thresh     = None,
                                D_out           = None,
                                overwrite_zarr  = False,
                                overwrite_png   = False,
                                smooth_FIA_days = 15):
        """
        Compute and plot a suite of fast ice metrics from processed simulation output.

        This wrapper method automates the loading, computation, and visualization of key
        fast ice metrics for a given simulation. It handles all three processing modes
        (raw, rolling mean, and boolean persistence) and generates Zarr-based metrics files
        as well as regional and hemispheric PNG plots.

        INPUTS:
           sim_name        : str, optional; simulation name (defaults to `self.sim_name`).
           ice_type        : str; the base fast ice type to process (e.g., "FI_B", "FI_BT").
           ispd_thresh     : float, optional; threshold for ice speed masking. Required for Zarr path construction.
           D_out           : Path or str, optional; output directory for Zarr and CSV metric files.
                             Defaults to simulation Zarr path under `"metrics/"`.
           overwrite_zarr  : bool, optional; if True, recompute metrics and overwrite existing Zarr files.
           overwrite_png   : bool, optional; if True, regenerate PNG figures even if they already exist.
           smooth_FIA_days : int, optional; smoothing window (in days) for plotting the FIA time series.

        OUTPUTS:
           None; results are written to disk and plotted using `SeaIcePlotter`.

        NOTES:
           + Calls `load_processed_cice()` for data loading, `compute_sea_ice_metrics()` for analysis,
             and `SeaIcePlotter` for figure generation.
           + This is the primary public method for computing FIA, FIP, and climatology comparisons.
           + Requires access to AF2020 climatology (`P_AF2020_cli_csv`) for observational benchmarking.

        """
        sim_name                 = sim_name if sim_name is not None else self.sim_name
        ispd_thresh              = ispd_thresh or self.ispd_thresh
        ispd_thresh_str          = f"{ispd_thresh:.1e}".replace("e-0", "e-")
        ice_type                 = ice_type or self.ice_type
        dt_range_str             = f"{self.dt0_str[:4]}-{self.dtN_str[:4]}"
        FIA_dict                 = {}
        af2020_df                = pd.read_csv(self.AF_FI_dict['P_AF2020_cli_csv'])
        FIA_dict["AF2020db_cli"] = self.interpolate_obs_fia(af2020_df)
        ice_types                = [ice_type, f"{ice_type}_roll", f"{ice_type}_bool"]
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
                DS, CICE_SO = self.load_processed_cice(ispd_thresh = ispd_thresh,
                                                       ice_type    = ice_type,
                                                       zarr_CICE   = True,
                                                       rolling     = roll,
                                                       slice_hem   = True)
                if i_type==f"{ice_type}_bool":
                    bool_mask          = self.boolean_fast_ice(DS['FI_mask'], dim="time", window=7, min_count=6)
                    DS_bool            = CICE_SO.where(bool_mask)
                    #DS_bool["FI_mask"] = DS["FI_mask"]
                    DS_bool["FI_mask"] = bool_mask
                    DS                 = DS_bool
                METS = self.compute_sea_ice_metrics(DS, sim_name, i_type, self.ispd_thresh_str, P_METS, P_sum, FIA_dict["AF2020db_cli"])
            FIA_dict[i_type] = METS['FIA']
        P_png = Path(self.D_graph, sim_name, f"FIA_FIP_{sim_name}_{ispd_thresh_str}_{dt_range_str}.png")
        self.plot_FIA_FIP_faceted(FIA_dict, METS['FIP'], P_png=P_png, plot_GI=True if self.use_gi else False)

    def compute_sea_ice_metrics(self, DS, sim_name, i_type, ispd_thresh_str, P_METS, P_sum, obs_clim):
        """

        Compute and persist diagnostic metrics describing fast ice coverage and seasonality.

        This method evaluates a suite of spatial and temporal diagnostics from a given fast ice dataset:
        + Area-integrated fast ice area/extent (FIA)
        + Fast ice persistence (FIP)
        + Onset timing, growth dynamics, duration, and spatial statistics
        + Optional RMSE against observational climatology (AF2020)

        Results are stored as a Zarr file and also exported as a CSV summary.

        INPUTS:
           DS              : xarray.Dataset; fast ice dataset including `aice` and `FI_mask`.
           sim_name        : str; simulation name used for file naming and CSV metadata.
           i_type          : str; fast ice group identifier (e.g., "FI_B", "FI_BT_bool").
           ispd_thresh_str : str; stringified ice speed threshold (e.g., "1.0e-3") for output naming.
           P_METS          : Path; path to output `.zarr` file storing full metrics dataset.
           P_sum           : Path; path to output `.csv` file storing scalar summary metrics.
           obs_clim        : xarray.DataArray or None; observational climatology of fast ice area (FIA)
                             for model comparison.

        OUTPUTS:
           xarray.Dataset; dataset containing 3D spatial metrics (FIP), time series metrics (FIA),
           and scalar diagnostics.

        NOTES:
        + FIA and FIP are computed directly from `aice` and `FI_mask`, then used to derive other metrics.
        + Summary includes:
            + Onset DOY, max growth, growth rate, season duration
            + FIP mean/std
            + Spatial distance extent
            + Optional RMSE to observational climatology
        + Automatically writes outputs to disk (Zarr + CSV).

        """
        METS = {}
        # 3D + 1D metrics
        FIA = self.compute_ice_area(DS['aice'], DS['tarea']).compute()
        FIP = self.compute_variable_aggregate(DS['aice']).compute()
        METS["FIA"] = FIA
        METS["FIP"] = FIP
        # Scalar / 1D metrics
        summary = {}
        summary["onset_doy"]    = self.compute_fia_onset_doy(FIA)
        summary["growth_rate"]  = self.compute_fia_growth_rate(FIA)
        summary["max_growth"]   = self.compute_fia_max_growth(FIA)
        summary["doy_max"]      = self.compute_doy_max(FIA)
        summary["duration"]     = self.compute_fast_ice_duration(FIA)
        fip_mean, fip_std       = self.compute_fip_spatial_stats(FIP)
        summary["FIP_mean"]     = fip_mean
        summary["FIP_std"]      = fip_std
        mean_dist, max_dist     = self.compute_fast_ice_distance_extent(DS['FI_mask'])
        summary["mean_FI_dist"] = mean_dist
        summary["max_FI_dist"]  = max_dist
        if obs_clim is not None:
            model_doy = FIA["time"].dt.dayofyear.values
            obs_vals = np.interp(model_doy, obs_clim.coords["doy"].values, obs_clim.values)
            summary["rmse_to_obs"] = self.compute_fia_rmse(FIA, xr.DataArray(obs_vals, coords=[("time", FIA["time"].data)]))
        # 🔁 Convert any dicts (e.g. duration) to xarray.DataArray
        for key, val in summary.items():
            if isinstance(val, dict):
                summary[key] = xr.DataArray(pd.Series(val), dims="year")
        # Combine all into a dataset
        DS_METS = xr.Dataset(summary)
        DS_METS["FIA"] = FIA
        DS_METS["FIP"] = FIP
        DS_METS.to_zarr(P_METS, mode="w", consolidated=True)
        self.logger.info(f"📊 Metrics written to {P_METS}")
        df = pd.DataFrame([DS_METS])
        df["sim_name"]    = sim_name
        df["ice_type"]    = i_type
        df["ispd_thresh"] = ispd_thresh_str
        df.to_csv(P_sum, index=False)
        self.logger.info(f"📊 Metrics summary written to {P_sum}")
        return DS_METS

    def compute_ice_area(self, SIC, GC_area, ice_area_scale=None, spatial_dim_names=None, add_grounded_iceberg_area=None, grounded_iceberg_area=None):
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
        IA = ((SIC * GC_area).sum(dim=spatial_dim_names)).persist()
        IA = IA + GI_total_area
        IA = IA/ice_area_scale
        return IA

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
        IV = (SIC * HI * GC_area).sum(dim=spatial_dim_names).persist()
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
        return float(valid.mean()), float(valid.std())

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
        return mean_dist, max_dist

    # INTER-COMPARISONS
    def compute_statistics(self, model, obs, dropna=True):
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