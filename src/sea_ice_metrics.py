import xarray as xr
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.ndimage import binary_dilation, distance_transform_edt    
    
class SeaIceMetrics:
    def __init__(self, **kwargs):
        return

    def extract_pack_ice(self, DS_FI, DS_HEM):
        """
        Extract the pack ice portion of the original dataset using a precomputed fast ice mask.

        This method uses the inverse of the fast ice mask (`PI_mask`) to isolate pack ice regions 
        from the full sea ice dataset. It assumes that `DS_FI` (a fast ice output) contains a valid 
        boolean `PI_mask` derived from previous masking logic.

        INPUTS:
            DS_FI  : xarray.Dataset; dataset containing the boolean `PI_mask` variable indicating pack ice regions.
            DS_HEM : xarray.Dataset; original unmasked sea ice dataset (e.g., the Southern Hemisphere slice of `iceh`).

        OUPUTS:
            xarray.Dataset; a masked version of `DS_HEM` where only pack ice regions are retained.

        """
        assert "PI_mask" in DS_FI, "PI_mask not found in dataset"
        self.logger.info("masking for pack ice on dataset")
        DS_PI = DS_HEM.where(DS_FI["PI_mask"])
        return DS_PI

    def compute_rolling_mean_on_dataset(self, ds, mean_period=None):
        """

        Apply a centered temporal rolling mean to a dataset.

        This method smooths temporal noise by computing a centered moving average across the time dimension,
        typically for use in rolling fast ice classification.

        INPUTS:
           ds          : xarray.Dataset; input dataset with a `time` dimension.
           mean_period : int, optional; rolling window size in days. Defaults to `self.mean_period`.

        OUTPUTS:
           xarray.Dataset; dataset with all variables averaged over the specified rolling window.

        """
        mean_period = mean_period if mean_period is not None else self.mean_period
        return ds.rolling(time=mean_period, center=True, min_periods=1).mean()

    def boolean_fast_ice(self, FI_mask, dim="time", window=7, min_count=6):
        """

        Generate a boolean (or binary-days) presence mask for fast ice using a rolling window threshold.

        This method applies a binary persistence filter to a fast ice mask. A cell is considered to
        contain fast ice if it meets or exceeds `min_count` valid (True) values within a rolling window
        of length `window`.

        INPUTS:
           FI_mask   : xarray.DataArray; Boolean time series mask of fast ice presence (True/False).
           dim       : str, optional; the name of the time dimension. Defaults to "time".
           window    : int, optional; size of the rolling window in days. Defaults to `self.bool_window` or 7.
           min_count : int, optional; minimum number of True values in the window to classify as fast ice.
                       Defaults to `self.bool_min_days` or 6.

        OUTPUTS:
           xarray.DataArray; Boolean mask where True indicates persistent fast ice presence.

        NOTES:
        + Used for defining stable fast ice coverage (e.g., in climatologies or seasonality studies).
        + This operation is memory-intensive and is automatically persisted with Dask.

        """
        window    = window    if window    is not None else self.bool_window
        min_count = min_count if min_count is not None else self.bool_min_days
        self.logger.info(f"🔁 Rolling boolean presence: window = {window}, min_count = {min_count}")
        FI_roll_mask = FI_mask.rolling({dim: window}, center=True).construct(f"{dim}_window").sum(dim=f"{dim}_window")
        FI_bool      = (FI_roll_mask >= min_count)
        return FI_bool

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
                roll        = i_type.endswith("_roll")
                DS, CICE_SO = self.load_processed_cice(ispd_thresh = ispd_thresh,
                                                       ice_type    = ice_type,
                                                       zarr_CICE   = True,
                                                       rolling     = roll,
                                                       slice_hem   = True)
                if i_type==f"{ice_type}_bool":
                    bool_mask          = self.boolean_fast_ice(DS['FI_mask'], dim="time", window=7, min_count=6)
                    DS_bool            = CICE_SO.where(bool_mask)
                    DS_bool["FI_mask"] = DS["FI_mask"]
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

    def compute_ice_area(self, SIC, GC_area, ice_area_scale=None):
        if self.use_gi:
            GI_total_area = self.compute_grounded_iceberg_area()
        else:
            GI_total_area = 0
        self.logger.info(f"{GI_total_area:0.2f} m^2 total circumpolar grounded iceberg area for {self.sim_name}")
        ice_area_scale = ice_area_scale if ice_area_scale is not None else self.FIC_scale
        self.logger.info(f"🧮 Spatially-integrating the product of sea ice concentrations and grid cell areas")
        IA = ((SIC * GC_area).sum(dim=("nj", "ni"))).persist()
        IA = IA + GI_total_area
        IA = IA/ice_area_scale
        return IA

    def compute_ice_volume(self, SIC, HI, GC_area, ice_volume_scale=1e12):
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
        GI_total_area = self.compute_grounded_iceberg_area() if self.use_gi else 0
        self.logger.info(f"{GI_total_area:0.2f} m² total grounded iceberg area for {self.sim_name}")
        self.logger.info("🧮 Computing ice volume: sum(SIC × HI × area)")
        IV = (SIC * HI * GC_area).sum(dim=("nj", "ni")).persist()
        IV += GI_total_area
        IV /= ice_volume_scale
        return IV


    def compute_variable_aggregate(self, da, time_coord_name='time'):
        return da.sum(dim=time_coord_name) / da[time_coord_name].sizes.get(time_coord_name, 1)

    # GROWTH METRICS
    def compute_fia_onset_doy(self, FIA, threshold_km2=50):
        t = FIA["time"].dt.dayofyear
        onset = FIA.where(FIA > threshold_km2, drop=True)
        return int(t.sel(time=onset.time[0]).item()) if onset.size > 0 else None

    def compute_fia_growth_rate(self, FIA, start_doy=71, end_doy=273):
        """
        Compute linear growth rate of FIA between start_doy and end_doy.
        Default start_doy and end_doy obtained from https://doi.org/10.5194/tc-15-5061-2021
        """
        t_doy = FIA["time"].dt.dayofyear
        mask  = (t_doy >= start_doy) & (t_doy <= end_doy)
        FIA_sub = FIA.where(mask, drop=True)
        if FIA_sub.time.size < 2:
            return np.nan
        days = (FIA_sub.time - FIA_sub.time[0]) / np.timedelta64(1, 'D')
        slope = np.polyfit(days, FIA_sub.values, 1)[0]
        return slope

    def compute_fia_max_growth(self, FIA):
        dFIA = FIA.diff("time")
        return dFIA.max().item()

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

    # SEASONALITY METRICS
    def compute_doy_max(self, FIA):
        idx = FIA.argmax("time")
        return int(FIA["time"].dt.dayofyear[idx].item())

    def compute_fast_ice_duration(self, FIA, threshold_km2=None):
        """
        Compute duration of the fast ice season.

        If `threshold_km2` is given: duration is the number of days with FIA > threshold.
        If `threshold_km2` is None: duration is defined as (DOY_max - DOY_min) per year.

        Parameters
        ----------
        FIA : xarray.DataArray
            Daily fast ice area with time coordinate.
        threshold_km2 : float or None
            Fixed threshold for defining start/end. If None, uses annual min/max to define season.

        Returns
        -------
        dict[int, int] or int
            Duration per year in days, or single int if only one year.
        """
        if "time" not in FIA.dims:
            raise ValueError("FIA must have 'time' coordinate.")
        if threshold_km2 is not None:
            mask = FIA > threshold_km2
            if mask.sum() == 0:
                return 0
            times = FIA["time"].where(mask, drop=True)
            return int((times[-1] - times[0]) / np.timedelta64(1, "D")) + 1
        durations = {}
        for yr, da in FIA.groupby("time.year"):
            if da.time.size < 2:
                durations[yr.item()] = 0
                continue
            t_doy = da["time"].dt.dayofyear
            i_min = da.argmin("time").item()
            i_max = da.argmax("time").item()
            t0, t1 = int(t_doy[i_min].item()), int(t_doy[i_max].item())
            durations[yr.item()] = abs(t1 - t0)
        if len(durations) == 1:
            return list(durations.values())[0]
        return durations

    def compute_fia_rmse(self, model_fia, obs_fia):
        if model_fia.sizes["time"] != obs_fia.size:
            raise ValueError("Time dimension mismatch between model and observations.")
        error = model_fia.values - obs_fia.values
        return float(np.sqrt(np.mean(error**2)))