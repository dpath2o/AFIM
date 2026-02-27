from __future__ import annotations
import xarray as xr
import numpy  as np
import pandas as pd
from pathlib  import Path
__all__ = ["SeaIceFast"]
class SeaIceFast:

    def __init__():
        return

    def coarsen_and_align_simulated_FI_to_observed_FI(self, sim_ds, obs_ds,
                                                      doy_vals=None,
                                                      method="mean",
                                                      obs_time_coord="time"):
        """
        Aggregate daily simulated fields onto observational fast-ice time windows.

        This method coarsens daily model output into the same temporal windows used by
        the AF2020 observational product (typically 24 bins defined by day-of-year start
        values). For each year and each DOY window, it aggregates model data over the
        corresponding date range and writes one output timestep per window.

        Parameters
        ----------
        sim_ds : xarray.Dataset
            Daily model output with a "time" coordinate.
        obs_ds : xarray.Dataset
            Observational dataset providing the target window definitions via its time
            coordinate (default name is `obs_time_coord`).
        doy_vals : list[int], optional
            Day-of-year start values defining observation windows. If None, uses
            `self.AF_FI_dict["DOY_vals"]`.
        method : {"mean","median"}, default "mean"
            Aggregation method to apply within each window.
        obs_time_coord : str, default "time"
            Name of the time coordinate in `obs_ds` used as the output time dimension.

        Returns
        -------
        xarray.Dataset
            Aggregated dataset with one timestep per observation window (per year),
            concatenated along `obs_time_coord`.

        Raises
        ------
        ValueError
            If no windows intersect the simulation data or `method` is unsupported.

        Notes
        -----
        - Window end dates are computed from successive DOY boundaries. The final window
        extends to the end of the year (accounting for leap years).
        - The output time stamps correspond to the start date of each window.
        """
        if doy_vals is None:
            doy_vals = self.AF_FI_dict["DOY_vals"]
        sim_time = sim_ds["time"].values
        obs_times = pd.to_datetime(obs_ds[obs_time_coord].values)
        years = np.unique(obs_times.year)
        grouped = []
        for year in years:
            for i, doy_start in enumerate(doy_vals):
                dt_start = pd.Timestamp(f"{year}-01-01") + pd.to_timedelta(doy_start - 1, unit="D")
                dt_end = pd.Timestamp(f"{year}-01-01") + pd.to_timedelta(
                    (doy_vals[i+1]-1 if i+1 < len(doy_vals) else 365 + int(pd.Timestamp(f"{year}-12-31").is_leap_year)), unit="D")
                period_mask = (sim_ds.time >= dt_start) & (sim_ds.time < dt_end)
                ds_window = sim_ds.sel(time=period_mask)
                if ds_window.time.size == 0:
                    self.logger.warning(f"🕳️ No model data in window {dt_start.date()} to {dt_end.date()}")
                    continue
                if method == "mean":
                    ds_agg = ds_window.mean(dim="time")
                elif method == "median":
                    ds_agg = ds_window.median(dim="time")
                else:
                    raise ValueError(f"Unsupported method: {method}")
                ds_agg = ds_agg.expand_dims({obs_time_coord: [dt_start]})
                grouped.append(ds_agg)
        if not grouped:
            raise ValueError("❌ No observation periods matched simulation data")
        ds_aligned = xr.concat(grouped, dim=obs_time_coord)
        self.logger.info(f"✅ Aligned model output to {len(ds_aligned[obs_time_coord])} obs windows")
        return ds_aligned

    def compute_fipdiff_stats_weighted(self, FIP: xr.Dataset,
                                       PIXEL                 : float = 5000.0,
                                       regs_dict             : dict | None = None,
                                       gi_mask_da            : xr.DataArray | None = None,
                                       clip01                : bool = True,
                                       threecat_percent_only : bool = True):
        regs_dict = regs_dict if regs_dict is not None else self.Ant_8sectors  # adjust name if needed
        for k in ("mod", "obs", "lon", "lat"):
            if k not in FIP:
                raise ValueError("FIP must contain 'mod', 'obs', 'lon', 'lat'.")
        mod = FIP["mod"].astype("float64")
        obs = FIP["obs"].astype("float64")
        if clip01:
            mod = mod.clip(0.0, 1.0)
            obs = obs.clip(0.0, 1.0)
        overlap      = xr.apply_ufunc(np.minimum, mod, obs)
        model_excess = xr.apply_ufunc(np.maximum, mod - obs, 0.0)
        obs_excess   = xr.apply_ufunc(np.maximum, obs - mod, 0.0)
        lon_m = self.normalise_longitudes(FIP["lon"], to="-180-180")
        lat   = FIP["lat"]
        cell_area_km2 = (float(PIXEL) * float(PIXEL)) / 1e6
        if gi_mask_da is not None:
            valid = (~gi_mask_da.astype(bool).broadcast_like(mod))
        else:
            valid = xr.ones_like(mod, dtype=bool)
        # loop over regional dictionary terms to create rows 
        rows = []
        for rname, reg_dict in regs_dict.items():
            geo_reg = reg_dict["plot_region"]  # must be (lon_min, lon_max, lat_min, lat_max) in -180..180
            R = self._region_mask(lon_m, lat, geo_reg) & valid
            A_km2 = float((overlap      * R).sum()) * cell_area_km2
            M_km2 = float((model_excess * R).sum()) * cell_area_km2
            O_km2 = float((obs_excess   * R).sum()) * cell_area_km2
            tot   = A_km2 + M_km2 + O_km2
            pct_A = 100.0 * A_km2 / tot if tot > 0 else np.nan
            pct_M = 100.0 * M_km2 / tot if tot > 0 else np.nan
            pct_O = 100.0 * O_km2 / tot if tot > 0 else np.nan
            MOD_km2 = float((mod * R).sum()) * cell_area_km2
            OBS_km2 = float((obs * R).sum()) * cell_area_km2
            pct_A_of_MOD = 100.0 * A_km2 / MOD_km2 if MOD_km2 > 0 else np.nan
            pct_M_of_MOD = 100.0 * M_km2 / MOD_km2 if MOD_km2 > 0 else np.nan
            pct_A_of_OBS = 100.0 * A_km2 / OBS_km2 if OBS_km2 > 0 else np.nan
            pct_O_of_OBS = 100.0 * O_km2 / OBS_km2 if OBS_km2 > 0 else np.nan
            rows.append(dict(region                    = rname,
                             pct_agreement             = pct_A,
                             pct_model_dominant        = pct_M,
                             pct_observation_dominant  = pct_O,
                             model_FIA_km2             = MOD_km2,
                             obs_FIA_km2               = OBS_km2,
                             agreement_overlap_km2     = A_km2,
                             model_excess_km2          = M_km2,
                             obs_excess_km2            = O_km2,
                             agreement_pct_of_model    = pct_A_of_MOD,
                             model_excess_pct_of_model = pct_M_of_MOD,
                             agreement_pct_of_obs      = pct_A_of_OBS,
                             obs_excess_pct_of_obs     = pct_O_of_OBS))
        # create dataframe from dictionary
        df      = pd.DataFrame(rows).set_index("region")
        A_sum   = df["agreement_overlap_km2"].sum()
        M_sum   = df["model_excess_km2"].sum()
        O_sum   = df["obs_excess_km2"].sum()
        tot_sum = A_sum + M_sum + O_sum
        # now sum together for pan-Ant. ... but what about NH??
        ant_row = dict(pct_agreement             = 100.0 * A_sum / tot_sum if tot_sum > 0 else np.nan,
                       pct_model_dominant        = 100.0 * M_sum / tot_sum if tot_sum > 0 else np.nan,
                       pct_observation_dominant  = 100.0 * O_sum / tot_sum if tot_sum > 0 else np.nan,
                       model_FIA_km2             = df["model_FIA_km2"].sum(),
                       obs_FIA_km2               = df["obs_FIA_km2"].sum(),
                       agreement_overlap_km2     = A_sum,
                       model_excess_km2          = M_sum,
                       obs_excess_km2            = O_sum,
                       agreement_pct_of_model    = np.nan,
                       model_excess_pct_of_model = np.nan,
                       agreement_pct_of_obs      = np.nan,
                       obs_excess_pct_of_obs     = np.nan)
        df_all = pd.concat([df, pd.DataFrame(ant_row, index=["ANT"])])
        if threecat_percent_only:
            return df_all[["pct_agreement", "pct_model_dominant", "pct_observation_dominant"]]
        else:
            return df_all[["model_FIA_km2", "agreement_overlap_km2", "model_excess_km2",
                           "agreement_pct_of_model", "model_excess_pct_of_model",
                           "obs_FIA_km2", "agreement_pct_of_obs", "obs_excess_km2", "obs_excess_pct_of_obs",
                           "pct_agreement", "pct_model_dominant", "pct_observation_dominant"]]

    def define_FIP_diff_names(self, 
                     AF2020            = False,
                     variable_name     = None,
                     lon_coord_name    = None,
                     lat_coord_name    = None,
                     time_coord_name   = None,
                     spatial_dim_names = None,
                     thresh_val        = None):
        """
        Standardise coordinate/variable naming for fast-ice persistence diagnostics.

        This is a lightweight helper that returns a consistent naming dictionary used by
        fast-ice mask and coordinate extraction routines for either:
          - AF2020 observational datasets (AF2020=True), or
          - model/CICE datasets (AF2020=False).

        Parameters
        ----------
        AF2020 : bool, default False
            If True, derive names and threshold from `self.AF_FI_dict`; otherwise use
            CICE/model conventions.
        variable_name, lon_coord_name, lat_coord_name, time_coord_name : str, optional
            Override defaults for variable and coordinate names.
        spatial_dim_names : tuple[str, str], optional
            Override spatial dim names. Defaults to `self.CICE_dict["spatial_dims"]`.
        thresh_val : float, optional
            Override the threshold value used to define "fast ice" presence.

        Returns
        -------
        dict
            Keys:
            - spatial_dims, time_dim, variable, latitude, longitude, thresh_val
        """
        name_space        = {}
        spatial_dim_names = spatial_dim_names if spatial_dim_names is not None else self.CICE_dict['spatial_dims']
        if AF2020:
            name_var   = variable_name   if variable_name   is not None else self.AF_FI_dict["variable_name"]
            name_lat   = lat_coord_name  if lat_coord_name  is not None else self.AF_FI_dict["lat_coord_name"]
            name_lon   = lon_coord_name  if lon_coord_name  is not None else self.AF_FI_dict["lon_coord_name"]
            name_dt    = time_coord_name if time_coord_name is not None else self.AF_FI_dict["time_coord_name"]
            thresh_val = thresh_val      if thresh_val      is not None else self.AF_FI_dict["threshold_value"]
        else:    
            name_var   = variable_name   if variable_name   is not None else "aice"
            name_lat   = lat_coord_name  if lat_coord_name  is not None else self.CICE_dict["tcoord_names"][1]
            name_lon   = lon_coord_name  if lon_coord_name  is not None else self.CICE_dict["tcoord_names"][0]
            name_dt    = time_coord_name if time_coord_name is not None else self.CICE_dict["time_dim"]
            thresh_val = thresh_val      if thresh_val      is not None else self.icon_thresh
        return {'spatial_dims': spatial_dim_names,
                'time_dim'    : name_dt,
                'variable'    : name_var,
                'latitude'    : name_lat,
                'longitude'   : name_lon,
                'thresh_val'  : thresh_val}

    def define_fast_ice_coordinates(self, ds,
                                    AF2020            = False,
                                    variable_name     = None,
                                    lon_coord_name    = None,
                                    lat_coord_name    = None,
                                    time_coord_name   = None,
                                    spatial_dim_names = None,
                                    thresh_val        = None):
        """
        Extract lat/lon/time arrays and standardised naming for fast-ice mask construction.

        Parameters
        ----------
        ds : xr.Dataset
            Dataset containing the required variable and coordinate fields.
        AF2020 : bool, default False
            Switch between AF2020 and CICE/model conventions.
        (other parameters)
            Passed through to `define_FIP_diff_names(...)`.

        Returns
        -------
        dict
            Keys:
            - datetimes : 1D datetime array
            - latitudes : 2D latitude array
            - longitudes : 2D longitude array
            - names : the naming dictionary returned by `define_FIP_diff_names`
        """
        names = self.define_FIP_diff_names(AF2020            = AF2020,
                                           variable_name     = variable_name,
                                           lon_coord_name    = lon_coord_name,
                                           lat_coord_name    = lat_coord_name,
                                           time_coord_name   = time_coord_name,
                                           spatial_dim_names = spatial_dim_names,
                                           thresh_val        = thresh_val)
        if "time" in ds[names['latitude']].dims:
            lat = ds[names['latitude']].isel(time=0).values
        else:
            lat = ds[names['latitude']].values
        if "time" in ds[names['longitude']].dims:
            lon = ds[names['longitude']].isel(time=0).values
        else:
            lon = ds[names['longitude']].values
        dt = ds[names['time_dim']].values
        return {'datetimes'  : dt,
                'latitudes'  : lat,
                'longitudes' : lon,
                'names'      : names}

    def define_fast_ice_mask_da(self, ds,
                                AF2020            = False,
                                variable_name     = None,
                                lon_coord_name    = None,
                                lat_coord_name    = None,
                                time_coord_name   = None,
                                spatial_dim_names = None,
                                thresh_val        = None):
        """
        Construct a binary fast-ice mask DataArray from a dataset.

        A grid cell is classified as fast ice when:
            ds[variable] >= thresh_val
        where names/threshold are determined by `define_FIP_diff_names()`.

        Parameters
        ----------
        ds : xr.Dataset
            Input dataset.
        AF2020 : bool, default False
            If True, use AF2020 conventions; else CICE/model.
        (other parameters)
            Passed through to `define_fast_ice_coordinates(...)`.

        Returns
        -------
        xr.DataArray
            Binary mask with dims given by `self.CICE_dict['three_dims']` and coordinates
            attached using the configured lat/lon/time coordinate names.

        Notes
        -----
        - This routine assumes the mask should always be expressed on the CICE naming
          conventions (coords from `self.CICE_dict`), even when sourced from AF2020.
        """
        coords  = self.define_fast_ice_coordinates(ds,
                                                    AF2020            = AF2020,
                                                    variable_name     = variable_name,
                                                    lon_coord_name    = lon_coord_name,
                                                    lat_coord_name    = lat_coord_name,
                                                    time_coord_name   = time_coord_name,
                                                    spatial_dim_names = spatial_dim_names,
                                                    thresh_val        = thresh_val)
        FI_mask = xr.where(ds[coords['names']['variable']] >= coords['names']['thresh_val'], 1.0, 0.0)
        return xr.DataArray(FI_mask,
                            dims   = self.CICE_dict['three_dims'],
                            coords = {self.CICE_dict['lat_coord_name'] : (self.CICE_dict['spatial_dims'], coords['latitudes']),
                                      self.CICE_dict['lon_coord_name'] : (self.CICE_dict['spatial_dims'], coords['longitudes']),
                                      self.CICE_dict['time_dim']       : (self.CICE_dict['time_dim']    , coords['datetimes'])})

    def define_fast_ice_persistence_da(self, fip, reG=False):
        """
        Wrap a persistence-like array into a consistently-labelled 2D DataArray.

        Parameters
        ----------
        fip : xr.DataArray
            Persistence field on either the native grid or a regridded regular grid.
        reG : bool, default False
            If True, attach coordinates from `define_regular_G(...)` (0.15° global SH);
            otherwise use fip's existing lat/lon coordinates.

        Returns
        -------
        xr.DataArray
            2D float32 DataArray with spatial dims `self.CICE_dict['spatial_dims']` and
            lat/lon coordinates attached.
        """
        if reG:
            G_dst = self.define_regular_G(0.15, region=[0,360,-90,0], spatial_dim_names=self.CICE_dict['spatial_dims'])
            lats  = G_dst[self.CICE_dict['lat_coord_name']].values
            lons  = G_dst[self.CICE_dict['lon_coord_name']].values
        else:
            lats = fip[self.CICE_dict['lat_coord_name']].values
            lons = fip[self.CICE_dict['lon_coord_name']].values
        return xr.DataArray(fip.astype('float32').values,
                            dims   = self.CICE_dict['spatial_dims'],
                            coords = {self.CICE_dict['lat_coord_name'] : (self.CICE_dict['spatial_dims'], lats),
                                    self.CICE_dict['lon_coord_name'] : (self.CICE_dict['spatial_dims'], lons)})

    def compute_fip_weight(self,
                           FIP   : xr.Dataset, 
                           mode  : str = "max",
                           t     : float = 0.1,
                           gamma : float = 1.0) -> xr.DataArray:
        """
        Compute an opacity/weight field from fast-ice persistence coverage.

        This is used for visualisation or for weighting a misfit field, producing
        values in [0, 1] where:
          - 0 => fully transparent / no weight (below threshold)
          - 1 => fully opaque / full weight (high persistence coverage)

        Parameters
        ----------
        FIP : xr.Dataset
            Must include `FIP["obs"]` (observed persistence fraction). May include
            `FIP["mod"]` (model persistence fraction).
        mode : {"max","mean","prod","obs"}, default "max"
            Coverage operator used before thresholding:
              - "max"  : max(obs, mod)
              - "mean" : 0.5*(obs + mod)
              - "prod" : sqrt(obs*mod) (geometric mean)
              - any other value falls back to obs only
        t : float, default 0.1
            Soft threshold. Values below t map to 0; values linearly ramp to 1 by 1.0.
        gamma : float, default 1.0
            Contrast shaping exponent (w ** gamma). gamma>1 increases contrast.

        Returns
        -------
        xr.DataArray
            Weight field in [0, 1], named "diff_weight", with descriptive attributes.

        Raises
        ------
        KeyError
            If FIP does not contain an "obs" variable.

        Notes
        -----
        - This routine is agnostic to the grid; it simply transforms persistence fractions.
        """
        if "obs" not in FIP:
            raise KeyError("FIP['obs'] is required")
        obs = FIP["obs"]
        mod = FIP.get("mod", None)
        if mode == "max" and mod is not None:
            cov = xr.apply_ufunc(np.maximum, obs, mod)
        elif mode == "mean" and mod is not None:
            cov = 0.5 * (obs + mod)
        elif mode == "prod" and mod is not None:
            cov = xr.apply_ufunc(np.sqrt, obs * mod)  # geometric mean
        else:
            cov = obs  # fall back to obs only
        # Soft ramp: 0 below t, then linear to 1 at 1.0
        w = (cov - t) / (1 - t)
        w = w.clip(0, 1)
        # Optional contrast shaping
        if gamma != 1.0:
            w = w ** gamma
        w.name = "diff_weight"
        w.attrs.update(dict(long_name="opacity weight", comment=f"weight ~ {mode}(obs,mod), threshold={t}, gamma={gamma}"))
        return w

