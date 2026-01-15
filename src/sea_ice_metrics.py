import gc, re, dask
import xarray      as xr
import numpy       as np
import pandas      as pd
from pathlib       import Path
class SeaIceMetrics:

    def __init__(self, **kwargs):
        return

    @staticmethod
    def _clean_zarr_chunks(ds):
        for var in ds.variables:
            if 'chunks' in ds[var].encoding:
                del ds[var].encoding['chunks']
        return ds.chunk(None)

    def fast_ice_metrics_data_dict(self, FI_mask, FI_data, A):
        """
        Utility function to create the dictionary for each FI type (dy, rl, bn).
        """
        return {"FI_mask"  : FI_mask,
                'aice'     : FI_data['aice'],
                'hi'       : FI_data['hi'],
                #'strength' : FI_data['strength'],
                'dvidtt'   : FI_data['dvidtt'],
                'dvidtd'   : FI_data['dvidtd'],
                'daidtt'   : FI_data['daidtt'],
                'daidtd'   : FI_data['daidtd'],
                'tarea'    : A}


    def pack_ice_metrics_data_dict(self, PI_mask, PI_data, A):
        # Utility function to create the dictionary for each PI type (dy, rl, bn).
        return {"PI_mask"  : PI_mask,
                'aice'     : PI_data['aice'],
                'hi'       : PI_data['hi'],
                #'strength' : PI_data['strength'],
                'dvidtt'   : PI_data['dvidtt'],
                'dvidtd'   : PI_data['dvidtd'],
                'daidtt'   : PI_data['daidtt'],
                'daidtd'   : PI_data['daidtd'],
                'tarea'    : A}

    def sea_ice_metrics_data_dict(self, SI_mask, SI_data, A):
        # Utility function to create the dictionary for SI (daily-only; no binary-day association).
        return {"SI_mask"  : SI_mask,
                'aice'     : SI_data['aice'],
                'hi'       : SI_data['hi'],
                #'strength' : SI_data['strength'],
                'dvidtt'   : SI_data['dvidtt'],
                'dvidtd'   : SI_data['dvidtd'],
                'daidtt'   : SI_data['daidtt'],
                'daidtd'   : SI_data['daidtd'],
                'tarea'    : A}

    def compute_sea_ice_metrics(self, da_dict,
                                ice_type       = None,
                                dt0_str        = None,
                                dtN_str        = None,
                                P_mets_zarr    = None,
                                ice_area_scale = None):
        """
        Compute sea ice metrics from a processed sea ice simulation and compare against observations.

        Parameters
        ----------
        da_dict : dict
            Dictionary containing xarray DataArrays or Datasets for the fields:
            'aice', 'hi', 'tarea', and a corresponding mask ('FI_mask', 'PI_mask', or 'SI_mask').
        ice_type : str, optional
            Type of ice to process (e.g., 'FI', 'FI_BT_bin', etc.). Will be cleaned internally.
        dt0_str : str, optional
            Start date of the analysis window in 'YYYY-MM-DD'. Defaults to self.dt0_str.
        dtN_str : str, optional
            End date of the analysis window in 'YYYY-MM-DD'. Defaults to self.dtN_str.
        P_mets_zarr : str or Path, optional
            Path to save Zarr archive of computed metrics. If None, uses default path.
        ice_area_scale : float, optional
            Scaling factor to apply to the area field. Defaults to self.FIC_scale.

        Returns
        -------
        xarray.Dataset
            Dataset containing time series, spatial, and summary metrics.
        """
        ice_type       = ice_type  or self.ice_type
        dt0_str        = dt0_str   or self.dt0_str
        dtN_str        = dtN_str   or self.dtN_str
        P_mets_zarr    = Path(P_mets_zarr) if P_mets_zarr else Path(self.D_metrics, f"{ice_type}_mets.zarr")
        ice_area_scale = ice_area_scale if ice_area_scale is not None else self.FIC_scale
        spatial_dim_names = self.CICE_dict['spatial_dims']
        self.logger.info(f"¡¡¡ COMPUTING ICE METRICS for {ice_type} !!!")
        ice_type_clean = re.sub(r'(_roll|_bin)?(_BT|_B|_Ta|_Tx)?$', '', ice_type)
        if ice_type_clean not in ("FI", "PI", "SI"):
            raise ValueError(f"Unexpected cleaned ice_type: {ice_type_clean}")
        self.logger.info(f"results will be written to {P_mets_zarr}")
        mask_name = f"{ice_type_clean}_mask"
        if mask_name in da_dict:
            I_mask = da_dict[mask_name]
        else:
            candidates = [k for k in ("FI_mask", "PI_mask", "SI_mask") if k in da_dict]
            if len(candidates) == 1:
                mask_name = candidates[0]
                I_mask = da_dict[mask_name]
                self.logger.warning(f"Mask key {ice_type_clean}_mask not found; inferring {mask_name} from da_dict.")
            else:
                raise KeyError(f"Could not resolve ice mask. Expected {ice_type_clean}_mask, or exactly one of FI_mask/PI_mask/SI_mask. Found: {candidates}")
        I_C       = da_dict['aice']
        I_T       = da_dict['hi']
        #I_S       = da_dict['strength']
        I_TVT     = da_dict['dvidtt']
        I_MVT     = da_dict['dvidtd']
        I_TAT     = da_dict['daidtt']
        I_MAT     = da_dict['daidtd']
        A         = da_dict['tarea']
        # --- Time Series Metrics ---
        IA   = self.compute_hemisphere_ice_area(I_C, A, ice_area_scale=ice_area_scale)
        IV   = self.compute_hemisphere_ice_volume(I_C, I_T, A)
        IT   = self.compute_hemisphere_ice_thickness(I_C, I_T, A)
        #IS   = self.compute_hemisphere_ice_strength(I_C, I_T, I_S)
        ITVR = self.compute_hemisphere_ice_volume_rate(I_C, I_TVT)
        IMVR = self.compute_hemisphere_ice_volume_rate(I_C, I_MVT)
        ITAR = self.compute_hemisphere_ice_area_rate(I_TAT, IA, A)
        IMAR = self.compute_hemisphere_ice_area_rate(I_MAT, IA, A)
        IP   = self.compute_hemisphere_variable_aggregate(I_C)
        self.logger.info("computing **ICE THICKNESS TEMPORAL-MEAN**")
        IHI = I_T.mean(dim=self.CICE_dict["time_dim"])
        #self.logger.info("computing **ICE STRENGTH TEMPORAL-SUM**; units mPa")
        #IST = (I_S/I_T).sum(dim=self.CICE_dict["time_dim"]) / 1e6
        self.logger.info("computing **ICE VOLUME TENDENCY (SPATIAL RATE)**; units m/yr")
        ITVR_YR = (I_TVT*1e2).mean(dim=self.CICE_dict["time_dim"]) / 3.65
        IMVR_YR = (I_MVT*1e2).mean(dim=self.CICE_dict["time_dim"]) / 3.65
        self.logger.info("computing **ICE AREA TENDENCY (SPATIAL RATE)**; units m/yr")
        ITAR_YR = (I_TAT*A).mean(dim=self.CICE_dict["time_dim"]) / 31_536_000
        IMAR_YR = (I_MAT*A).mean(dim=self.CICE_dict["time_dim"]) / 31_536_000
        self.logger.info("loading data into output dictionary...")
        METS = {f"{ice_type_clean}A"      : IA.load(),      #1D
                f"{ice_type_clean}V"      : IV.load(),      #1D
                f"{ice_type_clean}T"      : IT.load(),      #1D
                #f"{ice_type_clean}S"      : IS.load(),      #1D
                f"{ice_type_clean}TVR"    : ITVR.load(),    #1D
                f"{ice_type_clean}MVR"    : IMVR.load(),    #1D
                f"{ice_type_clean}TAR"    : ITAR.load(),    #1D
                f"{ice_type_clean}MAR"    : IMAR.load(),    #1D
                f"{ice_type_clean}P"      : IP.load(),      #2D
                f"{ice_type_clean}HI"     : IHI.load(),     #2D
                #f"{ice_type_clean}ST"     : IST.load(),     #2D
                f"{ice_type_clean}TVR_YR" : ITVR_YR.load(), #2D
                f"{ice_type_clean}MVR_YR" : IMVR_YR.load(), #2D
                f"{ice_type_clean}TAR_YR" : ITAR_YR.load(), #2D
                f"{ice_type_clean}MAR_YR" : IMAR_YR.load()} #2D
        # --- Skill Statistics ---
        try:
            if ice_type_clean == "FI":
                self.logger.info(f"loading FI obs: {self.AF_FI_dict['P_AF2020_FIA']}")
                IA_obs = xr.open_dataset(self.AF_FI_dict['P_AF2020_FIA'], engine="netcdf4")["AF2020"]
            elif ice_type_clean in ("PI", "SI"):
                NSIDC = self.compute_NSIDC_metrics()
                IA_obs = NSIDC['SIA']
            else:
                IA_obs = None
            IA_obs   = IA_obs.load()
            IA_skill = self.compute_skill_statistics(IA, IA_obs) if IA_obs is not None else {}
        except Exception as e:
            self.logger.warning(f"compute_skill_statistics failed: {e}")
            IA_skill = {}
        # --- Seasonal Statistics ---
        try:
            IA_seasonal = self.compute_seasonal_statistics(IA, stat_name=f"{ice_type_clean}A")
        except Exception as e:
            self.logger.warning(f"compute_seasonal_statistics failed: {e}")
            IA_seasonal = {}
        # --- Persistence Statistics ---
        if ice_type_clean == "FI":
            try:
                IP_stab = self.persistence_stability_index(I_mask, A)
            except Exception as e:
                self.logger.warning(f"persistence_stability_index failed: {e}")
                IP_stab = {}
            try:
                IP_dist = self.persistence_ice_distance_mean_max(IP)
            except Exception as e:
                self.logger.warning(f"persistence_ice_distance_mean_max failed: {e}")
                IP_dist = {}
        else:
            IP_stab, IP_dist = {}, {}
        # --- Build Output Dataset ---
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
        # --- Merge Metadata ---
        summary = {**{f"{ice_type_clean}A_{k}": v for k, v in IA_seasonal.items()},
                   **IP_stab, **IP_dist, **IA_skill, **self.sim_config}
        for k, v in summary.items():
            if k in self.sim_config:
                DS_METS.attrs[k] = v
            else:
                DS_METS[k] = xr.DataArray(v, dims=())
        # --- Save to Zarr ---
        if P_mets_zarr:
            DS_METS = self._clean_zarr_chunks(DS_METS)
            DS_METS.to_zarr(P_mets_zarr, mode="w", consolidated=True, zarr_format=2)
            self.logger.info(f"Metrics written to {P_mets_zarr}")
        return DS_METS

    def _subset_and_pad_time(self, ds: xr.Dataset, dt0_str: str, dtN_str: str, 
                             time_dim: str = "time") -> xr.Dataset:
        if time_dim not in ds.dims and time_dim not in ds.coords:
            # Nothing to do
            return ds
        if dt0_str is None or dtN_str is None:
            return ds
        t = ds[time_dim]
        if t.size == 0:
            return ds
        # Determine whether we're working with numpy datetime64 or cftime objects
        t0 = t.values[0]
        is_datetime64 = np.issubdtype(t.values.dtype, np.datetime64)
        # Build start/end in matching type
        if is_datetime64:
            start = np.datetime64(dt0_str)
            end   = np.datetime64(dtN_str)
        else:
            # cftime (dtype typically object)
            try:
                import cftime  # noqa: F401
            except Exception as e:
                raise RuntimeError("Dataset time coordinate appears to be cftime/object, but cftime is not available.") from e
            # Try to preserve calendar if present; default to 'standard'
            cal = t.encoding.get("calendar", None) or t.attrs.get("calendar", None) or "standard"
            start = xr.cftime_range(start=dt0_str, periods=1, calendar=cal)[0]
            end   = xr.cftime_range(start=dtN_str, periods=1, calendar=cal)[0]
        # Sanity
        if start > end:
            raise ValueError(f"Requested start > end: {dt0_str} > {dtN_str}")
        # Infer step from dataset time (fallback: 1 day)
        if t.size >= 2:
            step = t.values[1] - t.values[0]
        else:
            step = np.timedelta64(1, "D") if is_datetime64 else __import__("datetime").timedelta(days=1)
        # Build the requested time axis (inclusive)
        if is_datetime64:
            # ensure inclusive end (np.arange is end-exclusive)
            desired = np.arange(start, end + step, step)
        else:
            desired = []
            cur = start
            # robust loop for cftime + timedelta
            while cur <= end:
                desired.append(cur)
                cur = cur + step
            desired = np.array(desired, dtype=object)
        # Reindex onto the desired axis -> pads with NaNs where outside original coverage
        # (keeps original values where they overlap)
        ds_out = ds.reindex({time_dim: desired})
        return ds_out
    
    def load_computed_metrics(self,
                              fast_ice_class_method : str  = "binary-days",  # "raw", "rolling-mean", "binary-days"
                              BorC2T_type           : str  = None,
                              ice_type              : str  = None,
                              ispd_thresh           : str  = None,
                              zarr_directory        : str  = None,
                              clip_to_self          : bool = True,
                              time_dim              : str  = "time"):
        BorC2T_type = BorC2T_type    or self.BorC2T_type
        ice_type    = ice_type       or self.ice_type
        ispd_thresh = ispd_thresh    or self.ispd_thresh
        D_zarr      = zarr_directory or self.D_zarr
        self.define_classification_dir(ice_type = ice_type, D_zarr = D_zarr, ispd_thresh = ispd_thresh)
        if (ice_type == 'FI') or (ice_type =='PI'): 
            self.define_fast_ice_class_name(BorC2T_type = BorC2T_type, fast_ice_class_method = fast_ice_class_method)
            P_mets = self.D_class / f"{self.FI_class}_{self.metrics_name}.zarr"
        else:
            P_mets = self.D_class / f"{ice_type}_{self.metrics_name}.zarr"
        ds = xr.open_dataset(P_mets)
        if clip_to_self:
            ds = self._subset_and_pad_time(ds, self.dt0_str, self.dtN_str, time_dim=time_dim)
        return ds

    def compute_hemisphere_ice_area_rate(self, DAT, IA, A,
                                           spatial_dim_names  : list  = None):
        """
        """
        spatial_dim_names = spatial_dim_names if spatial_dim_names is not None else self.CICE_dict['spatial_dims']
        self.logger.info("Computing **DYNAMIC AREA TENDENCY** for hemisphere")
        self.logger.debug(f"  • Spatial dimension names           : {spatial_dim_names}")
        self.logger.debug("\n[Sea Ice Pressure Computation Steps]\n"
                          f"  1. multiply dynamic area tendency (1/s) by grid cell area (m^2) and sum\n"
                          f"  2. divide by ice area and volumetric scaling factor 1e9")
        return (DAT*A).sum(dim=spatial_dim_names) / IA / 1e9 # m/s

    def compute_hemisphere_ice_volume_rate(self, SIC, DVT,
                                           spatial_dim_names  : list  = None,
                                           sic_threshold      : float = None):
        """
        """
        sic_threshold     = sic_threshold     if sic_threshold     is not None else self.icon_thresh
        spatial_dim_names = spatial_dim_names if spatial_dim_names is not None else self.CICE_dict['spatial_dims']
        self.logger.info("Computing **DYNAMIC VOLUME TENDENCY** for hemisphere")
        self.logger.debug(f"  • Using SIC threshold               : {sic_threshold:.2f}\n"
                          f"  • Spatial dimension names           : {spatial_dim_names}")
        self.logger.debug("\n[Sea Ice Pressure Computation Steps]\n"
                          f"  1. Apply SIC mask (SIC > {sic_threshold:.2f})\n"
                          f"  2. Multiply by 100 cm/m, sum over spatial dims, then divide by seconds per day")
        mask = SIC > sic_threshold
        DVT  = DVT.where(mask)
        return (DVT*1e2).sum(dim=spatial_dim_names) / 8.64e4 #m/s

    def compute_hemisphere_ice_strength(self, SIC, HI, IS,
                                        spatial_dim_names  : list  = None,
                                        sic_threshold      : float = None,
                                        ice_strength_scale : float = 1e5):
        """
        Compute hemispheric mean sea ice pressure in hectopascals (hPa), based on masked concentration,
        thickness, strength, and area fields.

        This function calculates an area-weighted mean ice thickness (m) by dividing total
        ice volume by total ice-covered area. It then converts internal ice strength (N/m) 
        into an equivalent pressure field (Pa = N/m²), and finally into hectopascals (hPa). 
        Grid cells below a threshold sea ice concentration are masked from the analysis.

        INPUTS:
        -------------------
        SIC               : xr.DataArray; Sea ice concentration (unitless, typically 0–1), per grid cell.
        HI                : xr.DataArray; Sea ice thickness in meters, per grid cell.
        IS                : xr.DataArray; Internal ice strength in N/m, per grid cell.
        A                 : xr.DataArray; Grid cell area in m².
        spatial_dim_names : list, optional; Names of the spatial dimensions to reduce over (e.g., ["nj", "ni"]).
                            If not provided, defaults to values in self.CICE_dict['spatial_dims'].
        sic_threshold     : float, optional; Threshold for masking sea ice concentration (default: self.icon_thresh).

        OUTPUTS:
        -------
        P : xr.DataArray; Sea ice pressure in hectopascals (hPa), with masked cells excluded
            and computed as (IS / mean_thickness) / 100.

        Notes
        -----
        - The output pressure is an approximation derived from CICE's internal strength field and is not a direct model output.
        - Use with caution where sea ice thickness is very small or near-zero.
        """
        sic_threshold     = sic_threshold     if sic_threshold     is not None else self.icon_thresh
        spatial_dim_names = spatial_dim_names if spatial_dim_names is not None else self.CICE_dict['spatial_dims']
        self.logger.info("Computing **INTERNAL ICE PRESSURE** for hemisphere")
        self.logger.debug(f"  • Using SIC threshold               : {sic_threshold:.2f}\n"
                          f"  • Spatial dimension names           : {spatial_dim_names}")
        self.logger.debug("\n[Sea Ice Pressure Computation Steps]\n"
                          f"  1. Apply SIC mask (SIC > {sic_threshold:.2f})\n"
                          f"  2. Divide internal ice strength (N/m = kg*m/s^2) by thickness (m) to get ice pressure (kg/m*s^2 = Pa) then divide by 100 to get hPa")
        mask = SIC > sic_threshold
        HI   = HI.where(mask)
        IS   = IS.where(mask)
        return (IS / HI).sum(dim=spatial_dim_names) / ice_strength_scale

    def compute_sector_ice_area(self, I_mask, area_grid, sector_defs, GI_area=False):
        # Compute ice area per geographic sector and a domain total, from an arbitrary ice mask.
        # This generalises compute_sector_FIA while preserving the original FI-specific wrapper.
        if GI_area:
            self.compute_grounded_iceberg_area()
            GI_ttl_area = self.compute_grounded_iceberg_area() / 1e6
            self.logger.info(f"adding {GI_ttl_area} to ice area computation")
        else:
            GI_ttl_area = 0
        sector_names = list(sector_defs.keys())
        ia_values    = []
        for sec_name in sector_names:
            lon_min, lon_max, lat_min, lat_max = sector_defs[sec_name]["geo_region"]
            sec_mask = I_mask.where(
                (area_grid.lon >= lon_min) & (area_grid.lon <= lon_max) &
                (area_grid.lat >= lat_min) & (area_grid.lat <= lat_max)
            )
            tmp = (sec_mask * area_grid).sum(skipna=True)
            if hasattr(tmp, "compute"):
                ia_val = float(tmp.compute())
            else:
                ia_val = tmp.item() if hasattr(tmp, "item") else float(tmp)
            ia_values.append(ia_val)
        IA_da  = xr.DataArray(ia_values, dims=["sector"], coords={"sector": sector_names})
        IA_tot = sum(ia_values) + GI_ttl_area
        return IA_da, IA_tot

    def compute_sector_FIA(self, FI_mask, area_grid, sector_defs, GI_area=False):
        """
        Backwards-compatible wrapper for fast-ice area (FIA).

        For new workflows, prefer `compute_sector_ice_area`, which is ice-type agnostic.
        """
        return self.compute_sector_ice_area(FI_mask, area_grid, sector_defs, GI_area=GI_area)

    def compute_hemisphere_ice_area(self, SIC, A,
                                    ice_area_scale            = None,
                                    spatial_dim_names         = None,
                                    sic_threshold             = None,
                                    add_grounded_iceberg_area = None,
                                    grounded_iceberg_area     = None,
                                    region                    = None):
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
            
        A : xarray.DataArray
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
        spatial_dim_names         = spatial_dim_names         if spatial_dim_names         is not None else self.CICE_dict['spatial_dims']
        ice_area_scale            = ice_area_scale            if ice_area_scale            is not None else self.FIC_scale
        if region is not None:
            self.logger.info(f"Computing Ice **AREA** for {region} ")
        else:
            self.logger.info("Computing Ice **AREA** ")
        self.logger.debug(f"  • Using SIC threshold               : {sic_threshold:.2f}\n"
                          f"  • Spatial dimension names           : {spatial_dim_names}\n"
                          f"  • Ice area scaling factor           : {ice_area_scale:.2e} (to convert m² → km² or other units)\n"
                          f"  • Include grounded iceberg area     : {add_grounded_iceberg_area}")
        if add_grounded_iceberg_area:
            GI_total_area = (grounded_iceberg_area if grounded_iceberg_area is not None else self.compute_grounded_iceberg_area() )
        else:
            GI_total_area = 0
        self.logger.debug("\n[Sea Ice Area Computation Steps]\n"
                        f"  1. Apply SIC mask (SIC > {sic_threshold:.2f})\n"
                        f"  2. Multiply SIC × area and sum over {spatial_dim_names}\n"
                        f"  3. Add grounded iceberg area (GIA) = {GI_total_area:.2f} m²\n"
                        f"  4. Divide by scale factor {ice_area_scale:.2e} to convert units")
        mask = SIC > sic_threshold
        SIC  = SIC.where(mask)
        A    = A.where(mask)
        IA   = ((SIC * A).sum(dim=spatial_dim_names))
        return (IA + GI_total_area) / ice_area_scale 

    def compute_hemisphere_ice_volume(self, SIC, HI, A,
                                      add_grounded_iceberg_area = None,
                                      ice_volume_scale          = 1e12, 
                                      spatial_dim_names         = None,
                                      sic_threshold             = None,
                                      grounded_iceberg_area     = None):
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
            
        A : xarray.DataArray
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
        add_grounded_iceberg_area = add_grounded_iceberg_area if add_grounded_iceberg_area is not None else self.use_gi
        spatial_dim_names         = spatial_dim_names         if spatial_dim_names         is not None else self.CICE_dict['spatial_dims']
        sic_threshold             = sic_threshold             if sic_threshold             is not None else self.icon_thresh
        self.logger.info("Computing Ice **VOLUME** ")
        self.logger.debug(f"  • Using SIC threshold               : {sic_threshold:.2f}\n"
                          f"  • Spatial dimension names           : {spatial_dim_names}\n"
                          f"  • Ice volume scaling factor         : {ice_volume_scale:.2e} (to convert {1e6:.1e}-km^3)\n"
                          f"  • Include grounded iceberg area     : {add_grounded_iceberg_area}")
        if add_grounded_iceberg_area:
            GI_total_area = (grounded_iceberg_area if grounded_iceberg_area is not None else self.compute_grounded_iceberg_area() )
        else:
            GI_total_area = 0
        self.logger.debug("\n[Sea Ice Volume Computation Steps]\n"
                         f"  1. Apply SIC mask (SIC > {sic_threshold:.2f})\n"
                         f"  2. Multiply SIC × thickness × area and sum over {spatial_dim_names}\n"
                         f"  3. Add grounded iceberg area (GIA) = {GI_total_area:.2f} m²\n"
                         f"  4. Divide by scale factor {ice_volume_scale:.2e} to convert units")
        mask = SIC > sic_threshold
        HI   = HI.where(mask)
        SIC  = SIC.where(mask)
        A    = A.where(mask)
        IV   = (SIC * HI * A).sum(dim=spatial_dim_names)#.chunk({'time': -1}).compute()
        return (IV + GI_total_area) / ice_volume_scale 

    def compute_hemisphere_ice_thickness(self, SIC, HI, A, 
                                         spatial_dim_names = None,
                                         sic_threshold     = None):
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
        self.logger.info("Computing Ice **THICKNESS** ")
        self.logger.debug(f"  • Using SIC threshold               : {sic_threshold:.2f}\n"
                          f"  • Spatial dimension names           : {spatial_dim_names}")
        self.logger.debug("\n[Sea Ice Thickness Computation Steps]\n"
                        f"  1. Apply SIC mask (SIC > {sic_threshold:.2f})\n"
                        f"  2. Compute ice area: Multiply SIC × area and sum over {spatial_dim_names}\n"
                        f"  3. Compute ice volume: Multiply SIC × thickness × area and sum over {spatial_dim_names}\n"
                        f"  4. Divide by area / volume = thickness in metres")
        mask = SIC > sic_threshold
        HI   = HI.where(mask)
        SIC  = SIC.where(mask)
        A    = A.where(mask)
        IA   = (SIC * A).sum(dim=spatial_dim_names)
        IV   = (HI * A).sum(dim=spatial_dim_names)
        return (IV / IA)

    def compute_hemisphere_variable_aggregate(self, da, time_coord_name='time', flat_mean=False, da2=None, scale=None):
        """
        Compute a time-aggregated field from a DataArray across the specified time dimension.

        This method calculates the mean of the input variable over time, with an optional
        normalisation step based on the maximum value across the time dimension.
        Designed for Dask-backed DataArrays to support scalable computation.

        Parameters
        ----------
        da : xarray.DataArray
            Input data array to aggregate over time. Typically a geophysical variable such as
            sea ice thickness (m), strength (N/m), or velocity (m/s), with time as one of the dimensions.
        time_coord_name : str, optional
            Name of the time dimension to aggregate over. Default is 'time'.
        flat_mean : bool, optional
            If True, the data has simple mean over time, Default is False.

        Returns
        -------
        xarray.DataArray
            A time-aggregated version of the input array. If `normalise_max` is True,
            the result represents the normalised mean occurrence across the time dimension;
            otherwise, it is the simple time mean.

        Notes
        -----
        - If the specified `time_coord_name` is not found in the input DataArray,
        the original array is returned with a warning.
        - If `normalise_max` is enabled, the output represents a relative frequency-like
        measure, especially useful for thresholded or binary input fields (e.g., fast ice presence).
        - Assumes that all non-time dimensions are to be preserved in the output.
        """
        self.logger.info("Computing Ice **AGGREGATE OVER TIME** (percentage of occurence)")
        if time_coord_name not in da.dims:
            self.logger.warning(f"No time dimension '{time_coord_name}' found. Returning original array.")
            return da
        time_len = da[time_coord_name].sizes.get(time_coord_name, 1)
        if flat_mean:
            return da.mean(dim=time_coord_name)
        elif da2 and scale:
            return (da/da2).sum(dim=time_coord_name) / scale
        else:
            return da.sum(dim=time_coord_name) / time_len

    # SEASONAL METRICS
    @staticmethod
    def _ensure_odd(n: int) -> int:
        return n if n % 2 == 1 else n + 1

    @staticmethod
    def _safe_slope(y: np.ndarray, t_days: np.ndarray, scaling: float) -> float | None:
        """Least-squares slope without polyfit overhead; returns None if flat/NaN."""
        mask = np.isfinite(y) & np.isfinite(t_days)
        if mask.sum() < 2:
            return None
        y = y[mask]
        t = t_days[mask]
        t = t - t[0]  # improves conditioning
        t_mean = t.mean()
        y_mean = y.mean()
        denom = np.sum((t - t_mean) ** 2)
        if denom == 0:
            return None
        slope = np.sum((t - t_mean) * (y - y_mean)) / denom
        return float(slope * scaling)

    @staticmethod
    def _year_index_slices(times: np.ndarray) -> dict[int, np.ndarray]:
        years = pd.DatetimeIndex(times).year.values
        uniq = np.unique(years)
        return {int(y): np.where(years == y)[0] for y in uniq}

    @staticmethod
    def _drop_leap_day(times: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        # Remove Feb 29 for DOY-consistent filters
        di = pd.DatetimeIndex(times)
        keep = ~((di.month == 2) & (di.day == 29))
        return times[keep], y[keep]

    @staticmethod
    def _summary(x):
        x = [v for v in x if v is not None]
        return (float(np.mean(x)), float(np.std(x))) if len(x) else (None, None)

    def compute_seasonal_statistics(self, da,
                                    stat_name           : str = 'FIA',
                                    window              : int = 15,
                                    polyorder           : int = 2,
                                    min_onset_doy       : int = 50,
                                    min_retreat_doy     : int = 240,
                                    growth_range        : tuple[int, int] = (71, 273),
                                    retreat_early_range : tuple[int, int] = (273, 330),
                                    retreat_late_range  : tuple[tuple[int, int], tuple[int, int]] = ((331, 365), (1, 71)),
                                    scaling_factor      : float = 1e6,
                                    drop_leap_day       : bool = True,
                                    return_per_year     : bool = False):
        """
        Fast version that loads once, computes derivatives once per year, and uses NumPy slicing.

        Parameters
        ----------
        da : xr.DataArray  # 1D, daily, with 'time'
        window, polyorder, min_onset_doy, min_retreat_doy : as before
        growth_range, retreat_early_range, retreat_late_range : as before
        scaling_factor : float  # multiplies slopes for nicer units
        drop_leap_day : bool    # remove Feb 29 to keep DOY logic simple
        return_per_year : bool  # also return a per-year table if True

        Returns
        -------
        summary : dict[str, float | None]
        (optional) per_year : pandas.DataFrame  # one row per year (except last, by design)
        """
        from scipy.signal import savgol_filter
        from tqdm         import tqdm
        if 'time' not in da.dims:
            raise ValueError("Input DataArray must have a 'time' dimension")
        self.logger.info(f"computing SEASONAL STATISTICS for {stat_name} ... ")
        # 1) Load once (small 1D array) and optionally drop leap days
        self.logger.debug("   1. loading array ")
        da = da.load()  # safe for 1D daily series
        times = da.time.values
        y_full = da.values.astype(np.float32)
        self.logger.debug("   2. drop leap-days")
        if drop_leap_day:
            times, y_full = self._drop_leap_day(times, y_full)
        self.logger.debug("   3. compute DOY indeces once")
        doy_full = pd.DatetimeIndex(times).dayofyear.values
        idx_by_year = self._year_index_slices(times)
        years_sorted = sorted(idx_by_year.keys())
        self.logger.debug("   4. prepare statistical lists")
        per_year_rows = []
        max_list, min_list = [], []
        grow_list, ret_list = [], []
        ret_early_list, ret_late_list = [], []
        min_doy_list, max_doy_list = [], []
        onset_list = [] 
        duration_list = []
        self.logger.debug("   5. Ensure Savitzky–Golay (savgol) settings are valid")
        window = self._ensure_odd(window)
        if window < polyorder + 2:
            window = polyorder + 3 if (polyorder + 3) % 2 == 1 else polyorder + 4
        self.logger.debug("   6. compute derivatives once per-year")
        for i, year in enumerate(years_sorted[:-1]):
            self.logger.debug(f"     ....{year}")
            idx_this = idx_by_year[year]
            idx_next = idx_by_year[years_sorted[i + 1]]
            t_this = times[idx_this]
            y_this = y_full[idx_this]
            doy_this = doy_full[idx_this]
            # Skip tiny years
            if y_this.size < max(window, polyorder + 2):
                continue
            # Derivatives once per year
            try:
                dy_dt = savgol_filter(y_this, window_length=window, polyorder=polyorder, deriv=1, mode='interp')
                d2y_dt2 = savgol_filter(y_this, window_length=window, polyorder=polyorder, deriv=2, mode='interp')
            except ValueError:
                # window too long for this year, skip
                continue
            # --- Max/Min & their DOY ---
            if np.all(~np.isfinite(y_this)):
                # skip year with all NaN
                continue
            with np.errstate(invalid='ignore'):
                max_idx = int(np.nanargmax(y_this))
                min_idx = int(np.nanargmin(y_this))
            max_val = float(y_this[max_idx])
            min_val = float(y_this[min_idx])
            max_doy = int(doy_this[max_idx])
            min_doy = int(doy_this[min_idx])
            max_list.append(max_val); min_list.append(min_val)
            max_doy_list.append(max_doy); min_doy_list.append(min_doy)
            # --- Onset via derivative (first positive slope after min_onset_doy) ---
            onset_candidates = np.flatnonzero((dy_dt > 0) & (doy_this >= min_onset_doy))
            onset_doy = int(doy_this[onset_candidates[0]]) if onset_candidates.size else None
            # --- Retreat via inflection (last negative curvature after min_retreat_doy) ---
            retreat_candidates = np.flatnonzero((d2y_dt2 < 0) & (doy_this >= min_retreat_doy))
            retreat_doy = int(doy_this[retreat_candidates[-1]]) if retreat_candidates.size else None
            # --- Duration (onset → retreat, handle wrap) ---
            if onset_doy is not None and retreat_doy is not None:
                duration = (retreat_doy - onset_doy) if (retreat_doy > onset_doy) else (365 - onset_doy + retreat_doy)
                duration_list.append(int(duration))
                onset_list.append(int(onset_doy))
            # --- Slopes (growth / retreat early / retreat late / overall retreat) ---
            # Build masks once (NumPy, no xarray .where)
            doy_mask = doy_this  # alias
            g0, g1 = growth_range
            re0, re1 = retreat_early_range
            rl0a, rl1a = retreat_late_range[0]
            rl0b, rl1b = retreat_late_range[1]
            # growth season (this year)
            g_mask = (doy_mask >= g0) & (doy_mask <= g1)
            y_grow = y_this[g_mask]
            t_grow = (t_this[g_mask] - t_this[g_mask][0]) / np.timedelta64(1, 'D') if y_grow.size else np.array([])
            slope_grow = self._safe_slope(y_grow, t_grow, scaling_factor)
            if slope_grow is not None:
                grow_list.append(slope_grow)
            # retreat early (this year)
            e_mask = (doy_mask >= re0) & (doy_mask <= re1)
            y_early = y_this[e_mask]
            t_early = (t_this[e_mask] - t_this[e_mask][0]) / np.timedelta64(1, 'D') if y_early.size else np.array([])
            slope_early = self._safe_slope(y_early, t_early, scaling_factor)
            if slope_early is not None:
                ret_early_list.append(slope_early)
            # retreat late (next year)
            t_next = times[idx_next]
            y_next = y_full[idx_next]
            doy_next = doy_full[idx_next]
            l_mask_a = (doy_mask >= rl0a) & (doy_mask <= rl1a)
            l_mask_b = (doy_next >= rl0b) & (doy_next <= rl1b)
            y_late = np.concatenate([y_this[l_mask_a], y_next[l_mask_b]])
            t_late = np.concatenate([(t_this[l_mask_a] - t_this[l_mask_a][0]) / np.timedelta64(1, 'D') if l_mask_a.any() else np.array([]),
                                      # Continue time axis across year boundary
                                     ((t_next[l_mask_b] - t_this[l_mask_a][0]) / np.timedelta64(1, 'D')) if l_mask_b.any() and l_mask_a.any() else
                                     ((t_next[l_mask_b] - t_next[l_mask_b][0]) / np.timedelta64(1, 'D')) if l_mask_b.any() else np.array([])]) if (l_mask_a.any() or l_mask_b.any()) else np.array([])
            slope_late = self._safe_slope(y_late, t_late, scaling_factor)
            if slope_late is not None:
                ret_late_list.append(slope_late)
            # overall retreat (early + late concatenated)
            if y_early.size or y_late.size:
                y_ret = np.concatenate([y_early, y_late])
                # stitch time so it's strictly increasing
                t_ret = np.concatenate([t_early if y_early.size else np.array([]), (t_late + (t_early[-1] + 1) if y_early.size and y_late.size else t_late)])
                slope_ret = self._safe_slope(y_ret, t_ret, scaling_factor)
                if slope_ret is not None:
                    # Keep sign convention negative (ablation)
                    ret_list.append(-abs(slope_ret))
            if return_per_year:
                per_year_rows.append({'year'                : year,
                                      'max'                 : max_val,   'min'         : min_val,
                                      'doy_max'             : max_doy,   'doy_min'     : min_doy,
                                      'onset_doy'           : onset_doy, 'retreat_doy' : retreat_doy,
                                      'duration_days'       : duration if (onset_doy is not None and retreat_doy is not None) else None,
                                      'growth_slope'        : slope_grow,
                                      'retreat_early_slope' : slope_early,
                                      'retreat_late_slope'  : slope_late,
                                      'retreat_slope'       : slope_ret if slope_ret is not None else None,})
        self.logger.debug("   7. compile summary and return")
        summary = {"Maximum Mean"       : self._summary(max_list)[0],       "Maximum Std"       : self._summary(max_list)[1],
                   "Minimum Mean"       : self._summary(min_list)[0],       "Minimum Std"       : self._summary(min_list)[1],
                   "Growth Mean"        : self._summary(grow_list)[0],      "Growth Std"        : self._summary(grow_list)[1],
                   "Retreat Mean"       : self._summary(ret_list)[0],       "Retreat Std"       : self._summary(ret_list)[1],
                   "Retreat Early Mean" : self._summary(ret_early_list)[0], "Retreat Early Std" : self._summary(ret_early_list)[1],
                   "Retreat Late Mean"  : self._summary(ret_late_list)[0],  "Retreat Late Std"  : self._summary(ret_late_list)[1],
                   "Duration-days Mean" : self._summary(duration_list)[0],  "Duration-days Std" : self._summary(duration_list)[1],
                   "DOY Min Mean"       : self._summary(min_doy_list)[0],   "DOY Min Std"       : self._summary(min_doy_list)[1],
                   "DOY Max Mean"       : self._summary(max_doy_list)[0],   "DOY Max Std"       : self._summary(max_doy_list)[1],
                   "DOY Onset Mean"     : self._summary(onset_list)[0],     "DOY Onset Std"     : self._summary(onset_list)[1],}
        if return_per_year:
            return summary, pd.DataFrame(per_year_rows).set_index('year')
        return summary

    # PERSISTENCE METRICS
    def persistence_stability_index(self, I_mask, A,
                                    persistence_threshold: float = 0.8,
                                    winter_months        : tuple = (5, 6, 7, 8, 9, 10),
                                    area_scale           : float = 1e9) -> dict:
        """
        Compute FIPSI (Fast-Ice Persistence Stability Index) from a binary fast-ice mask.

        Definition
        ----------
        FIPSI = (area of grid cells whose WINTER fast-ice persistence ≥ threshold)
                / (area of grid cells that have fast ice at least once in WINTER),
        where persistence is the fraction of winter days (e.g., May–Oct) that a cell is fast.

        Parameters
        ----------
        I_mask : xarray.DataArray or xarray.Dataset
            Binary fast-ice mask (0/1) with a 'time' dimension; if a Dataset is given,
            a variable named 'FI_mask' is used.
        A : xarray.DataArray or xarray.Dataset
            Grid-cell area on the same grid; if a Dataset is given, the first of
            {'tarea','area','TAREA'} is used. If A has a 'time' dim, the first slice is used.
        persistence_threshold : float, optional
            Minimum winter persistence to classify a cell as persistent (default 0.8).
        winter_months : tuple[int], optional
            Months to include for the winter persistence (default May–Oct).

        Returns
        -------
        dict:
            {'persistence_stability_index': float,   # FIPSI in [0,1] or NaN if undefined
             'area_persistent_winter': float,        # area of persistent winter FI (units of A)
             'area_ever_winter': float               # area with any winter FI (units of A)}

        Notes
        -----
        • Uses the class's `spatial_dims` to sum areas.
        • Assumes `I_mask` is truly binary (0/1). For classification outputs that are boolean, they will be cast to float internally.
        • Designed to be Dask-friendly: only two area reductions are computed eagerly.
        """
        self.logger.info("Computing FIPSI: persistent-winter area / ever-winter area")
        I_mask = self._as_da_mask(I_mask).astype("float32")
        A      = self._as_da_area(A)
        # Ensure A has only spatial dims
        if "time" in A.dims:
            A = A.isel(time=0)
            if "time" in A.coords:
                # safer than drop_vars which can produce a Dataset
                A = A.drop_vars("time", errors="ignore")
            A = A.squeeze(drop=True)
        spat_dims = self.CICE_dict["spatial_dims"]
        # restrict to winter months and form persistence (fraction in [0,1])
        in_winter = I_mask["time"].dt.month.isin(list(winter_months))
        FI_winter = I_mask.where(in_winter)
        prstnc    = FI_winter.mean(dim="time", skipna=True)
        # masks for footprints
        persistent_mask = prstnc >= float(persistence_threshold)       # persistent winter FI
        ever_mask       = FI_winter.max(dim="time", skipna=True) > 0   # ever-winter FI
        # # area-weighted sums over spatial dims
        prstnc_A = (A.where(persistent_mask)).sum(dim=spat_dims, skipna=True)
        ttl_A    = (A.where(ever_mask)).sum(dim=spat_dims, skipna=True)
        with dask.config.set(scheduler="threads"):   # or "synchronous"
            prstnc_A = prstnc_A.compute()
            ttl_A    = ttl_A.compute()
        prstnc_A = self._to_float_scalar(prstnc_A)
        ttl_A    = self._to_float_scalar(ttl_A)
        psi      = (prstnc_A / ttl_A) if (np.isfinite(ttl_A) and ttl_A > 0.0) else np.nan
        return {"persistence_stability_index": float(psi),
                "area_persistent_winter"     : prstnc_A/area_scale,
                "area_ever_winter"           : ttl_A/area_scale}

    def _prepare_BAS_coast(self,
                           path_coast_shape  : str   =  None,
                           crs_out           : str   = "EPSG:3031",
                           target_spacing_km : float = 1.0):
        """
        Load BAS high-res coastline polygons and build a KD-Tree of densified coastline points in a metric CRS (default EPSG:3031).
        Caches results on `self` to avoid repeated work.

        source: https://add-catalogue.data.bas.ac.uk/records/9b4fab56-4999-4fe1-b17f-1466e41151c4.html

        Parameters
        ----------
        shp_path : str
            Path to the polygon shapefile (e.g., ".../add_coastline_high_res_polygon_v7_9.shp").
        crs_out : str
            Target projected CRS for metric distances (e.g., EPSG:3031).
        target_spacing_km : float
            Approx spacing along the coastline (km) for densification.

        Sets
        ----
        self._coast_kdtree  : scipy.spatial.cKDTree instance in projected coords
        self._coast_xy_proj : (x_coast, y_coast) ndarray tuple (meters)
        self._coast_crs     : str of the projected CRS
        """
        import geopandas      as gpd
        from scipy.spatial    import cKDTree
        from shapely.geometry import LineString, Polygon, MultiPolygon
        from shapely.ops      import unary_union
        P_shp = path_coast_shape if path_coast_shape is not None else self.BAS_dict["P_Ant_Cstln"]
        self.logger.info(f"Loading coastline polygons: {P_shp}")
        gdf = gpd.read_file(P_shp)
        # Reproject to metric CRS (EPSG:3031 by default)
        if gdf.crs is None:
            self.logger.warning("Input CRS is None; assuming EPSG:4326 (lon/lat).")
            gdf = gdf.set_crs("EPSG:4326", allow_override=True)
        gdf = gdf.to_crs(crs_out)
        # Merge polygons to avoid duplicate edges (optional but helpful)
        self.logger.info("Merging polygons and extracting boundaries...")
        geom = unary_union(gdf.geometry.values)  # MultiPolygon | Polygon
        if isinstance(geom, Polygon):
            polygons = [geom]
        elif isinstance(geom, MultiPolygon):
            polygons = list(geom.geoms)
        else:
            # In case the file already has lines, just take them
            polygons = []
        # Convert polygon exteriors to lines (ignore interiors/holes for coastline)
        lines = []
        for poly in polygons:
            try:
                lines.append(LineString(poly.exterior.coords))
            except Exception:
                pass
        # If there are already line features in the shapefile, add them too
        lines.extend([g for g in gdf.geometry.values if g.geom_type.lower().endswith("linestring")])
        if not lines:
            raise ValueError("No coastline linework could be derived from the polygons/lines in the shapefile.")
        # Densify lines at ~target spacing
        spacing_m = max(100.0, float(target_spacing_km) * 1000.0)  # clip at >=100 m
        pts_x, pts_y = [], []
        for ln in lines:
            try:
                length = ln.length
                if not np.isfinite(length) or length <= 0:
                    continue
                n = max(2, int(np.ceil(length / spacing_m)) + 1)
                # sample [0,1] along the line
                for f in np.linspace(0.0, 1.0, n):
                    p = ln.interpolate(f, normalized=True)
                    xy = p.coords[0]
                    pts_x.append(xy[0])
                    pts_y.append(xy[1])
            except Exception:
                continue
        x_coast = np.asarray(pts_x, dtype=float)
        y_coast = np.asarray(pts_y, dtype=float)
        if x_coast.size == 0:
            raise RuntimeError("Coastline densification produced zero points.")
        self.logger.info(f"Built coastline point set: {x_coast.size:,} pts @ ~{target_spacing_km} km spacing.")
        kdt = cKDTree(np.column_stack([x_coast, y_coast]))
        # Cache
        self._coast_kdtree  = kdt
        self._coast_xy_proj = (x_coast, y_coast)
        self._coast_crs     = crs_out

    def persistence_ice_distance_mean_max(self, ice_prstnc,
                                          persistence_threshold : float = 0.8,
                                          path_coast_shape      : str = None,
                                          crs_out               : str = "EPSG:3031"):
        """
        Mean/max distance from the USNIC/ADD coastline for persistent fast ice.
        Uses a densified coastline KD-Tree in a metric CRS (default EPSG:3031).

        Assumes `ice_prstnc` is a persistence field in [0..1] computed from a
        binary fast-ice mask over the austral-winter window (e.g., May–Oct).
        """
        from pyproj import Transformer
        P_shp = path_coast_shape if path_coast_shape is not None else self.BAS_dict["P_Ant_Cstln"]
        spat_dims = self.CICE_dict["spatial_dims"]
        self.logger.info("Computing persistence distance (USNIC coastline KD-Tree)")
        with dask.config.set(scheduler="threads"):
            ice_prstnc = ice_prstnc.compute()  # instead of .load()
        if not getattr(self, "bgrid_loaded", False):
            self.load_cice_grid(slice_hem=True)
        Gt       = self.G_t
        lon_name = 'TLON' if 'TLON' in Gt else ('lon' if 'lon' in Gt else None)
        lat_name = 'TLAT' if 'TLAT' in Gt else ('lat' if 'lat' in Gt else None)
        if lon_name is None or lat_name is None:
            raise KeyError("Could not find lon/lat fields (TLON/TLAT or lon/lat) in self.G_t.")
        grid_shape = tuple(ice_prstnc.sizes[d] for d in spat_dims)
        if grid_shape != Gt[lon_name].shape:
            self.logger.warning("Grid mismatch detected; reloading b-grid with slice_hem=True")
            self.load_cice_grid(slice_hem=True)
            Gt = self.G_t
            if tuple(ice_prstnc.sizes[d] for d in spat_dims) != Gt[lon_name].shape:
                raise ValueError("ice_prstnc and grid coordinate shapes still mismatch after reload.")
        # Prepare coastline KD-tree once (unchanged)
        need_kdtree = (not hasattr(self, "_coast_kdtree")) or (getattr(self, "_coast_crs", None) != crs_out)
        if need_kdtree:
            self._prepare_BAS_coast(P_shp, crs_out=crs_out, target_spacing_km=1.0)
        # Build persistent mask from fully-realized array
        prst = (ice_prstnc >= float(persistence_threshold)).values
        if not np.any(prst):
            self.logger.warning("No persistent fast-ice cells found at this threshold.")
            return {"persistence_mean_distance": float('nan'),
                    "persistence_max_distance":  float('nan')}
        lon = np.asarray(Gt[lon_name].values)
        lat = np.asarray(Gt[lat_name].values)
        transformer = Transformer.from_crs("EPSG:4326", crs_out, always_xy=True)
        x_all, y_all = transformer.transform(lon, lat)
        iy, ix = np.where(prst)  # spat_dims order (nj, ni)
        x_p = x_all[iy, ix]
        y_p = y_all[iy, ix]
        if x_p.size == 0:
            return {"persistence_mean_distance": float('nan'),
                    "persistence_max_distance":  float('nan')}
        d_m, _ = self._coast_kdtree.query(np.column_stack([x_p, y_p]), k=1)
        d_km = d_m / 1000.0
        return {"persistence_mean_distance": float(np.nanmean(d_km)),
                "persistence_max_distance":  float(np.nanmax(d_km))}

    # INTER-COMPARISONS
    @staticmethod
    def _skill_stats(x, y):
        from sklearn.metrics import mean_squared_error, mean_absolute_error
        from scipy.stats     import pearsonr
        valid = (~np.isnan(x)) & (~np.isnan(y))
        x, y  = x[valid], y[valid]
        return {"Bias"     : np.mean(x - y),
                "RMSE"     : np.sqrt(mean_squared_error(y, x)),
                "MAE"      : mean_absolute_error(y, x),
                "Corr"     : pearsonr(x, y)[0],
                "SD_Model" : np.std(x),
                "SD_Obs"   : np.std(y)}

    @staticmethod
    def _safe_load_array_to_memory(x):
        # --- Force local computation to avoid serializing Dask HLGs to the cluster ---
        # Keep memory light and avoid shipping graphs to distributed workers
        try:
            return x.astype("float32").compute(scheduler="single-threaded")
        except Exception:
            return x.compute(scheduler="single-threaded") if hasattr(x, "compute") else x
        
    def compute_skill_statistics(self, mod, obs, min_points=100):
        """
        Compute statistical skill metrics between model and observational time series.
        If time alignment fails and `climatology_if_mismatch` is True, computes climatological
        skill instead (e.g., using day-of-year climatologies).

        Parameters
        ----------
        mod : xarray.DataArray; Time series of modelled values.
        obs : xarray.DataArray; Time series of observed values (must be aligned with model, or climatology is used).

        Returns
        -------
        dict; Dictionary of skill metrics.
        """
        mod = self._safe_load_array_to_memory(mod)
        obs = self._safe_load_array_to_memory(obs)
        if hasattr(mod, 'time') and hasattr(obs, 'time'):
            dt_mod  = mod.time.values
            dt_obs  = obs.time.values
            dt_xsct = np.intersect1d(dt_mod, dt_obs)
            if len(dt_xsct) > min_points:
                self.logger.warning(f"more than {min_points} date-times intersected, therefore do stats on time series")
                mod_data = mod.sel(time=dt_xsct).values
                obs_data = obs.sel(time=dt_xsct).values
                return self._skill_stats(mod_data, obs_data)
            else:
                self.logger.warning("date-time mismatch detected — falling back to climatology comparison.")
                mod_clim = mod.groupby('time.dayofyear').mean('time')
                obs_clim = obs.groupby('time.dayofyear').mean('time')
                doy_xsct = np.intersect1d(mod_clim.dayofyear.values, obs_clim.dayofyear.values)
                mod_data = mod_clim.sel(dayofyear=doy_xsct).values
                obs_data = obs_clim.sel(dayofyear=doy_xsct).values
                return self._skill_stats(mod_data, obs_data)
        self.logger.warning("Time dimension not found or could not align. Returning empty stats.")
        return {"Bias": np.nan, "RMSE": np.nan, "MAE": np.nan, "Corr": np.nan,
                "SD_Model": np.nan, "SD_Obs": np.nan}

    def antarctic_year_labels(self, time_da: xr.DataArray) -> xr.DataArray:
        """Return an 'ayear' label for Jul→Jun years (e.g., 2000.5 dates -> ayear=2000)."""
        year  = time_da.dt.year
        month = time_da.dt.month
        ayear = xr.where(month <= 6, year - 1, year)
        return ayear.rename("ayear")

    def extrema_means_antarctic_year(self, da_ts: xr.DataArray):
        """
        da_ts: 1-D (time,) fast-ice area time series (km^2).
        Returns dict with:
        - max_mean: mean of per-AY maxima
        - min_mean: mean of per-AY minima
        - max_ts  : per-AY maxima time series (indexed by 'ayear')
        - min_ts  : per-AY minima time series (indexed by 'ayear')
        """
        # for tiny 1D series, it's fine to fully chunk or even compute up front
        # da_ts = da_ts.chunk({"time": -1})
        ay = self.antarctic_year_labels(da_ts["time"])
        g  = da_ts.groupby(ay)
        max_ts = g.max(skipna=True)                  # dims: ayear
        min_ts = g.min(skipna=True)
        return {
            "max_mean": max_ts.mean("ayear", skipna=True),
            "min_mean": min_ts.mean("ayear", skipna=True),
            "max_ts": max_ts,
            "min_ts": min_ts,
        }
