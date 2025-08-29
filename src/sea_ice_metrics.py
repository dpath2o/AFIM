import gc, re, dask
import xarray      as xr
import numpy       as np
import pandas      as pd
from pathlib       import Path

class SeaIceMetrics:

    def __init__(self, **kwargs):
        return

    def fast_ice_metrics_wrapper(self,
                                 sim_name       = None,
                                 dt0_str        = None,
                                 dtN_str        = None,
                                 ice_type       = None,
                                 ispd_thresh    = None,
                                 cice_vars      = None,
                                 drop_vars      = None,
                                 drop_coords    = None,
                                 ice_area_scale = None,
                                 overwrite_zarr = None,
                                 overwrite_png  = None):
        """
        Compute, load, and plot fast ice metrics (FIA and FIP) for a given simulation.

        This method orchestrates the full pipeline for computing fast ice area (FIA)
        and fast ice persistence (FIP) from CICE model output, comparing them with
        observational data, and saving output figures and Zarr files.

        Parameters
        ----------
        sim_name : str, optional
            Simulation name. Defaults to `self.sim_name`.
        dt0_str, dtN_str : str, optional
            Start and end dates (inclusive), formatted as 'YYYY-MM-DD'.
            Defaults to `self.dt0_str` and `self.dtN_str`.
        ice_type : str, optional
            Type of fast ice mask to use (e.g., 'FI_BT'). Defaults to `self.ice_type`.
        ispd_thresh : float, optional
            Ice speed threshold for defining fast ice (in m/s).
        cice_vars : list of str, optional
            CICE variables to load. Default: ['aice', 'hi', 'strength', 'tarea']
        drop_vars, drop_coords : list of str, optional
            Variables and coordinates to drop from loaded datasets.
        ice_area_scale : float, optional
            Scale factor for converting grid cell area. Default from `self.FIC_scale`.
        overwrite_zarr, overwrite_png : bool, optional
            If True, overwrite saved Zarr and PNG files respectively.

        Output
        ------
        - Saves Zarr and CSV metric files.
        - Generates faceted FIA+FIP PNG figure for three variants of ice mask.

        Notes
        -----
        - Evaluates raw, rolling, and binary variants of the input ice type.
        - Mask names are assumed to follow CICE convention with suffixes.
        - See Also: `compute_sea_ice_metrics`, `plot_FIA_FIP_faceted`
        """
        sim_name       = sim_name or self.sim_name
        ice_type       = ice_type or self.ice_type
        dt0_str        = dt0_str or self.dt0_str
        dtN_str        = dtN_str or self.dtN_str
        ispd_thresh    = ispd_thresh    if ispd_thresh is not None else self.ispd_thresh
        ice_area_scale = ice_area_scale if ice_area_scale is not None else self.FIC_scale
        cice_vars      = cice_vars      if cice_vars is not None else ['aice', 'hi', 'strength', 'tarea']
        drop_vars      = drop_vars      if drop_vars is not None else self.CICE_dict["drop_vars"]
        drop_coords    = drop_coords    if drop_coords is not None else self.CICE_dict["drop_coords"]
        overwrite_zarr = overwrite_zarr if overwrite_zarr is not None else self.overwrite_zarr_group
        overwrite_png  = overwrite_png  if overwrite_png is not None else self.ow_fig
        # dependent variables
        ispd_thresh_str = f"{ispd_thresh:.1e}".replace("e-0", "e-")
        dt_range_str    = f"{dt0_str[:4]}-{dtN_str[:4]}"
        ice_type_clean  = re.sub(r'(_roll|_bin)?(_BT|_B|_Ta|_Tx)?$', '', ice_type)
        ice_types       = [ice_type, f"{ice_type}_roll", f"{ice_type}_bin"]
        # loop over three ice_types
        for i_type in ice_types:
            P_METS    = Path(self.D_metrics, f"{i_type}_mets.zarr")
            mask_name = f"{ice_type_clean}_mask"
            if not P_METS.exists() or overwrite_zarr:
                self.logger.info(f"Computing metrics: {P_METS}")
                roll   = i_type.endswith("_roll")
                I_mask = self.load_classified_ice(rolling     = roll,
                                                  ice_type    = ice_type,
                                                  ispd_thresh = ispd_thresh,
                                                  dt0_str     = dt0_str,
                                                  dtN_str     = dtN_str,
                                                  variables   = mask_name)[mask_name]
                CICE_SO = self.load_cice_zarr(slice_hem=True, variables=cice_vars)
                if i_type.endswith("_bin") and ice_type_clean == "FI":
                    FI_bin       = self.binary_days_fast_ice_classification(I_mask).chunk(self.CICE_dict["FI_chunks"])
                    I            = CICE_SO.where(FI_bin)
                    I[mask_name] = FI_bin
                else:
                    I = CICE_SO.where(I_mask)
                I      = I.chunk(self.CICE_dict["FI_chunks"])
                I_dict = {mask_name  : I_mask.persist(),
                          'aice'     : I['aice'].drop_vars(drop_coords).persist(),
                          'hi'       : I['hi'].drop_vars(drop_coords).persist(),
                          'strength' : I['strength'].drop_vars(drop_coords).persist(),
                          'dvidtt'   : I['dvidtt'].drop_vars(drop_coords).persist(),
                          'tarea'    : CICE_SO['tarea'].isel(time=0).drop_vars(drop_coords).persist()}
                self.compute_sea_ice_metrics(I_dict,
                                            ice_type       = ice_type_clean,
                                            dt0_str        = dt0_str,
                                            dtN_str        = dtN_str,
                                            P_mets_zarr    = P_METS,
                                            ice_area_scale = ice_area_scale)
        # plot fast ice area and persistence
        self.pygmt_FIA_FIP_panel(sim_name      = sim_name,
                                 ispd_thresh   = ispd_thresh,
                                 overwrite_fig = overwrite_png,
                                 plot_GI       = self.use_gi)

    @staticmethod
    def _clean_zarr_chunks(ds):
        for var in ds.variables:
            if 'chunks' in ds[var].encoding:
                del ds[var].encoding['chunks']
        return ds.chunk(None)

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
            'aice', 'hi', 'tarea', and a corresponding mask ('FI_mask' or 'PI_mask').
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
        if ice_type_clean not in ("FI", "PI"):
            raise ValueError(f"Unexpected cleaned ice_type: {ice_type_clean}")
        self.logger.info(f"results will be written to {P_mets_zarr}")
        mask_name = f"{ice_type_clean}_mask"
        I_mask    = da_dict[mask_name]
        I_C       = da_dict['aice']
        I_T       = da_dict['hi']
        I_S       = da_dict['strength']
        I_TVT     = da_dict['dvidtt']
        I_MVT     = da_dict['dvidtd']
        I_TAT     = da_dict['daidtt']
        I_MAT     = da_dict['daidtd']
        A         = da_dict['tarea']
        # --- Time Series Metrics ---
        IA   = self.compute_hemisphere_ice_area(I_C, A, ice_area_scale=ice_area_scale)
        IV   = self.compute_hemisphere_ice_volume(I_C, I_T, A)
        IT   = self.compute_hemisphere_ice_thickness(I_C, I_T, A)
        IS   = self.compute_hemisphere_ice_strength(I_C, I_T, I_S)
        ITVR = self.compute_hemisphere_ice_volume_rate(I_C, I_TVT)
        IMVR = self.compute_hemisphere_ice_volume_rate(I_C, I_MVT)
        ITAR = self.compute_hemisphere_ice_area_rate(I_TAT, IA, A)
        IMAR = self.compute_hemisphere_ice_area_rate(I_MAT, IA, A)
        IP   = self.compute_hemisphere_variable_aggregate(I_C)
        self.logger.info("computing **ICE THICKNESS TEMPORAL-MEAN**")
        IHI = I_T.mean(dim=self.CICE_dict["time_dim"])
        self.logger.info("computing **ICE STRENGTH TEMPORAL-SUM**; units mPa")
        IST = (I_S/I_T).sum(dim=self.CICE_dict["time_dim"]) / 1e6
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
                f"{ice_type_clean}S"      : IS.load(),      #1D
                f"{ice_type_clean}TVR"    : ITVR.load(),    #1D
                f"{ice_type_clean}MVR"    : IMVR.load(),    #1D
                f"{ice_type_clean}TAR"    : ITAR.load(),    #1D
                f"{ice_type_clean}MAR"    : IMAR.load(),    #1D
                f"{ice_type_clean}P"      : IP.load(),      #2D
                f"{ice_type_clean}HI"     : IHI.load(),     #2D
                f"{ice_type_clean}ST"     : IST.load(),     #2D
                f"{ice_type_clean}TVR_YR" : ITVR_YR.load(), #2D
                f"{ice_type_clean}MVR_YR" : IMVR_YR.load(), #2D
                f"{ice_type_clean}TAR_YR" : ITAR_YR.load(), #2D
                f"{ice_type_clean}MAR_YR" : IMAR_YR.load()} #2D
        # --- Skill Statistics ---
        try:
            if ice_type_clean == "FI":
                self.logger.info(f"loading FI obs: {self.AF_FI_dict['P_AF2020_FIA']}")
                IA_obs = xr.open_dataset(self.AF_FI_dict['P_AF2020_FIA'], engine="netcdf4")["AF2020"]
            elif ice_type_clean == "PI":
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
        try:
            IP_stab = self.persistence_stability_index(IP, A)
        except Exception as e:
            self.logger.warning(f"persistence_stability_index failed: {e}")
            IP_stab = {}
        try:
            IP_dist = self.persistence_ice_distance_mean_max(IP)
        except Exception as e:
            self.logger.warning(f"persistence_ice_distance_mean_max failed: {e}")
            IP_dist = {}
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
            DS_METS.to_zarr(P_mets_zarr, mode="w", consolidated=True)
            self.logger.info(f"Metrics written to {P_mets_zarr}")
        return DS_METS

    def load_computed_metrics(self,
                              spatial_grid_type : str   = "bin", #bin, roll, None (which is 'daily')
                              ice_type          : str   = None,
                              ispd_thresh       : str   = None):
        ice_type    = ice_type       if ice_type    is not None else self.ice_type
        ispd_thresh = ispd_thresh    if ispd_thresh is not None else self.ispd_thresh
        ispd_thresh_str = f"{ispd_thresh:.1e}".replace("e-0", "e-")
        D_mets      = Path(self.D_zarr, f"ispd_thresh_{ispd_thresh_str}", "metrics")
        P_mets      = Path(D_mets, f"{ice_type}_{spatial_grid_type}_mets.zarr") if spatial_grid_type is not None else Path(D_mets, f"{ice_type}_mets.zarr")
        return xr.open_dataset(P_mets)

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
            self.logger.info("Computing Ice **AREA** for Hemisphere")
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
        IA   = ((SIC * A).sum(dim=spatial_dim_names))#.chunk({'time': -1}).compute()
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
        self.logger.info("Computing Ice **VOLUME** for Hemisphere")
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
        self.logger.info("Computing Ice **THICKNESS** for Hemisphere")
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
                                    stat_name : str = 'FIA',
                                    window: int = 15,
                                    polyorder: int = 2,
                                    min_onset_doy: int = 50,
                                    min_retreat_doy: int = 240,
                                    growth_range: tuple[int, int] = (71, 273),
                                    retreat_early_range: tuple[int, int] = (273, 330),
                                    retreat_late_range: tuple[tuple[int, int], tuple[int, int]] = ((331, 365), (1, 71)),
                                    scaling_factor: float = 1e6,
                                    drop_leap_day: bool = True,
                                    return_per_year: bool = False):
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
    def persistence_stability_index(self, ice_prstnc, A, persistence = 0.8):
        """
        Compute the Hemispheric Persistence Stability Index.

        This function quantifies the spatial stability of sea ice presence over time by comparing the
        area of persistent ice coverage to the total area where ice is ever present. It returns a single
        index value representing the ratio of consistently persistent ice to all intermittently
        present ice across the hemisphere.

        Parameters
        ----------
        ice_prstnc : xarray.DataArray
            Time-varying persistence of ice presence (0/1 or probability).
        A : xarray.DataArray
            Grid cell area (e.g., `TAREA`) on the same spatial grid.
        persistence : float, optional
            Minimum persistence threshold to define "persistent" ice (default is 0.8).

        Returns
        -------
        dict; containing a single scalar value:
            - 'persistence_stability_index': float; ratio of persistent ice area to total ice footprint.

        Notes
        -----
        - The persistent footprint is the summed area of cells with mean persistence >= threshold.
        - The total footprint includes all cells with nonzero presence over time.
        - Both areas are computed over the provided grid cell area `A`.
        - The result is computed using Dask for efficient parallel execution.
        - Requires that the input arrays are aligned and dimensionally compatible.
        """
        self.logger.info("Computing Hemispheric Stability Index: (persistent_footprint / total_footprint)")
        ice_prstnc = ice_prstnc.load()
        if 'time' in A.dims:
            A = A.isel(time=0).drop_vars('time')
        prstnc_mask     = ice_prstnc >= persistence      
        always_mask     = ice_prstnc > 0                 
        prstnc_A        = (A.where(prstnc_mask)).sum(dim=self.CICE_dict["spatial_dims"])
        ttl_A           = (A.where(always_mask)).sum(dim=self.CICE_dict["spatial_dims"])
        prstnc_A, ttl_A = dask.compute(prstnc_A, ttl_A)
        psi             = (prstnc_A / ttl_A) if ttl_A > 0 else np.nan
        return {"persistence_stability_index": psi}

    def persistence_ice_distance_mean_max(self, ice_prstnc, persistence_min=0.8):
        """
        Estimate how far sea ice (e.g., fast or pack ice) extends from the coast.

        This function calculates the mean and maximum distance of persistent sea ice presence
        from the coastline using a Euclidean distance transform. Rather than assuming constant
        grid spacing, it leverages the native grid cell area (`TAREA`) to estimate physical
        distances in kilometers. Grid cells are included if their ice persistence exceeds
        `persistence_min`.

        Parameters
        ----------
        ice_prstnc : xarray.DataArray
            Time-averaged presence of sea ice (e.g., fast ice mask mean over time).
        persistence_min : float, optional
            Minimum persistence threshold to include a grid cell in the analysis (default is 0.5).

        Returns
        -------
        dict
            Dictionary containing:
            - 'persistence_mean_dist' : float
                Mean distance (km) from coastline for cells with persistent ice.
            - 'persistence_max_dist' : float
                Maximum such distance (km).

        Notes
        -----
        - Requires the B-grid (`self.G_t`) to be loaded with `kmt_mod` and `area` fields.
        - Automatically reloads the B-grid with hemisphere slicing if needed.
        - Identifies coastal ocean grid cells using binary dilation.
        - Uses `TAREA` (in m^2) to derive an approximate local grid spacing.
        """
        from scipy.ndimage import distance_transform_edt, binary_dilation
        spat_dims = self.CICE_dict["spatial_dims"]
        self.logger.info("Computing persistence distance statistics mean & max")
        ice_prstnc = ice_prstnc.load()
        # Ensure bgrid is loaded and check for shape match
        if not getattr(self, "bgrid_loaded", False):
            self.load_bgrid(slice_hem=True)
        kmt = self.G_t['kmt_mod'].values
        if kmt.shape != tuple(ice_prstnc.sizes[dim] for dim in spat_dims):
            kmt = None
            self.logger.warning("    grid shape mismatch detected: reloading bgrid with slice_hem=True")
            self.load_bgrid(slice_hem=True)
            kmt = self.G_t['kmt_mod'].values
        land_mask   = (kmt == 0)
        sea_mask    = ~land_mask
        cst_mask    = sea_mask & binary_dilation(land_mask)
        dist_idx    = distance_transform_edt(~cst_mask)
        tarea       = self.G_t['area'].values
        grid_dx_km  = np.sqrt(tarea) / 1000.0
        cst_dist_km = dist_idx * grid_dx_km
        dist_da     = xr.DataArray(data   = cst_dist_km,
                                   dims   = spat_dims,
                                   coords = {dim: ice_prstnc[dim] for dim in spat_dims},
                                   name   = "distance_to_coast")
        self.logger.info(f"    applying persistence minimum filter {persistence_min} to ice persistence data array")
        prstnc_filt = (ice_prstnc >= persistence_min)
        self.logger.info("   comptuing distances...")
        ice_dists   = dist_da.where(prstnc_filt).compute()
        mean_dist   = ice_dists.mean().values
        max_dist    = ice_dists.max().values
        return {f'persistence_mean_distance': mean_dist,
                f'persistence_max_distance': max_dist}

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
                mod_grp  = mod.groupby('time.dayofyear')#.persist()
                mod_clim = mod_grp.mean('time')#.compute()
                obs_grp  = obs.groupby('time.dayofyear')#.persist()
                obs_clim = obs_grp.mean('time')#.compute()
                doy_mod  = mod_clim.dayofyear.values
                doy_obs  = obs_clim.dayofyear.values
                doy_xsct = np.intersect1d(doy_mod, doy_obs)
                mod_data = mod_clim.sel(dayofyear=doy_xsct).values
                obs_data = obs_clim.sel(dayofyear=doy_xsct).values
                return self._skill_stats(mod_data, obs_data)
        self.logger.warning("Time dimension not found or could not align. Returning empty stats.")
        return {"Bias"    : np.nan,
                "RMSE"    : np.nan,
                "MAE"     : np.nan,
                "Corr"    : np.nan,
                "SD_Model": np.nan,
                "SD_Obs"  : np.nan,}

