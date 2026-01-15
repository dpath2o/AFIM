from __future__ import annotations
import os
import xarray   as xr
import pandas   as pd
import numpy    as np
from pathlib    import Path
from typing     import Iterable, List, Optional, Tuple
__all__ = ["SeaIceObservations"]
class SeaIceObservations:
    """
    Observational sea-ice I/O and derived diagnostics for AFIM workflows.

    This class centralises reading, standardisation, and metric generation for several
    observational products used throughout AFIM/SeaIceToolbox analyses:

    - NSIDC (e.g., G02202 v4) sea-ice concentration:
        * local file discovery and loading over time windows (daily or monthly layouts),
        * masking of flagged concentration values,
        * computation of hemispheric sea-ice area (SIA) and extent (SIE) time series.

    - AF2020 (Fraser et al. 2020) landfast sea-ice:
        * reading FIA CSV time series and constructing DOY climatologies,
        * regridding original AF2020 raster products to the CICE T-grid with xESMF,
        * convenience regional fast-ice area time series from AF2020-like gridded datasets.

    - ESA-CCI / AWI CryoSat/Envisat/Sentinel sea-ice thickness (SIT):
        * robust path discovery across multiple institutional directory layouts,
        * indexing and de-duplication (one file per day or month),
        * L2P swath → daily gridded SIT via pyresample neighbour queries with optional
          Gaussian weighting and explicit uncertainty propagation,
        * L3 gridded monthly SIT assembly (no resampling), with optional SIC thresholding,
          QC flag masking, and hemispheric weighted means.

    Expected configuration and dependencies
    -------------------------------------
    This class is intended to be used inside the AFIM stack and expects several
    attributes to exist on `self`:

    Required attributes (typical)
    - logger : logging.Logger
        Used for progress and warnings.
    - dt0_str, dtN_str : str
        Default analysis window in "YYYY-MM-DD" format.
    - hemisphere_dict : dict
        Must include 'abbreviation' (e.g., "SH" or "NH").
    - NSIDC_dict : dict
        Paths and dataset conventions for NSIDC loading and metrics. Typical keys:
            * "D_original"
            * "G02202_v4_file_format"
            * "file_versions" (dict ver -> "YYYY-MM-DD" start date, or None)
            * "SIC_name" (e.g., "cdr_seaice_conc")
            * "cdr_seaice_conc_flags" (list of flag integer codes)
            * "P_cell_area" (cell area file path)
            * "projection_string" (Proj4 string)
    - AF_FI_dict : dict
        Paths and naming for AF2020. Typical keys:
            * "P_AF2020_csv"
            * "D_AF2020_db_org"
            * "AF_reG_weights"
            * "variable_name", "lat_coord_name", "lon_coord_name", "time_coord_name"
            * "threshold_value"
            * "P_reG_reg_weights" (optional)
    - CICE_dict : dict
        CICE grid conventions and coordinate names (for mapping masks/regridding):
            * 'spatial_dims', 'three_dims'
            * 'lon_coord_name', 'lat_coord_name'
            * 'tcoord_names' (lon, lat)
            * 'time_dim'
            * 'P_reG_reg_weights' (optional)
    - Sea_Ice_Obs_dict : dict
        Root directories for institutions:
            * 'ESA-CCI', 'AWI'
    - icon_thresh : float
        Sea-ice concentration threshold for extent masking (fraction).
    - SIC_scale : float
        Unit scaling between dataset SIC and area weighting.
    - normalise_longitudes(lon, to=...) : callable
        Longitude wrapping helper used across multiple routines.
    - define_cice_grid(...) : callable
        Must be available for AF2020 conservative regridding to the CICE T-grid.
    - load_cice_grid(...) / define_cice_grid(...) conventions
        Must be consistent with your xESMF usage.

    Notes on time and units
    -----------------------
    - All internal timestamp comparisons for file filtering are performed in UTC.
    - NSIDC SIA/SIE are computed using an explicit cell-area product and your
      configured concentration scaling.
    - SIT uncertainty propagation uses standard-error style aggregation:
        SE_cell  = sqrt(sum(w^2 * unc^2)) / sum(w)
        SE_hemi  = sqrt(sum(W^2 * SE_cell^2)) / sum(W)
      where W is cos(latitude) and w is neighbour weight.

    Implementation note (RTD/autodoc)
    --------------------------------
    Keep docstrings on public methods comprehensive and consistent with your JSON keys.
    Methods with a leading underscore are treated as internal helpers.
    """

    def __init__(self, **kwargs):
        """
        Initialise the observations helper.

        Parameters
        ----------
        **kwargs
            Optional configuration injected into `self` by the surrounding workflow.

        Notes
        -----
        - The shown implementation is a no-op; in production you typically assign kwargs
          onto `self` or inherit from a base class that handles configuration.
        """
        return
    
    ####################################################################################################
    ##                                              NSIDC         
    ####################################################################################################
    def load_local_NSIDC(self, dt0_str=None, dtN_str=None, local_directory=None, monthly_files=False):
        """
        Load NSIDC sea-ice concentration NetCDF files from a local directory tree.

        The loader constructs expected filenames for each day (or month) between
        `dt0_str` and `dtN_str`, selects the appropriate file version based on your
        configured version start dates, and uses `xarray.open_mfdataset` to open all
        available files in one dataset.

        Parameters
        ----------
        dt0_str, dtN_str : str, optional
            Start/end dates in "YYYY-MM-DD" format. If omitted, defaults to
            `self.dt0_str` and `self.dtN_str`.
        local_directory : str or pathlib.Path, optional
            Root directory containing NSIDC files. If omitted, defaults to
            `Path(self.NSIDC_dict["D_original"], self.hemisphere, freq_str)`.
        monthly_files : bool, default False
            If True, treat the dataset as monthly and search one file per month.
            If False, treat the dataset as daily.

        Returns
        -------
        xr.Dataset
            Multi-file dataset combined "by_coords".

        Raises
        ------
        FileNotFoundError
            If no matching files are found within the requested window.

        Notes
        -----
        - A small preprocess hook promotes `time` as a coordinate when datasets contain
          a `tdim` dimension.
        - Missing files are logged as warnings and skipped.
        """
        from datetime import datetime
        def promote_time(ds):
            ds = ds[["cdr_seaice_conc"]] if "cdr_seaice_conc" in ds else ds
            if "time" in ds.variables and "tdim" in ds.dims:
                ds = ds.swap_dims({"tdim": "time"})
                ds = ds.set_coords("time")
            return ds
        f_fmt    = self.NSIDC_dict["G02202_v4_file_format"]
        dt0_str  = dt0_str if dt0_str is not None else self.dt0_str
        dtN_str  = dtN_str if dtN_str is not None else self.dtN_str
        freq_str = 'monthly'             if monthly_files  else 'daily'
        D_local  = Path(local_directory) if local_directory else Path(self.NSIDC_dict["D_original"], self.hemisphere, freq_str)
        dt0      = datetime.strptime(dt0_str, '%Y-%m-%d')
        dtN      = datetime.strptime(dtN_str, '%Y-%m-%d')
        f_vers_parsed = [ (ver, datetime.strptime(date_str, "%Y-%m-%d") if date_str else datetime.min)
                          for ver, date_str in self.NSIDC_dict["file_versions"].items() ]
        f_vers_parsed.sort(key=lambda x: x[1])  # Sort by date
        date_range = pd.date_range(start=dt0, end=dtN, freq='MS' if monthly_files else 'D')
        P_NSIDCs   = []
        for dt_ in date_range:
            f_version = next( (ver for ver, d in reversed(f_vers_parsed) if dt_ >= d), f_vers_parsed[0][0] )
            F_        = f_fmt.format(freq = freq_str,
                                     hem  = self.hemisphere_dict['abbreviation'].lower(),
                                     date = dt_.strftime('%Y%m%d'),
                                     ver  = f_version)
            P_        = Path(D_local,F_)
            if P_.exists():
                P_NSIDCs.append(P_)
            else:
                self.logger.warning(f"Missing file: {P_}")
        if not P_NSIDCs:
            raise FileNotFoundError(f"No NSIDC files found in {D_local} between {dt0_str} and {dtN_str}")
        self.logger.info(f"Using xarray to load {len(P_NSIDCs)} files in {D_local}. Last file: {P_NSIDCs[-1].name}")
        return xr.open_mfdataset(P_NSIDCs,
                                 combine    = "by_coords",
                                 parallel   = True,
                                 preprocess = promote_time)

    def NSIDC_coordinate_transformation(self, ds):
        """
        Add geographic longitude/latitude coordinates to an NSIDC grid.

        NSIDC concentration products are commonly stored in a Polar Stereographic
        projection with Cartesian coordinates (e.g., `xgrid`, `ygrid`, in meters).
        This routine converts the projected coordinates to WGS84 longitude/latitude
        and attaches them to the dataset.

        Parameters
        ----------
        ds : xr.Dataset
            Input dataset containing projected coordinate variables. Expected variables:
            - xgrid : (x) or (y, x) coordinate(s) in meters
            - ygrid : (y) or (y, x) coordinate(s) in meters

        Returns
        -------
        xr.Dataset
            Dataset with added variables:
            - lon : (y, x) longitude in degrees east
            - lat : (y, x) latitude in degrees north

        Raises
        ------
        KeyError
            If xgrid/ygrid variables are missing.
        AssertionError
            If xgrid and ygrid shapes are incompatible for transformation.
        ImportError
            If pyproj is not installed.

        Notes
        -----
        - Uses `self.NSIDC_dict["projection_string"]` (Proj4) for the source CRS.
        - Uses EPSG specified in `self.sea_ice_dic["projection_wgs84"]` (typically 4326)
          for the target CRS.
        - This method assumes x/y describe the native NSIDC grid in meters.
        """
        from pyproj import CRS, Transformer
        x = ds['xgrid'].values
        y = ds['ygrid'].values
        assert x.shape == y.shape, "xgrid and ygrid must have the same shape"
        crs_nsidc = CRS.from_proj4(self.NSIDC_dict["projection_string"])
        crs_wgs84 = CRS.from_epsg(self.sea_ice_dic["projection_wgs84"])
        trans     = Transformer.from_crs(crs_proj, crs_wgs84, always_xy=True)
        lon, lat  = transformer.transform(x, y)
        ds['lat'] = (('y', 'x'), lat)
        ds['lon'] = (('y', 'x'), lon)
        return ds

    def compute_NSIDC_metrics(self,
                              dt0_str         = None,
                              dtN_str         = None,
                              local_load_dir  = None,
                              monthly_files   = False,
                              P_zarr          = None,
                              overwrite       = False):
        """
        Compute and persist NSIDC hemispheric sea-ice area (SIA) and extent (SIE).

        This routine loads NSIDC concentration for the specified window, masks out
        flagged concentration values, applies a concentration threshold to define
        "ice-covered" cells, and then computes hemispheric aggregates using the NSIDC
        cell-area product.

        Parameters
        ----------
        dt0_str, dtN_str : str, optional
            Start/end dates in "YYYY-MM-DD" format. Defaults to `self.dt0_str` and
            `self.dtN_str`.
        local_load_dir : str or pathlib.Path, optional
            Directory containing the NSIDC source NetCDFs. If omitted, uses your
            configured `self.NSIDC_dict["D_original"]` hemisphere sub-tree.
        monthly_files : bool, default False
            If True, compute metrics using monthly inputs; otherwise daily inputs.
        P_zarr : str or pathlib.Path, optional
            Output Zarr path. Defaults to `<local_load_dir>/zarr`.
        overwrite : bool, default False
            If True, recompute even if the Zarr store already exists.

        Returns
        -------
        xr.Dataset
            Dataset with time series:
            - SIA(time) : sea-ice area (units determined by your scaling)
            - SIE(time) : sea-ice extent (same area units)
            plus a time coordinate.

        Raises
        ------
        FileNotFoundError
            If no NSIDC source files are found for the requested window.
        KeyError
            If required variables/paths are missing from your NSIDC configuration.

        Notes
        -----
        - Flags are defined in `self.NSIDC_dict["cdr_seaice_conc_flags"]` and removed
          by setting affected concentrations to NaN before aggregation.
        - Area is loaded from `self.NSIDC_dict["P_cell_area"]`.
        - Extent mask is `aice > self.icon_thresh`.
        """
        flags    = self.NSIDC_dict["cdr_seaice_conc_flags"]
        SIC_name = self.NSIDC_dict["SIC_name"]
        dt0_str  = dt0_str if dt0_str is not None else self.dt0_str
        dtN_str  = dtN_str if dtN_str is not None else self.dtN_str
        if monthly_files:
            freq_str = 'monthly'
        else:
            freq_str = 'daily'
        D_local = local_load_dir if local_load_dir is not None else Path(self.NSIDC_dict["D_original"], self.hemisphere, freq_str)
        P_zarr  = P_zarr         if P_zarr         is not None else Path(D_local, "zarr" )
        if not os.path.exists(P_zarr) or overwrite:
            NSIDC = self.load_local_NSIDC(dt0_str         = dt0_str,
                                          dtN_str         = dtN_str,
                                          local_directory = D_local)
            aice  = NSIDC[SIC_name]
            for flag in flags:
                aice = xr.where(aice==flag/100, np.nan, aice)
            area     = xr.open_dataset(self.NSIDC_dict["P_cell_area"]).cell_area / self.SIC_scale
            mask     = aice > self.icon_thresh
            SIA      = (aice * area).where(mask.notnull()).sum(dim=['y', 'x'], skipna=True)
            SIE      = (mask * area).sum(dim=['y', 'x'], skipna=True)
            NSIDC_ts = xr.Dataset({'SIA' : (('time'),SIA.values),
                                   'SIE' : (('time'),SIE.values)},
                                  coords = {'time' : (('time'), SIA.time.values)})
            self.logger.info(f"**writing** NSIDC time series zarr file {P_zarr}")
            NSIDC_ts.to_zarr(P_zarr)
            return NSIDC_ts
        else:
            self.logger.info(f"**loading** previously created NSIDC time series zarr file {P_zarr}")
            return xr.open_zarr(P_zarr, consolidated=True)

    ####################################################################################################
    ##                                              AF2020         
    ####################################################################################################
    def load_AF2020_FIA_summary(self, 
                                repeat_climatology: bool = True,
                                start             : str  = "1994-01-01",
                                end               : str  = "1999-12-31") -> xr.Dataset:
        """
        Load AF2020 landfast sea-ice area (FIA) time series and build a DOY climatology.

        This method reads the AF2020 FIA CSV, standardises sector names, produces a
        circumpolar daily time series, and computes a day-of-year (DOY) climatology.
        Optionally, the climatology is repeated over an arbitrary target time range
        for direct comparison against model series.

        Parameters
        ----------
        repeat_climatology : bool, default True
            If True, returns a repeated climatology series (`FIA_clim_repeat`) over the
            specified [start, end] window. If False, returns only the raw FIA series
            and DOY climatology.
        start, end : str, default ("1994-01-01", "1999-12-31")
            Date window used when repeating the climatology (daily frequency).

        Returns
        -------
        xr.Dataset
            Variables:
            - FIA_obs(time, region) : observed FIA time series (typically "circumpolar")
            - FIA_clim(doy, region) : DOY climatology
            - FIA_clim_repeat(time, region) : climatology mapped to [start, end] (optional)

        Raises
        ------
        FileNotFoundError
            If the CSV path does not exist.
        ValueError
            If expected columns required for climatology construction are missing.

        Notes
        -----
        - Units are converted from km² to 10³ km² by dividing by 1000.
        - Leap-day handling: this implementation drops DOY 366 when repeating unless
          explicitly handled; ensure consistency with your model calendar.
        """
        from scipy.interpolate import interp1d
        csv_path = self.AF_FI_dict['P_AF2020_csv']
        df = pd.read_csv(csv_path)
        regions = ["circumpolar", "indian_ocean", "western_pacific", "ross_sea", "amundsen_sea", "weddell_sea"]
        df = df.rename(columns={'Circumpolar': 'circumpolar',
                                'IOsector'   : 'indian_ocean',
                                'WPOsector'  : 'western_pacific',
                                'RSsector'   : 'ross_sea',
                                'BASsector'  : 'amundsen_sea',
                                'WSsector'  : 'weddell_sea'})
        df = df.dropna(subset=["Year", "Month_start", "Day_of_month_start"])
        df["date"] = pd.to_datetime({"year": df["Year"].astype(int),
                                    "month": df["Month_start"].astype(int),
                                    "day": df["Day_of_month_start"].astype(int)}, errors="coerce")
        df = df.dropna(subset=["date"])
        df["doy"] = df["date"].dt.dayofyear
        df = df.set_index("date")
        df = df[regions + ["doy"]]
        fia_obs = df[regions] / 1000.0  # km² → 1000 km²
        fia_obs_xr = fia_obs.to_xarray().to_array("region").transpose("date", "region")
        if 'doy' not in df.columns or 'circumpolar' not in df.columns:
            raise ValueError("Expected 'doy' and 'circumpolar' columns for climatology generation.")
        doy_vals = np.arange(1, 366)
        x = df.groupby("doy")["circumpolar"].mean().dropna().index.values
        y = df.groupby("doy")["circumpolar"].mean().dropna().values / 1000.0
        interp_func = interp1d(x, y, kind='linear', fill_value="extrapolate")
        clim_interp = pd.DataFrame({"circumpolar": interp_func(doy_vals)}, index=doy_vals)
        fia_clim_doy = xr.DataArray(clim_interp["circumpolar"].values[:, None],  # <-- Add [:, None] to make it 2D
                                    coords={"doy": doy_vals, "region": ["circumpolar"]},
                                    dims=["doy", "region"])

        if repeat_climatology:
            model_dates = pd.date_range(start, end, freq="D")
            doy_map = model_dates.dayofyear
            valid_mask = doy_map < 366
            model_dates = model_dates[valid_mask]
            doy_map = doy_map[valid_mask]
            repeated = fia_clim_doy.sel(doy=xr.DataArray(doy_map, dims="time")).values
            fia_clim_re = xr.DataArray(repeated,
                                       coords = {"time": model_dates, "region": ["circumpolar"]},
                                       dims   = ["time", "region"])
        else:
            fia_clim_re = None
        ds = xr.Dataset({"FIA_obs": fia_obs_xr.sel(region="circumpolar"),
                        "FIA_clim": fia_clim_doy})
        if fia_clim_re is not None:
            ds["FIA_clim_repeat"] = fia_clim_re
        self.logger.info(f"AF2020 FIA data loaded from {csv_path}")
        return ds

    def AF2020_interpolate_FIA_clim(self, csv_df, full_doy=np.arange(1, 366)):
        """
        Interpolate an AF2020 DOY climatology to a complete 1..365 DOY axis.

        Parameters
        ----------
        csv_df : pandas.DataFrame
            DataFrame containing at least:
            - 'DOY_start' : day-of-year index values
            - 'Circumpolar' : FIA values (km²)
        full_doy : array-like, default np.arange(1, 366)
            Target DOY coordinate to interpolate onto.

        Returns
        -------
        xr.DataArray
            Interpolated climatology with coordinate:
            - doy : day-of-year
            Values are in 10³ km².
        """
        from scipy.interpolate import interp1d
        df = csv_df.copy()
        if 'Circumpolar' not in df.columns or 'DOY_start' not in df.columns:
            raise ValueError("Expected columns 'DOY_start' and 'circumpolar' in CSV.")
        x = df['DOY_start'].values
        y = df['Circumpolar'].values/1e3
        if len(x) != len(y):
            raise ValueError("x and y arrays must be equal in length.")
        interp_func = interp1d(x, y, kind='linear', fill_value="extrapolate")
        return xr.DataArray(interp_func(full_doy), coords=[("doy", full_doy)])

    def AF2020_clim_to_model_time(self, model: xr.DataArray, clim: xr.DataArray) -> xr.DataArray:
        """
        Expand a DOY-based climatology to match a model time axis.

        Parameters
        ----------
        model : xr.DataArray
            Model time series containing a `time` coordinate (datetime-like).
        clim : xr.DataArray
            Climatology indexed by `doy` (and optionally `region`).

        Returns
        -------
        xr.DataArray
            Climatology values mapped onto `model.time` using each date's day-of-year.

        Notes
        -----
        - Leap days: if the observational climatology does not define DOY=366, callers
          should drop Feb 29 from the model calendar or define a mapping strategy
          (e.g., nearest/pad) consistent with the analysis.
        """
        # Compute day-of-year for each model date (drop Feb 29 if not in obs)
        doy = xr.DataArray(model['time'].dt.dayofyear, coords={"time": model.time}, dims="time")
        # Match each model day to a climatological value
        matched = clim.sel(doy=doy, method="nearest")  # or method="pad"
        matched.coords["time"] = model.time  # Align with model time
        return matched

    def build_AF2020_grid_corners(self, lat_c, lon_c):
        """
        Construct (ny+1, nx+1) corner arrays from (ny, nx) center arrays.

        This helper builds approximate corner coordinates by padding center fields
        and averaging surrounding values. Corner arrays are commonly required for
        conservative remapping in xESMF.

        Parameters
        ----------
        lat_c, lon_c : np.ndarray
            2D arrays (ny, nx) of center latitude/longitude in degrees.

        Returns
        -------
        (lat_b, lon_b) : tuple[np.ndarray, np.ndarray]
            2D arrays (ny+1, nx+1) of corner latitude/longitude.

        Notes
        -----
        - This approach assumes a locally smooth grid. For strongly curvilinear grids,
          corner construction should ideally use native corner metadata when available.
        - Output longitudes are wrapped to [0, 360).
        """

        ny, nx = lat_c.shape
        lat_padded = np.pad(lat_c, ((1, 1), (1, 1)), mode='edge')
        lon_padded = np.pad(lon_c, ((1, 1), (1, 1)), mode='edge')
        lat_b = 0.25 * (lat_padded[0:-1, 0:-1] + lat_padded[0:-1, 1:] +
                        lat_padded[1:, 0:-1] + lat_padded[1:, 1:])
        lon_b = 0.25 * (lon_padded[0:-1, 0:-1] + lon_padded[0:-1, 1:] +
                        lon_padded[1:, 0:-1] + lon_padded[1:, 1:])
        lon_b = ((lon_b + 360) % 360)
        return lat_b, lon_b

    def _load_AF2020_FI_org_and_save_zarr(self, P_zarr):
        """
        Load original AF2020 FastIce NetCDF rasters, compute persistence, regrid to CICE T-grid, and write Zarr.

        Workflow
        --------
        1) Open all AF2020 "FastIce_70_*.nc" rasters (mosaic-window products).
        2) Build an xESMF conservative regridder from the AF2020 grid to the CICE T-grid,
           reusing a persistent weights file when available.
        3) Convert the AF2020 categorical time-series variable to a binary fast-ice mask
           using a threshold (>=4), then compute a time-mean persistence fraction.
        4) Regrid persistence to the CICE grid and save a compact dataset to Zarr.

        Parameters
        ----------
        P_zarr : pathlib.Path
            Destination Zarr store path.

        Returns
        -------
        xr.Dataset
            Dataset containing:
            - FI(nj, ni) : fast-ice persistence fraction on CICE T-grid
            - FI_t_alt(t_FI_obs) : auxiliary date strings from the original product
            plus lon/lat coordinates on the CICE T-grid.

        Notes
        -----
        - Weight file path is taken from `self.AF_FI_dict["AF_reG_weights"]`.
        - Target grid is built via `self.define_cice_grid(grid_type="t", ...)`.
        - Conservative regridding requires valid grid corners on both source and target.
        """
        import xesmf as xe
        D_obs  = Path(self.AF_FI_dict['D_AF2020_db_org'])
        P_orgs = sorted(D_obs.glob("FastIce_70_*.nc"))
        self.logger.info(f"Loading original AF2020 NetCDF files: {P_orgs}")
        FI_obs        = xr.open_mfdataset(P_orgs, engine='netcdf4', combine='by_coords')
        FI_obs        = FI_obs.chunk({'time': 32})
        lat_c         = FI_obs.latitude
        lon_c         = FI_obs.longitude
        lat_b, lon_b  = self.build_AF2020_grid_corners(lat_c.values, lon_c.values)
        G_obs         = xr.Dataset({"lat"   : (("Y", "X"), lat_c),
                                    "lon"   : (("Y", "X"), lon_c),
                                    "lat_b" : (("Y_b", "X_b"), lat_b),
                                    "lon_b" : (("Y_b", "X_b"), lon_b)})
        G_t           = self.define_cice_grid(grid_type='t', mask=False, build_grid_corners=True)
        P_weights     = self.AF_FI_dict["AF_reG_weights"]
        reuse_weights = os.path.exists(P_weights)
        self.logger.info(f"Regridding observational data using weights: {P_weights} (reuse={reuse_weights})")
        reG_obs = xe.Regridder(G_obs, G_t,
                               method            = "conservative",
                               periodic          = False,
                               ignore_degenerate = True,
                               reuse_weights     = reuse_weights,
                               filename          = P_weights)
        fast_ice_mask = xr.where(FI_obs['Fast_Ice_Time_series'] >= 4, 1.0, 0.0)
        FI_frac       = fast_ice_mask.sum(dim='time') / fast_ice_mask['time'].size
        FI_frac       = FI_frac.where(FI_frac > 0)  # set 0s to NaN (optional: mask ocean)
        FI_frac_reG   = reG_obs(FI_frac)
        ds_out = xr.Dataset(data_vars = {"FI"       : FI_frac_reG,
                                         "FI_t_alt" : ("t_FI_obs", FI_obs.date_alt, {"long_name"   : FI_obs.date_alt.attrs.get("long_name", ""),
                                                                                     "description" : FI_obs.date_alt.attrs.get("description", "")})},
                            coords = {"t_FI_obs" : ("t_FI_obs", FI_obs.time, {"description": "Start date of 15- or 20-day image mosaic window.",
                                                                              "units": "days since 2000-1-1 0:0:0"}),
                                      "lon"      : (("nj", "ni"), G_t['lon'].values, {"units": "degrees_east"}),
                                      "lat"      : (("nj", "ni"), G_t['lat'].values, {"units": "degrees_north"})})
        self.logger.info(f"Saving regridded dataset to Zarr: {P_zarr}")
        ds_out.to_zarr(P_zarr, mode="w", consolidated=True)
        return ds_out

    def AF2020_regional_fast_ice_area(self, ds: xr.Dataset,
                                        region    = (52.5, 97.5, -70.5, -64),
                                        lon_name  = "longitude",
                                        lat_name  = "latitude",
                                        area_name = "area",
                                        fast_name = "Fast_Ice_Time_series",
                                        threshold = 4,
                                        lon_wrap  = "0-360"):
        """
        Compute a regional fast-ice area time series from a gridded observation product.

        This routine is designed for AF2020-like datasets that provide:
        - curvilinear or projected 2-D lon/lat fields, and
        - a per-cell area field (m²), and
        - a fast-ice indicator variable over time.

        The region is defined by a geographic bounding box:
            region = (lon_min, lon_max, lat_min, lat_max)

        Longitude seam-crossing is supported (e.g., lon_min=350, lon_max=20 in 0–360
        convention). The computation is performed as:

            FIA(t) = Σ_{cells in region} [ fast(t,cell) * area(cell) ]

        where `fast(t,cell)` is a boolean mask defined by `ds[fast_name] >= threshold`.

        Parameters
        ----------
        ds : xarray.Dataset
            Dataset containing at minimum:
            - ds[lon_name]  : 2-D longitude field with a time dimension (time, Y, X) or
                                broadcastable; only time=0 is used for the static grid.
            - ds[lat_name]  : 2-D latitude field with a time dimension (time, Y, X) or
                                broadcastable; only time=0 is used for the static grid.
            - ds[area_name] : 2-D grid-cell area in m² with a time dimension (time, Y, X)
                                or broadcastable; only time=0 is used.
            - ds[fast_name] : fast-ice indicator variable with a time dimension.
        region : tuple[float, float, float, float], default (52.5, 97.5, -70.5, -64)
            Geographic region bounds as (lon_min, lon_max, lat_min, lat_max), in degrees.
        lon_name : str, default "longitude"
            Name of the longitude variable/coordinate in `ds`.
        lat_name : str, default "latitude"
            Name of the latitude variable/coordinate in `ds`.
        area_name : str, default "area"
            Name of the per-cell area variable/coordinate in `ds` (assumed m²).
        fast_name : str, default "Fast_Ice_Time_series"
            Name of the fast-ice indicator variable in `ds`.
        threshold : float, default 4
            Threshold applied to `ds[fast_name]` to produce a boolean fast-ice mask.
            Cells with `ds[fast_name] >= threshold` are treated as fast ice.
        lon_wrap : {"0-360", "-180-180"}, default "0-360"
            Longitude convention used for the region definition and for wrapping `lon2d`.

        Returns
        -------
        xarray.DataArray
            1-D time series of regional fast-ice area in **1e3 km²** (i.e., thousands of km²),
            named `"FIA_1e3-km2"`. The time coordinate is inherited from `ds[fast_name]`.

        Notes
        -----
        - Lon/lat/area are read from `time=0` only to avoid unnecessary time broadcasting.
        - The regional mask is optionally cropped to the minimal bounding index range to
        reduce computation.
        - If `region` crosses the longitude seam (lon_min > lon_max), the longitude mask
        is built using an OR condition to include both sides.

        Examples
        --------
        >>> fia = obs.AF2020_regional_fast_ice_area(ds, region=(350, 20, -75, -65), lon_wrap="0-360")
        >>> fia.plot()
        """
        lon_min, lon_max, lat_min, lat_max = region
        # --- 1) Use static lon/lat & area from first time slice (avoid time broadcast)
        lon2d = ds[lon_name].isel(time=0)
        lat2d = ds[lat_name].isel(time=0)
        area2d = ds[area_name].isel(time=0)  # assumed m^2
        # --- 2) Normalize longitudes to match region convention
        def wrap(lon):
            if lon_wrap == "0-360":
                return ((lon % 360) + 360) % 360
            else:
                return ((lon + 180.0) % 360.0) - 180.0
        lon2d = wrap(lon2d)
        if lon_wrap == "0-360":
            lon_min, lon_max = lon_min % 360, lon_max % 360
        else:
            lon_min = ((lon_min + 180) % 360) - 180
            lon_max = ((lon_max + 180) % 360) - 180
        # --- 3) Build region mask (handle seam-crossing gracefully)
        if lon_min <= lon_max:
            mask_lon = (lon2d >= lon_min) & (lon2d <= lon_max)
        else:
            mask_lon = (lon2d >= lon_min) | (lon2d <= lon_max)  # across seam
        mask_lat = (lat2d >= lat_min) & (lat2d <= lat_max)
        region_mask = mask_lon & mask_lat                      # (Y, X)
        # --- 4) Convert area to 1000-km^2
        area_km2 = area2d / 1e9
        # OPTIONAL: crop to bounding indices to reduce work
        jdim, idim = lat2d.dims  # expect ("Y","X")
        j_any = region_mask.any(dim=idim)
        i_any = region_mask.any(dim=jdim)
        if bool(j_any.any()) and bool(i_any.any()):
            j_idx = np.where(j_any.values)[0]
            i_idx = np.where(i_any.values)[0]
            j0, j1 = int(j_idx.min()), int(j_idx.max()) + 1
            i0, i1 = int(i_idx.min()), int(i_idx.max()) + 1
            region_mask = region_mask.isel({jdim: slice(j0, j1), idim: slice(i0, i1)})
            area_km2    = area_km2.isel({jdim: slice(j0, j1), idim: slice(i0, i1)})
        # --- 5) Fast-ice presence per time (binary), then area sum over region
        fast = (ds[fast_name] >= threshold)
        # align cropped mask/area to fast dims (time, Y, X)
        region_mask3 = region_mask.broadcast_like(fast.isel(time=0))
        area_km2_3   = area_km2.broadcast_like(fast.isel(time=0))
        # area at each time = sum( fast * region_mask * area )
        fia_ts = (fast.where(region_mask3) * area_km2_3).sum(dim=(jdim, idim))
        fia_ts.name = "FIA_1e3-km2"
        fia_ts.attrs["description"] = f"Fast-ice area (1e3-km^2) in region {region}, threshold>={threshold}"
        return fia_ts

    ####################################################################################################
    ##                                          CCI SIT         
    ####################################################################################################
    def build_sea_ice_satellite_paths(self,
                                      dt0_str      : str,
                                      dtN_str      : str,
                                      hemisphere   : Optional[str] = None,
                                      institutions : Optional[Iterable[str]] = ("ESA","AWI"),
                                      sensors      : Optional[Iterable[str]] = None,          # ["cryosat2","envisat","sentinel3a","sentinel3b"]
                                      levels       : Optional[Iterable[str]] = None,          # ESA: ["L2P","L3C"]; AWI: ["l2p_release","l3cp_release"]
                                      versions     : Optional[Iterable[str]] = None,          # ESA only; and only "v2.0" or "v3.0"
                                      root_esa     : Optional[str] = None, 
                                      root_awi     : Optional[str] = None) -> List[Path]:
        """
        Build a de-duplicated, time-filtered list of NetCDF files for satellite sea-ice
        products across ESA and/or AWI directory layouts.

        Time filtering:
            - Returns files whose date/month lies within the closed interval
            [dt0_str, dtN_str], interpreted in UTC.
            - Daily products (e.g., ESA L2P) are filtered by date parsed from the filename
            (YYYYMMDD). Monthly products (e.g., ESA L3C or AWI monthly layout) are
            filtered by overlap between the (year, month) and [dt0, dtN].

        Supported directory patterns (walked automatically):
            ESA “flat” layout:
                <root_esa>/<L2P|L3C>/<sensor>/<hemisphere>/*.nc
                (hemisphere directory may be 'NH'/'SH' or lowercase)
                - L2P: filenames contain YYYYMMDD
                - L3C: filenames often contain YYYYMM

            ESA versioned layout:
                <root_esa>/<L2P|L3C>/<sensor>/<version>/<HEMISPHERE>/<YYYY>/<MM>/*.nc
                (HEMISPHERE can be 'NH'/'SH' or lowercase; <MM> subdirs may be absent)

            AWI release layout:
                <root_awi>/<l2p_release|l3cp_release>/<hemisphere>/<sensor>/<YYYY>[/<MM>]/*.nc

        Parameters
        ----------
        dt0_str, dtN_str
            Start and end of the desired time window (inclusive). Any pandas-parsable
            string is accepted (e.g., "2002-01-01", "2012-12-31"). Interpreted as UTC.
        hemisphere
            Hemisphere selector. If None, uses `self.hemisphere_dict['abbreviation']`.
            Case is normalized; both 'sh'/'nh' and 'SH'/'NH' are supported in dir scans.
        institutions
            Iterable of institutions to include (case-insensitive). Valid values include
            "ESA" and/or "AWI". If None or empty, both are searched.
        sensors
            Optional iterable of sensor names to filter (case-insensitive, matched to
            directory names under the institution layout).
        levels
            Optional iterable of level strings. ESA: "L2P", "L3C".
            AWI: "l2p_release", "l3cp_release". If "l2p"/"l3c" are provided for AWI
            they are mapped to "l2p_release"/"l3cp_release".
        versions
            ESA only: one or more version directory names (e.g., "v2.0", "v3.0").
            Case-insensitive. Ignored for AWI.
        root_esa, root_awi
            Root directories for ESA and AWI datasets.

        Returns
        -------
        List[Path]
            Sorted unique list of file paths within the requested time window that
            satisfy the filters.

        Notes
        -----
        - Filename parsing assumes dates follow YYYYMMDD (daily) or YYYYMM (monthly)
        within 2000–2099. Files without such tokens are kept only when they live in
        a year/month directory that overlaps [dt0, dtN].
        - If `institutions` is None or empty, both ESA and AWI are scanned.
        - All timestamps are treated as UTC to avoid tz-comparison issues.

        Examples
        --------
        >>> paths = self.build_sea_ice_satellite_paths(
        ...     "2002-01-01", "2002-12-31",
        ...     hemisphere="sh",
        ...     institutions=["ESA", "AWI"],
        ...     sensors=["cryosat2", "envisat"],
        ...     levels=["L2P","L3C","l2p_release"]
        ... )
        >>> len(paths), paths[:3]
        """
        dt0_str = dt0_str if dt0_str is not None else self.dt0_str
        dtN_str = dtN_str if dtN_str is not None else self.dtN_str
        t0 = pd.Timestamp(dt0_str, tz="UTC")
        t1 = pd.Timestamp(dtN_str, tz="UTC")
        if t1 < t0:
            raise ValueError("dtN_str date must be >= dt0_str date")
        hem = hemisphere if hemisphere is not None else self.hemisphere_dict["abbreviation"]
        hem_upper = str(hem).upper()
        hem_lower = str(hem).lower()
        hem_dirs_esa = [hem_upper, hem_lower]  # try both cases for ESA
        root_esa = root_esa if root_esa is not None else self.Sea_Ice_Obs_dict['ESA-CCI']
        root_awi = root_awi if root_awi is not None else self.Sea_Ice_Obs_dict['AWI']
        root_esa = Path(root_esa)
        root_awi = Path(root_awi)
        # filters (case-normalized)
        insts = set(s.upper() for s in (self._norm_list(institutions) or []))
        sensor_set = set(s.lower() for s in (self._norm_list(sensors) or [])) if sensors else None
        version_set = set(v.lower() for v in (self._norm_list(versions) or [])) if versions else None
        level_set = set(l.lower() for l in (self._norm_list(levels) or [])) if levels else None
        out: List[Path] = []
        # ---------------- ESA ----------------
        if not insts or "ESA" in insts:
            esa_levels = ["L2P", "L3C"] if not level_set else ([L.upper() for L in level_set if L.upper() in ("L2P", "L3C")] or ["L2P", "L3C"])
            for L in esa_levels:
                Ldir = Path(root_esa,L)
                if not Ldir.exists():
                    continue
                # <root>/<L>/<sensor>/...
                for sdir in sorted(p for p in Ldir.iterdir() if p.is_dir()):
                    sensor = sdir.name.lower()
                    if sensor_set and sensor not in sensor_set:
                        continue
                    # A) flat: .../L2P|L3C/<sensor>/<hem>/*.nc   (hem 'nh'/'sh' any case)
                    for hem_dir in hem_dirs_esa:
                        flat_dir = Path(sdir,hem_dir)
                        if flat_dir.exists():
                            for fp in flat_dir.glob("*.nc"):
                                if L == "L2P":
                                    ts = self._parse_yyyymmdd(fp.name)
                                    if ts is None or ts < t0 or ts > t1:
                                        continue
                                else:  # L3C
                                    ym = self._parse_yyyymm(fp.name)
                                    if ym and not self._month_overlap(ym[0], ym[1], t0, t1):
                                        continue
                                out.append(fp)
                    # B) versioned: .../L2P|L3C/<sensor>/<ver>/<HEM>/YYYY/MM/*.nc
                    for vdir in sorted(p for p in sdir.iterdir() if p.is_dir()):
                        vname = vdir.name.lower()
                        # skip if this 'vdir' is actually the hemisphere dir from flat layout
                        if vname in {"nh", "sh"}:
                            continue
                        if version_set and vname not in version_set:
                            continue
                        for hem_dir in hem_dirs_esa:
                            hdir = Path(vdir,hem_dir)
                            if not hdir.exists():
                                continue
                            for ydir in sorted(hdir.glob("[12][0-9][0-9][0-9]")):
                                y = int(ydir.name)
                                mdirs = list(sorted(ydir.glob("[01][0-9]")))
                                if mdirs:
                                    for mdir in mdirs:
                                        m = int(mdir.name)
                                        if not self._month_overlap(y, m, t0, t1):
                                            continue
                                        for fp in sorted(mdir.glob("*.nc")):
                                            if L == "L2P":
                                                ts = self._parse_yyyymmdd(fp.name)
                                                if ts is None or ts < t0 or ts > t1:
                                                    continue
                                            out.append(fp)
                                else:
                                    for fp in sorted(ydir.glob("*.nc")):
                                        if L == "L2P":
                                            ts = self._parse_yyyymmdd(fp.name)
                                            if ts is None or ts < t0 or ts > t1:
                                                continue
                                        else:
                                            ym = self._parse_yyyymm(fp.name)
                                            if ym and not self._month_overlap(ym[0], ym[1], t0, t1):
                                                continue
                                        out.append(fp)
        # ---------------- AWI ----------------
        if not insts or "AWI" in insts:
            if not level_set:
                awi_levels: List[str] = ["l2p_release", "l3cp_release"]
            else:
                # map friendly ESA-like inputs to AWI release names
                canonical = []
                for L in level_set:
                    l = L.lower()
                    if l in ("l2p_release", "l3cp_release"):
                        canonical.append(l)
                    elif l == "l2p":
                        canonical.append("l2p_release")
                    elif l == "l3c":
                        canonical.append("l3cp_release")
                awi_levels = sorted(set(canonical)) or ["l2p_release", "l3cp_release"]
            LHEM = hem_lower  # AWI uses lower-case hemisphere dirs
            for Lraw in awi_levels:
                Ldir = Path(root_awi,Lraw,LHEM)
                if not Ldir.exists():
                    continue
                for sdir in sorted(p for p in Ldir.iterdir() if p.is_dir()):
                    sensor = sdir.name.lower()
                    if sensor_set and sensor not in sensor_set:
                        continue
                    for ydir in sorted(sdir.glob("[12][0-9][0-9][0-9]")):
                        y = int(ydir.name)
                        mdirs = list(sorted(ydir.glob("[01][0-9]")))
                        if mdirs:
                            for mdir in mdirs:
                                m = int(mdir.name)
                                if not self._month_overlap(y, m, t0, t1):
                                    continue
                                out.extend(sorted(mdir.glob("*.nc")))
                        else:
                            # year-only dir; filter by yyyymm in filenames when available
                            for fp in sorted(ydir.glob("*.nc")):
                                ym = self._parse_yyyymm(fp.name)
                                if ym and not self._month_overlap(ym[0], ym[1], t0, t1):
                                    continue
                                out.append(fp)
        # de-dup & sort
        return sorted(set(out))

    # ------------------------- map grids (LAEA/EASE-like) ------------------------- #
    # ------------------------- file path indexing & QA helpers ------------------------- #
    @staticmethod
    def laea_area_def(hem   : str,
                      nx    : int,
                      ny    : int,
                      dx_km : float,
                      lon0  : float = 0.0):
        """
        Build a Lambert Azimuthal Equal-Area grid centered on the pole.

        Parameters
        ----------
        hem : {'sh','nh'}
            Hemisphere.
        nx, ny : int
            Grid shape (x, y).
        dx_km : float
            Grid spacing in kilometers (projected meters).
        lon0 : float, default 0.0
            Central longitude (°E). ESA products commonly use 0.

        Returns
        -------
        pyresample.geometry.AreaDefinition
        """
        from pyresample import geometry
        if geometry is None:
            raise ImportError("pyresample is required: `pip install pyresample`")
        hem       = str(hem).lower()
        lat0      = -90.0 if hem == "sh" else 90.0
        proj_id   = f"laea_{hem}_{int(dx_km)}km_{nx}x{ny}"
        proj_dict = {"proj"    : "laea",
                     "lat_0"   : lat0,
                     "lon_0"   : lon0,
                     "a"       : 6378137.0,
                     "b"       : 6356752.314245,
                     "units"   : "m",
                     "no_defs" : True}
        # extent centered on the pole
        half_w      = (nx * dx_km * 1000.0) / 2.0
        half_h      = (ny * dx_km * 1000.0) / 2.0
        area_extent = (-half_w, -half_h, half_w, half_h)
        area_id     = proj_id
        description = f"LAEA {hem.upper()} {dx_km}km {nx}x{ny} lon0={lon0}"
        return geometry.AreaDefinition(area_id, description, proj_id, proj_dict, nx, ny, area_extent)

    @staticmethod
    def ease2_sh_50km():
        """Convenience grid roughly matching AWI SH L3C (216×216 @ ~50 km, lon0=0)."""
        return SeaIceObservations.laea_area_def("sh", nx=216, ny=216, dx_km=50.0, lon0=0.0)

    @staticmethod
    def ease2_nh_25km():
        """Convenience grid for NH (432×432 @ ~25 km, lon0=0)."""
        return SeaIceObservations.laea_area_def("nh", nx=432, ny=432, dx_km=25.0, lon0=0.0)

    
    def index_satellite_paths(self, paths: Iterable[Path]) -> pd.DataFrame:
        """
        Build a tidy index of candidate satellite files with parsed dates for QA/QC and selection.

        This function inspects each file name and classifies it as:
        - "daily"   : contains a YYYYMMDD token (date parsed as that day)
        - "monthly" : contains a YYYYMM token (date set to first of that month)
        - "unknown" : no parseable token (date set to NaT)

        Parameters
        ----------
        paths : iterable of pathlib.Path
            Input file paths.

        Returns
        -------
        pandas.DataFrame
            A stable-sorted table with columns:
            - path  : Path
            - kind  : {"daily","monthly","unknown"}
            - date  : UTC Timestamp (daily) or first-of-month (monthly) or NaT
            - year, month, day : ints (or None)

        Notes
        -----
        Sorting is stable and deterministic to make downstream "first"/"last" selection
        repeatable.
        """
        rows = []
        for p in map(Path, paths):
            name = p.name
            ts = self._parse_yyyymmdd(name)
            if ts is not None:
                rows.append({"path": p, "kind": "daily",
                             "date": ts, "year": ts.year, "month": ts.month, "day": ts.day})
                continue
            ym = self._parse_yyyymm(name)
            if ym:
                y, m = ym
                rows.append({"path": p, "kind": "monthly",
                             "date": pd.Timestamp(year=y, month=m, day=1, tz="UTC"),
                             "year": y, "month": m, "day": None})
            else:
                rows.append({"path": p, "kind": "unknown",
                             "date": pd.NaT, "year": None, "month": None, "day": None})
        df = pd.DataFrame(rows)
        # Stable sort so 'first'/'last' are deterministic
        df = df.sort_values(by          = ["kind", "date", "year", "month", "day", "path"],
                            na_position = "last",
                            kind        = "mergesort").reset_index(drop=True)
        return df

    def dedup_daily_one_file_per_day(self, df: pd.DataFrame, prefer: str = "last") -> pd.DataFrame:
        """
        Select exactly one daily file per UTC day from an indexed path table.

        Parameters
        ----------
        df : pandas.DataFrame
            Output of `index_satellite_paths(...)`.
        prefer : {"first","last"}, default "last"
            If multiple files map to the same day, select the first/last after sorting.

        Returns
        -------
        pandas.DataFrame
            Subset of `df` restricted to kind == "daily" with unique 'date'.

        Notes
        -----
        - Non-daily entries are discarded.
        - The returned frame preserves the same columns as the input `df`.
        """
        dfd = df[df["kind"] == "daily"].copy()
        if dfd.empty:
            return dfd
        dfd  = dfd.sort_values(["date", "path"])
        pick = {"first": "first", "last": "last"}[prefer]
        dfd  = dfd.groupby("date", as_index=False, sort=True).agg(pick)
        return dfd
    
    def dedup_monthly_one_file_per_month(self, df: pd.DataFrame, prefer: str = "last") -> pd.DataFrame:
        """
        Select exactly one file per (year, month) from an indexed path table.

        This is intended for monthly products, but will also accept daily files that encode
        a YYYYMM token.

        Parameters
        ----------
        df : pandas.DataFrame
            Output of `index_satellite_paths(...)`.
        prefer : {"first","last"}, default "last"
            If multiple files map to the same (year, month), select the first/last after sorting.

        Returns
        -------
        pandas.DataFrame
            Subset with unique (year, month) combinations.

        Notes
        -----
        Rows without parsed (year, month) are dropped.
        """
        dfm = df[df["kind"].isin(["monthly","daily"])].copy()
        if dfm.empty:
            return dfm
        # For daily entries we keep (year, month) for grouping
        dfm["yymm"] = list(zip(dfm["year"], dfm["month"]))
        dfm = dfm.dropna(subset=["yymm"])
        dfm = dfm.sort_values(["year","month","path"])
        pick = {"first":"first","last":"last"}[prefer]
        dfm = dfm.groupby(["year","month"], as_index=False, sort=True).agg(pick)
        return dfm

    # ------------------------- LEVEL 2 (DAILY UN-GRIDDED) ------------------------- #
    def _l2p_weights(self,
                     index_array,                 # masked or ndarray; shape (k,n_out) OR (n_out,k)
                     dist_array,                  # masked or ndarray; same shape/orientation
                     valid_input_index,           # (Nvalid_in,)
                     valid_output_index,          # (n_out,) ints OR bool mask of length ny*nx
                     swath_SIT,                   # (Nswath,)
                     swath_unc,                   # (Nswath,)
                     out_shape : Tuple[int, int],  # (ny, nx)
                     roi_m     : float,
                     sigma_m   : Optional[float] = None):
        """
        Combine neighbour swath samples into gridded cell statistics using uniform or Gaussian weights.

        This helper is designed to consume the outputs of:
            pyresample.kd_tree.get_neighbour_info(...)

        and reduce swath SIT and uncertainty into per-grid-cell:
        - weighted mean SIT
        - propagated standard error (SE)
        - weighted mean of per-sample uncertainty
        - count of valid contributing neighbours

        Parameters
        ----------
        index_array : array-like (masked or ndarray)
            Neighbour indices from `get_neighbour_info`. Shape is either (k, n_out) or
            (n_out, k), where k is the number of neighbours and n_out is the number of
            valid output cells.
        dist_array : array-like (masked or ndarray)
            Neighbour distances (meters) aligned with `index_array`.
        valid_input_index : array-like of int
            Indices selecting the valid subset of swath input points used by the KD-tree.
        valid_output_index : array-like of int or bool
            Output-cell selector from `get_neighbour_info`. If boolean, length must equal
            ny*nx; if integer, interpreted as flat indices into the target grid.
        swath_SIT : array-like of float
            Swath sea-ice thickness samples corresponding to the full swath arrays
            (pre-filtering). Units: meters.
        swath_unc : array-like of float
            Swath uncertainty samples (per-sample standard error or equivalent). Units: meters.
        out_shape : tuple[int, int]
            Target grid shape as (ny, nx).
        roi_m : float
            Radius of influence used for neighbour selection (meters).
        sigma_m : float or None, default None
            If provided, apply Gaussian weights:
                w = exp(-d^2 / (2*sigma^2))
            Otherwise, use uniform weights for valid neighbours.

        Returns
        -------
        tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray, numpy.ndarray]
            (SIT, SIT_unc, u_cell, n_cell) where:
            - SIT      : (ny, nx) float32, weighted mean thickness [m]
            - SIT_unc  : (ny, nx) float32, propagated SE [m]
            - u_cell   : (ny, nx) float32, weighted mean of per-sample uncertainty [m]
            - n_cell   : (ny, nx) int32, number of valid neighbours per cell

        Notes
        -----
        - Neighbours outside `roi_m` are excluded.
        - Masked neighbour entries are treated as non-existent.
        - SE propagation uses: sqrt(sum(w^2 * u^2)) / sum(w).
        """
        ny, nx = out_shape
        # reduce sources to KD-tree's valid input subset
        src_s = np.asarray(swath_SIT, dtype=np.float64)[valid_input_index]
        src_u = np.asarray(swath_unc, dtype=np.float64)[valid_input_index]
        # normalise valid_output_index → integer flat indices
        voi            = np.asarray(valid_output_index)
        out_idx        = np.flatnonzero(voi) if voi.dtype == bool else voi.astype(np.int64, copy=False)
        n_out_expected = out_idx.size
        # fill masked arrays to ndarrays + build "exists" mask
        idx          = np.ma.filled(index_array, 0).astype(np.int64, copy=False)
        dst          = np.ma.filled(dist_array,  np.inf).astype(np.float64, copy=False)
        neigh_exists = ~np.ma.getmaskarray(index_array)
        # orient to (k, n_out) i.e. neighbours along axis 0
        if idx.shape[1] == n_out_expected:
            pass
        elif idx.shape[0] == n_out_expected:
            idx, dst, neigh_exists = idx.T, dst.T, neigh_exists.T
        else:
            raise ValueError(f"Neighbour array shape {idx.shape} incompatible with n_out={n_out_expected}")
        k, n_out = idx.shape
        if n_out != n_out_expected:
            raise ValueError(f"Internal n_out mismatch: {n_out} vs {n_out_expected}")
        # gather neighbour values from reduced swath
        v_s = np.take(src_s, idx, mode="clip")    # (k, n_out)
        v_u = np.take(src_u, idx, mode="clip")    # (k, n_out)
        # validity
        within_roi = np.isfinite(dst) & (dst <= roi_m)
        val_ok     = neigh_exists & within_roi & np.isfinite(v_s) & np.isfinite(v_u)
        # weights
        if sigma_m is None:
            w = np.where(val_ok, 1.0, 0.0)
        else:
            w = np.where(val_ok, np.exp(-(dst**2) / (2.0 * sigma_m**2)), 0.0)
        # reductions along neighbour axis
        wsum     = w.sum(axis=0)
        sitnum   = (w * v_s).sum(axis=0)
        umean    = (w * v_u).sum(axis=0)
        senum    = ((w**2) * (v_u**2)).sum(axis=0)
        n_cnt    = val_ok.sum(axis=0).astype("int32")
        sit_mean = np.divide(sitnum, wsum, out=np.full_like(wsum, np.nan, dtype=float), where=wsum>0)
        u_mean   = np.divide(umean,  wsum, out=np.full_like(wsum,  np.nan, dtype=float), where=wsum>0)
        se_mean  = np.divide(np.sqrt(senum), wsum, out=np.full_like(wsum, np.nan, dtype=float), where=wsum>0)
        # scatter back to flat arrays
        nflat        = ny * nx
        SIT_flat     = np.full(nflat, np.nan, dtype="float32")
        SIT_unc_flat = np.full(nflat, np.nan, dtype="float32")
        U_cell_flat  = np.full(nflat, np.nan, dtype="float32")
        N_cell_flat  = np.zeros(nflat, dtype="int32")
        SIT_flat[out_idx]     = sit_mean.astype("float32", copy=False)
        SIT_unc_flat[out_idx] = se_mean.astype("float32", copy=False)
        U_cell_flat[out_idx]  = u_mean.astype("float32", copy=False)
        N_cell_flat[out_idx]  = n_cnt
        return (SIT_flat.reshape(ny, nx),
                SIT_unc_flat.reshape(ny, nx),
                U_cell_flat.reshape(ny, nx),
                N_cell_flat.reshape(ny, nx))

    # ------------------------- main driver ------------------------- #
    def l2p_to_daily_grid_pyresample(self,
                                     paths             : List[Path],
                                     hem               : Optional[str] = None,
                                     area_def = None,
                                     roi_km            : float = 75.0,
                                     neighbours        : int = 16,
                                     epsilon           : float = 0.0,
                                     gaussian_sigma_km : Optional[float] = None,
                                     lat_cut           : Optional[float] = None,
                                     time_at_noon      : bool = True) -> xr.Dataset:
        """     
        Resample ESA/AWI **L2P swath** files to a **daily gridded** product on an
        equal-area polar grid using pyresample nearest-neighbour sets.

        The algorithm:
        1) Group all valid swath samples by UTC day (per file) using `time` variable.
        2) For each day:
           - Build a SwathDefinition from lon/lat points.
           - Use `pyresample.kd_tree.get_neighbour_info` to find up to `neighbours`
             candidates within `roi_km` of each grid cell.
           - Combine neighbour SIT & uncertainty with uniform or Gaussian weights.
        3) Assemble (time, y, x) arrays and compute hemispheric daily means with
           cos(latitude) weights.

        Parameters
        ----------
        paths : list of Path
            L2P NetCDF files. Each must include variables:
            `sea_ice_thickness`, `sea_ice_thickness_uncertainty`, `lon`, `lat`, `time`.
        hem : {'sh','nh'}, optional
            Hemisphere. Defaults to `self.hemisphere_dict['abbreviation']` when present,
            else 'sh'.
        area_def : pyresample.geometry.AreaDefinition, optional
            Target equal-area grid. If None, uses `ease2_sh_50km()` for SH or
            `ease2_nh_25km()` for NH.
        roi_km : float, default 75.0
            Radius of influence for neighbour search (kilometres).
        neighbours : int, default 16
            Max neighbours per grid cell from KD-tree.
        epsilon : float, default 0.0
            Search tolerance passed to pyresample (trade accuracy vs. speed).
        gaussian_sigma_km : float or None, default None
            If set, apply Gaussian weights with this sigma (km); otherwise uniform.
        lat_cut : float or None, default None
            Optional extra swath mask: for SH keep `lat <= lat_cut`; for NH keep
            `lat >= lat_cut`. Useful for excluding marginal seas.
        time_at_noon : bool, default True
            If True, time coordinate is day+12:00Z; else midnight.

        Returns
        -------
        xarray.Dataset
            Dimensions: time, y, x
            Coordinates: lon(y,x), lat(y,x), time
            Data variables:
                SIT(time,y,x)       [m]
                SIT_unc(time,y,x)   [m] propagated standard error
                u_cell(time,y,x)    [m] per-cell uncertainty mean (neighbour-weighted)
                n_cell(time,y,x)    [#] count of contributing neighbours per cell
                SIT_hem(time)       [m] hemispheric mean
                SIT_hem_unc(time)   [m] propagated SE for hemispheric mean
                SIT_hem_u(time)     [m] hemispheric mean of per-cell uncertainty
                SIT_hem_n(time)     [#] total contributing neighbours (sum over grid)

        Notes
        -----
        - Hemispheric means use cos(lat) weights on the target grid after resampling.
        - Propagated SE follows: sqrt(sum(w^2 * se^2)) / sum(w).
        - Files missing required variables are skipped (warning).
        """
        from pyresample import geometry, kd_tree
        if geometry is None or kd_tree is None:
            raise ImportError("pyresample is required: `pip install pyresample`")
        if not paths:
            raise ValueError("No L2P paths provided")
        hem = (hem or self.hemisphere_dict.get("abbreviation", "sh")).lower()
        if area_def is None:
            area_def = self.ease2_sh_50km() if hem == "sh" else self.ease2_nh_25km()
        # grid lon/lat for metadata & hemispheric weights
        lons2d, lats2d = area_def.get_lonlats()
        ny, nx = area_def.shape
        # collect swath samples by UTC day
        by_day: Dict[pd.Timestamp, Dict[str, List[np.ndarray]]] = {}
        for fp in sorted(map(Path, paths)):
            try:
                with xr.open_dataset(fp, decode_times=True, chunks={}) as ds:
                    required = ("sea_ice_thickness", "sea_ice_thickness_uncertainty", "lon", "lat", "time")
                    if not all(v in ds.variables for v in required):
                        # skip quietly if this file is not SIT L2P
                        continue
                    sit_da = self._get_first(ds, ["sea_ice_thickness", "SIT", "sit"])
                    unc_da = self._get_first(ds, ["sea_ice_thickness_uncertainty", "SIT_unc", "uncertainty"])
                    lon_da = self._get_first(ds, ["lon", "longitude"])
                    lat_da = self._get_first(ds, ["lat", "latitude"])
                    tim_da = self._get_first(ds, ["time", "acquisition_time"])
                    if any(v is None for v in (sit_da, unc_da, lon_da, lat_da, tim_da)):
                        # Skip quietly if not an SIT L2P file
                        continue
                    sit = sit_da.values.astype("float64", copy=False)
                    unc = unc_da.values.astype("float64", copy=False)
                    lon = lon_da.values.astype("float64", copy=False)
                    lon = self._match_lon_range(lon, lons2d)
                    lat = lat_da.values.astype("float64", copy=False)
                    tt  = pd.to_datetime(tim_da.values)
            except Exception:
                # unreadable file: skip
                continue
            valid = np.isfinite(sit) & np.isfinite(unc) & np.isfinite(lon) & np.isfinite(lat)
            if lat_cut is not None:
                valid &= (lat <= lat_cut) if hem == "sh" else (lat >= lat_cut)
            if not np.any(valid):
                self.logger.info(f"[skip] {fp.name}: finite(sit)={np.isfinite(sit).sum()} "
                                 f"finite(lon)={np.isfinite(lon).sum()} finite(lat)={np.isfinite(lat).sum()}")
                continue
            sit, unc, lon, lat, tt = sit[valid], unc[valid], lon[valid], lat[valid], tt[valid]
            days = pd.to_datetime(tt.floor("D").to_numpy())
            for d in np.unique(days):
                sel = (days == d)
                if not np.any(sel):
                    continue
                dkey = pd.Timestamp(d, tz="UTC")
                rec = by_day.setdefault(dkey, dict(sit=[], unc=[], lon=[], lat=[]))
                rec["sit"].append(sit[sel])
                rec["unc"].append(unc[sel])
                rec["lon"].append(lon[sel])
                rec["lat"].append(lat[sel])
        # empty output if nothing usable
        if not by_day:
            return xr.Dataset(
                coords=dict(time = ("time", np.array([], dtype="datetime64[ns]")),
                            y    = ("y", np.arange(ny, dtype=int)),
                            x    = ("x", np.arange(nx, dtype=int)),
                            lon  = (("y","x"), lons2d),
                            lat  = (("y","x"), lats2d)))
        roi_m = roi_km * 1000.0
        sigma_m = None if gaussian_sigma_km is None else gaussian_sigma_km * 1000.0
        # allocate outputs
        days_sorted = sorted(by_day.keys())
        T           = len(days_sorted)
        SIT         = np.full((T, ny, nx), np.nan, dtype="float32")
        SIT_unc     = np.full((T, ny, nx), np.nan, dtype="float32")
        U_cell      = np.full((T, ny, nx), np.nan, dtype="float32")
        N_cell      = np.zeros((T, ny, nx), dtype="int32")
        # time coordinate
        if time_at_noon:
            tcoord = [pd.Timestamp(d.date()).tz_localize("UTC") + pd.Timedelta(hours=12) for d in days_sorted]
        else:
            tcoord = [pd.Timestamp(d.date()).tz_localize("UTC") for d in days_sorted]
        tcoord = pd.to_datetime(tcoord).tz_localize(None).to_numpy()
        # resample per day
        for it, d in enumerate(days_sorted):
            s     = np.concatenate(by_day[d]["sit"])
            u     = np.concatenate(by_day[d]["unc"])
            lo    = np.concatenate(by_day[d]["lon"])
            la    = np.concatenate(by_day[d]["lat"])
            swath = geometry.SwathDefinition(lons=lo, lats=la)
            valid_in_idx, valid_out_idx, index_array, dist_array = kd_tree.get_neighbour_info(source_geo_def      = swath,
                                                                                              target_geo_def      = area_def,
                                                                                              radius_of_influence = roi_m,
                                                                                              neighbours          = neighbours,
                                                                                              epsilon             = epsilon)
            if getattr(valid_out_idx, "size", 0) == 0:
                continue
            sit2d, se2d, u2d, n2d = self._l2p_weights(index_array        = index_array,
                                                      dist_array         = dist_array,
                                                      valid_input_index  = valid_in_idx,
                                                      valid_output_index = valid_out_idx,
                                                      swath_SIT          = s,
                                                      swath_unc          = u,
                                                      out_shape          = (ny, nx),
                                                      roi_m              = roi_m,
                                                      sigma_m            = sigma_m)
            SIT[it,:,:]     = sit2d
            SIT_unc[it,:,:] = se2d
            U_cell[it,:,:]  = u2d
            N_cell[it,:,:]  = n2d
        # hemispheric daily aggregates with cos(lat) weights
        W       = np.cos(np.deg2rad(lats2d)).astype("float64")          # (y,x)
        valid   = np.isfinite(SIT) & np.isfinite(W)[None, :, :]
        S       = np.where(valid, SIT,     0.0).astype("float64", copy=False)
        U       = np.where(valid, SIT_unc, 0.0).astype("float64", copy=False)
        Uc      = np.where(valid, U_cell,  0.0).astype("float64", copy=False)
        W3      = W[None, :, :]
        wsum    = (W3 * valid).sum(axis=(1, 2))                   # (time,)
        sit_num = (S  * W3).sum(axis=(1, 2))
        ume_num = (Uc * W3).sum(axis=(1, 2))
        se_num  = ((W3**2) * (U**2)).sum(axis=(1, 2))
        SIT_hem     = sit_num / np.where(wsum > 0, wsum, np.nan)
        SIT_hem_u   = ume_num / np.where(wsum > 0, wsum, np.nan)
        SIT_hem_unc = np.sqrt(se_num) / np.where(wsum > 0, wsum, np.nan)
        SIT_hem_n   = N_cell.sum(axis=(1, 2)).astype("int32")
        ds = xr.Dataset(data_vars = dict(SIT         = (("time","y","x"), SIT),
                                         SIT_unc     = (("time","y","x"), SIT_unc),
                                         u_cell      = (("time","y","x"), U_cell),
                                         n_cell      = (("time","y","x"), N_cell),
                                         SIT_hem     = (("time",), SIT_hem.astype("float32")),
                                         SIT_hem_unc = (("time",), SIT_hem_unc.astype("float32")),
                                         SIT_hem_u   = (("time",), SIT_hem_u.astype("float32")),
                                         SIT_hem_n   = (("time",), SIT_hem_n),),
                        coords = dict(time = ("time", tcoord),
                                      y    = ("y", np.arange(ny, dtype=int)),
                                      x    = ("x", np.arange(nx, dtype=int)),
                                      lon  = (("y","x"), lons2d.astype("float64")),
                                      lat  = (("y","x"), lats2d.astype("float64"))),
                        attrs = dict(note = "L2P swath resampled to target grid via pyresample neighbour info; "
                                            "Gaussian weighting if gaussian_sigma_km set; SE propagated as sqrt(sum(w^2 u^2))/sum(w)."))
        ds["SIT"].attrs.update(dict(standard_name="sea_ice_thickness", units="m"))
        ds["SIT_unc"].attrs.update(dict(standard_name="sea_ice_thickness standard_error", units="m"))
        return ds

    # ------------------------- one-call convenience ------------------------- #
    def make_daily_gridded_SIT(self,
                               dt0_str           : str,
                               dtN_str           : str,
                               hemisphere        : Optional[str] = None,
                               institutions      : Optional[Iterable[str]] = ("ESA", "AWI"),
                               sensors           : Optional[Iterable[str]] = None,
                               levels            : Optional[Iterable[str]] = None,
                               versions          : Optional[Iterable[str]] = None,
                               root_esa          : Optional[Path] = None,
                               root_awi          : Optional[Path] = None,
                               roi_km            : float = 75.0,
                               neighbours        : int = 16,
                               epsilon           : float = 0.0,
                               gaussian_sigma_km : Optional[float] = None,
                               lat_cut           : Optional[float] = None,
                               time_at_noon      : bool = True,
                               prefer            : str = "last",
                               area_def          = None) -> xr.Dataset:
        """
        Create a daily gridded sea-ice thickness (SIT) dataset over a time window in one call.

        This convenience wrapper performs:
        1) File discovery via `build_sea_ice_satellite_paths(...)` across ESA and/or AWI layouts.
        2) Indexing and de-duplication to one candidate file per UTC day.
        3) Swath-to-grid resampling with `l2p_to_daily_grid_pyresample(...)`.

        Parameters
        ----------
        dt0_str, dtN_str : str
            Start/end of the time window (inclusive). Any pandas-parsable date string.
        hemisphere : {"sh","nh"} or None, optional
            Hemisphere selector; defaults to `self.hemisphere_dict['abbreviation']` when present.
        institutions : iterable of str or None, default ("ESA","AWI")
            Institutions to search (case-insensitive).
        sensors, levels, versions : iterable of str or None, optional
            Optional filters passed to `build_sea_ice_satellite_paths(...)`.
        root_esa, root_awi : pathlib.Path or None, optional
            Root directories. If None, falls back to configured defaults (module globals or class config).
        roi_km : float, default 75.0
            Radius of influence for neighbour search (km).
        neighbours : int, default 16
            Maximum neighbours per target grid cell.
        epsilon : float, default 0.0
            KD-tree epsilon tolerance (speed/accuracy trade-off).
        gaussian_sigma_km : float or None, optional
            If set, apply Gaussian distance weights with this sigma (km). If None, uniform weights.
        lat_cut : float or None, optional
            Additional latitude filter applied to swath points (hemisphere-aware).
        time_at_noon : bool, default True
            If True, output time coordinate is day+12:00Z; otherwise midnight.
        prefer : {"first","last"}, default "last"
            If multiple files map to a day, choose first/last after deterministic sorting.
        area_def : pyresample.geometry.AreaDefinition or None, optional
            Target grid. If None, uses internal EASE/LAEA-like convenience grids.

        Returns
        -------
        xarray.Dataset
            Daily gridded dataset with dimensions (time, y, x) including:
            - SIT, SIT_unc, u_cell, n_cell (gridded)
            - SIT_hem, SIT_hem_unc, SIT_hem_u, SIT_hem_n (hemispheric aggregates)

        See Also
        --------
        build_sea_ice_satellite_paths
        index_satellite_paths
        dedup_daily_one_file_per_day
        l2p_to_daily_grid_pyresample
        """
        # prefer configured roots if not supplied
        if root_esa is None:
            root_esa = globals().get("ROOT_ESA_DEF", None)
        if root_awi is None:
            root_awi = globals().get("ROOT_AWI_DEF", None)
        paths = self.build_sea_ice_satellite_paths(dt0_str      = dt0_str,
                                                   dtN_str      = dtN_str,
                                                   hemisphere   = hemisphere,
                                                   institutions = institutions,
                                                   sensors      = sensors,
                                                   levels       = levels,
                                                   versions     = versions,
                                                   root_esa     = root_esa,
                                                   root_awi     = root_awi)
        df          = self.index_satellite_paths(paths)
        df_daily    = self.dedup_daily_one_file_per_day(df, prefer=prefer)
        daily_paths = df_daily["path"].tolist()
        return self.l2p_to_daily_grid_pyresample(paths             = daily_paths,
                                                 hem               = hemisphere,
                                                 area_def          = area_def,
                                                 roi_km            = roi_km,
                                                 neighbours        = neighbours,
                                                 epsilon           = epsilon,
                                                 gaussian_sigma_km = gaussian_sigma_km,
                                                 lat_cut           = lat_cut,
                                                 time_at_noon      = time_at_noon)

    # ------------------------- LEVEL 3 (MONTLY GRIDDED) ------------------------- #
    def l3c_paths_to_monthly_grid(self,
                                  paths              : list[Path],
                                  hem                : str   | None = None,
                                  sic_thresh_percent : float | None = 15.0,
                                  mask_strategy      : str  = "none",    # {"none","quality","status","both"}
                                  include_snow       : bool = True,
                                  include_flags      : bool = True,
                                  min_obs            : int  = 1,               # if an "n_observations" variable exists
                                  time_midpoint      : bool = True) -> xr.Dataset:
        """
        Assemble ESA L3C / AWI l3cp_release files (already on LAEA grids) into a
        single (time, y, x) dataset with hemispheric monthly means.

        This function mirrors the L2P driver logic, but **no resampling** is done:
        it reads gridded variables straight from the files and standardizes dims.

        Parameters
        ----------
        paths : list[Path]
            Monthly L3 NetCDF paths (ESA L3C or AWI l3cp_release).
        hem : {'sh','nh'}, optional
            Hemisphere hint (controls sign in optional lat_cut if you add one later).
            Defaults to `self.hemisphere_dict['abbreviation']` if present.
        sic_thresh_percent : float or None, default 15.0
            If set, mask cells with sea_ice_concentration < threshold (percent).
        mask_strategy : {"none","quality","status","both"}, default "none"
            - "quality" keeps cells with quality_flag == 0
            - "status"  keeps cells with status_flag  == 0
            - "both"    requires both conditions
        include_snow : bool, default True
            Include snow_depth and snow_depth_uncertainty when present.
        include_flags : bool, default True
            Include flag fields (quality_flag, status_flag, region_code) when present.
        min_obs : int, default 1
            If a count variable exists (e.g., "n_observations"), require count >= min_obs.
        time_midpoint : bool, default True
            If time_bnds variable exists, use midpoint as the time coordinate.

        Returns
        -------
        xarray.Dataset
            Dimensions: time, y, x
            Coordinates: lon(y,x), lat(y,x), time
            Data variables:
                SIT(time,y,x)       [m]
                SIT_unc(time,y,x)   [m] per-cell uncertainty (SE)
                SIC(time,y,x)       [%] if present
                snow_depth(time,y,x) [m] if requested/present
                snow_depth_unc(time,y,x) [m] if requested/present
                quality_flag/status_flag/region_code if requested/present
                SIT_hem(time)       [m] hemispheric mean
                SIT_hem_unc(time)   [m] propagated SE = sqrt(sum(w^2 * se^2))/sum(w)
                SIT_hem_u(time)     [m] hemispheric mean of per-cell uncertainty
                SIT_hem_n(time)     [#] count of valid contributing cells
        """
        if not paths:
            raise ValueError("No L3 paths provided")
        hem = (hem or self.hemisphere_dict.get("abbreviation", "sh")).lower()
        # Store per-file outputs here then concat along time at the end
        out_list: list[xr.Dataset] = []
        lat_template = None
        lon_template = None
        for fp in sorted(map(Path, paths)):
            try:
                ds = xr.open_dataset(fp, decode_times=True)
            except Exception:
                continue
            # Flexible getters
            def G(names):
                return self._get_first(ds, names)
            # Dim normalization (yc/xc or y/x)
            ydim = "yc" if "yc" in ds.dims else ("y" if "y" in ds.dims else None)
            xdim = "xc" if "xc" in ds.dims else ("x" if "x" in ds.dims else None)
            if ydim is None or xdim is None:
                ds.close()
                continue
            ds_work = ds.rename({ydim: "y", xdim: "x"})
            # Coordinates
            lat2d = G(["lat"])
            lon2d = G(["lon"])
            if lat2d is None or lon2d is None:
                ds.close()
                continue
            # Persist a lat/lon template (assumed constant across files)
            lat_template = lat2d.values if lat_template is None else lat_template
            lon_template = lon2d.values if lon_template is None else lon_template
            # Variables of interest
            SIT   = G(["sea_ice_thickness", "SIT", "sit"])
            SIT_U = G(["sea_ice_thickness_uncertainty", "SIT_unc", "sit_unc", "uncertainty"])
            SIC   = G(["sea_ice_concentration", "sic"])
            SD    = G(["snow_depth", "snow_thickness"])
            SDU   = G(["snow_depth_uncertainty", "snow_thickness_uncertainty"])
            QF    = G(["quality_flag"])
            SF    = G(["status_flag"])
            RC    = G(["region_code"])
            NOBS  = G(["n_observations", "num_observations", "count"])
            # Required: SIT & SIT_U at least
            if SIT is None or SIT_U is None:
                ds.close(); continue
            # Time coordinate (often size 1 per file, but robust to >1)
            time = G(["time"])
            tbnd = G(["time_bnds", "time_bounds"])
            if time is None:
                ds.close(); continue
            tvals = pd.to_datetime(time.values)
            if tvals.ndim == 0:
                tvals = tvals[None]
            if time_midpoint and (tbnd is not None):
                try:
                    tb = xr.DataArray(tbnd).compute().values  # (time, nv)
                    tmid = np.mean(tb, axis=1)
                    tvals = pd.to_datetime(tmid, unit="s", origin="unix")
                except Exception:
                    # fall back to time center
                    pass
            # Build a mask
            base_mask = np.isfinite(SIT) & np.isfinite(SIT_U) & np.isfinite(lat2d) & np.isfinite(lon2d)
            if SIC is not None and sic_thresh_percent is not None:
                base_mask &= np.isfinite(SIC) & (SIC >= float(sic_thresh_percent))
            if NOBS is not None and isinstance(min_obs, (int, float)) and min_obs > 1:
                base_mask &= (NOBS >= int(min_obs))
            if mask_strategy in ("quality", "both") and QF is not None:
                base_mask &= (QF == 0)
            if mask_strategy in ("status", "both") and SF is not None:
                base_mask &= (SF == 0)
            mask_strategy = (mask_strategy or "none").lower()
            # Ensure 3D shape: (time, y, x)
            def _as3(da):
                v = da.transpose(..., "yc", "xc").values
                if v.ndim == 2:
                    v = v[None, ...]
                return v
            sit3  = _as3(SIT)
            sunc3 = _as3(SIT_U)
            mask3 = base_mask
            if mask3.ndim == 2:
                mask3 = mask3[None, ...]
            if SIC is not None:
                sic3 = _as3(SIC)
            if SD is not None and include_snow:
                sd3 = _as3(SD)
            if SDU is not None and include_snow:
                sdu3 = _as3(SDU)
            if QF is not None and include_flags:
                qf3 = _as3(QF)
            if SF is not None and include_flags:
                sf3 = _as3(SF)
            if RC is not None and include_flags:
                rc3 = _as3(RC)
            # Apply mask: keep NaNs outside valid cells
            S  = np.where(mask3, sit3,  np.nan).astype("float32", copy=False)
            U  = np.where(mask3, sunc3, np.nan).astype("float32", copy=False)
            if SIC is not None:
                SICv = np.where(mask3, sic3, np.nan).astype("float32", copy=False)
            if SD is not None and include_snow:
                SDv  = np.where(mask3, sd3,  np.nan).astype("float32", copy=False)
            if SDU is not None and include_snow:
                SDUv = np.where(mask3, sdu3, np.nan).astype("float32", copy=False)
            # Hemispheric weighting
            W2      = np.cos(np.deg2rad(lat2d.values)).astype("float64")
            W3      = W2[None, :, :]
            valid   = np.isfinite(S)
            wsum    = (W3 * valid).sum(axis=(1, 2))
            sit_num = (np.where(valid, S, 0.0)  * W3).sum(axis=(1, 2))
            u_num   = (np.where(valid, U, 0.0)  * W3).sum(axis=(1, 2))
            se_num  = ((W3**2) * np.where(valid, U, 0.0)**2).sum(axis=(1, 2))
            SIT_hem     = sit_num / np.where(wsum > 0, wsum, np.nan)
            SIT_hem_u   = u_num   / np.where(wsum > 0, wsum, np.nan)
            SIT_hem_unc = np.sqrt(se_num) / np.where(wsum > 0, wsum, np.nan)
            SIT_hem_n   = valid.sum(axis=(1, 2)).astype("int32")
            # ---- NEW: concentration-weighted hemispheric mean (SIC weights) ----
            SIT_wgt     = np.full((S.shape[0],), np.nan, dtype="float64")
            SIT_wgt_u   = np.full((S.shape[0],), np.nan, dtype="float64")
            SIT_wgt_unc = np.full((S.shape[0],), np.nan, dtype="float64")
            SIT_wgt_n   = np.zeros((S.shape[0],), dtype="int32")
            if SIC is not None:
                # Convert SIC to fraction [0..1] for correct uncertainty propagation.
                sic_units = (SIC.attrs.get("units","") or "").lower()
                sic3_raw = _as3(SIC)
                if ("percent" in sic_units) or ("%" in sic_units):
                    Wsic = sic3_raw.astype("float64") / 100.0
                else:
                    # heuristic fallback
                    vmax = np.nanmax(sic3_raw)
                    Wsic = (sic3_raw / 100.0) if np.isfinite(vmax) and vmax > 1.5 else sic3_raw
                    Wsic = Wsic.astype("float64")
                # Honor the same base mask used for S/U (keeps SIC >= thresh, flags, etc.)
                # mask3 is (time,y,x); if it's 2D, lift it to 3D
                m3 = mask3
                if getattr(m3, "ndim", 2) == 2:
                    m3 = m3[None, ...]
                Wsic = np.where(m3, Wsic, np.nan)
                # Weighted mean thickness DOES NOT require uncertainty to be finite
                v_sit = np.isfinite(S) & np.isfinite(Wsic) & (Wsic > 0)
                wsum  = (np.where(v_sit, Wsic, 0.0)).sum(axis=(1, 2))
                num   = (np.where(v_sit, S * Wsic, 0.0)).sum(axis=(1, 2))
                SIT_wgt   = num / np.where(wsum > 0, wsum, np.nan)
                SIT_wgt_n = v_sit.sum(axis=(1, 2)).astype("int32")
                # For the uncertainty aggregates, only use pixels where uncertainty is finite
                v_u      = v_sit & np.isfinite(U)
                wsum_u   = (np.where(v_u, Wsic, 0.0)).sum(axis=(1, 2))
                num_u    = (np.where(v_u, U * Wsic, 0.0)).sum(axis=(1, 2))
                se_num_u = ((np.where(v_u, Wsic, 0.0) ** 2) * (np.where(v_u, U, 0.0) ** 2)).sum(axis=(1, 2))
                SIT_wgt_u   = num_u / np.where(wsum_u > 0, wsum_u, np.nan)
                SIT_wgt_unc = np.sqrt(se_num_u) / np.where(wsum_u > 0, wsum_u, np.nan)
            # Assemble per-file dataset (time may be >1)
            T = S.shape[0]
            ds_out = xr.Dataset(data_vars = dict(SIT         = (("time","y","x"), S),
                                                 SIT_unc     = (("time","y","x"), U),
                                                 SIT_hem     = (("time",), SIT_hem.astype("float32")),
                                                 SIT_hem_unc = (("time",), SIT_hem_unc.astype("float32")),
                                                 SIT_hem_u   = (("time",), SIT_hem_u.astype("float32")),
                                                 SIT_hem_n   = (("time",), SIT_hem_n),
                                                 SIT_wgt     = (("time",), SIT_wgt.astype("float32")),
                                                 SIT_wgt_unc = (("time",), SIT_wgt_unc.astype("float32")),
                                                 SIT_wgt_u   = (("time",), SIT_wgt_u.astype("float32")),
                                                 SIT_wgt_n   = (("time",), SIT_wgt_n)),
                                coords = dict(time = ("time", pd.to_datetime(tvals).tz_localize(None).to_numpy()),
                                              y    = ("y", np.arange(lat2d.sizes["yc"], dtype=int)),
                                              x    = ("x", np.arange(lat2d.sizes["xc"], dtype=int)),
                                              lon  = (("y","x"), lon2d.values.astype("float64")),
                                              lat  = (("y","x"), lat2d.values.astype("float64"))))
            if SIC is not None:
                ds_out["SIC"] = (("time","y","x"), SICv)
                ds_out["SIC"].attrs.update(dict(standard_name="sea_ice_area_fraction", units="percent"))
            if include_snow and SD is not None:
                ds_out["snow_depth"] = (("time","y","x"), SDv)
                ds_out["snow_depth"].attrs.update(dict(standard_name="surface_snow_thickness_where_sea_ice", units="m"))
            if include_snow and SDU is not None:
                ds_out["snow_depth_unc"] = (("time","y","x"), SDUv)
                ds_out["snow_depth_unc"].attrs.update(dict(standard_name="surface_snow_thickness_where_sea_ice standard_error", units="m"))
            if include_flags and QF is not None:
                ds_out["quality_flag"] = (("time","y","x"), qf3.astype("int8"))
            if include_flags and SF is not None:
                ds_out["status_flag"] = (("time","y","x"), sf3.astype("int8"))
            if include_flags and RC is not None:
                ds_out["region_code"] = (("time","y","x"), rc3.astype("int8"))
            out_list.append(ds_out)
            ds.close()
        if not out_list:
            # Empty shell with lat/lon if we managed to see one file’s coords
            return xr.Dataset(coords = dict(time = ("time", np.array([], dtype="datetime64[ns]")),
                                            y    = ("y", np.arange(lat_template.shape[0], dtype=int) if lat_template is not None else np.array([], dtype=int)),
                                            x    = ("x", np.arange(lat_template.shape[1], dtype=int) if lat_template is not None else np.array([], dtype=int)),
                                            lon  = (("y","x"), lon_template if lon_template is not None else np.zeros((0,0))),
                                            lat  = (("y","x"), lat_template if lat_template is not None else np.zeros((0,0))),))
        ds_all = xr.concat(out_list, dim="time").sortby("time")
        ds_all["SIT"].attrs.update(dict(standard_name="sea_ice_thickness", units="m"))
        ds_all["SIT_unc"].attrs.update(dict(standard_name="sea_ice_thickness standard_error", units="m"))
        ds_all.attrs.update(dict(note="Assembled L3 (ESA L3C / AWI l3cp_release) monthly grids with optional QC & SIC threshold; hemispheric means use cos(lat) weights."))
        return ds_all

    def make_monthly_gridded_SIT_L3(self,
                                    dt0_str            : str,
                                    dtN_str            : str,
                                    hemisphere         : str | None = None,
                                    institutions       : Iterable[str] | None = ("ESA","AWI"),
                                    sensors            : Iterable[str] | None = None,
                                    levels             : Iterable[str] | None = ("L3C","l3cp_release"),
                                    versions           : Iterable[str] | None = None,
                                    root_esa           : Path | None = None,
                                    root_awi           : Path | None = None,
                                    prefer             : str = "last",
                                    sic_thresh_percent : float | None = 15.0,
                                    mask_strategy      : str = "none",         # {"none","quality","status","both"}
                                    include_snow       : bool = True,
                                    include_flags      : bool = True,
                                    min_obs            : int = 1,
                                    time_midpoint      : bool = True,) -> xr.Dataset:
        """
        Create a monthly gridded SIT dataset from Level-3 products (ESA L3C / AWI l3cp_release).

        This convenience wrapper performs:
        1) File discovery via `build_sea_ice_satellite_paths(...)`.
        2) Indexing and de-duplication to one file per (year, month).
        3) Assembly into a single (time, y, x) dataset using `l3c_paths_to_monthly_grid(...)`.

        Parameters
        ----------
        dt0_str, dtN_str : str
            Start/end of the time window (inclusive).
        hemisphere : {"sh","nh"} or None, optional
            Hemisphere selector; defaults to `self.hemisphere_dict['abbreviation']` when present.
        institutions, sensors, levels, versions : iterable of str or None, optional
            Filters passed to `build_sea_ice_satellite_paths(...)`. By default, `levels`
            targets monthly L3 products (e.g., ("L3C","l3cp_release")).
        root_esa, root_awi : pathlib.Path or None, optional
            Root directories. If None, falls back to configured defaults (module globals or class config).
        prefer : {"first","last"}, default "last"
            If multiple files map to a month, choose first/last after deterministic sorting.
        sic_thresh_percent : float or None, default 15.0
            If set, mask cells with SIC below this threshold (percent).
        mask_strategy : {"none","quality","status","both"}, default "none"
            Optional QC masking using quality_flag and/or status_flag when present.
        include_snow : bool, default True
            Include snow depth fields when present.
        include_flags : bool, default True
            Include flag fields (quality/status/region code) when present.
        min_obs : int, default 1
            If a per-cell observation count exists, require count >= min_obs.
        time_midpoint : bool, default True
            If time bounds are present, set time coordinate to the midpoint.

        Returns
        -------
        xarray.Dataset
            Monthly gridded dataset with dimensions (time, y, x), including:
            - SIT, SIT_unc, optional SIC/snow/flags
            - SIT_hem* aggregates and optional SIC-weighted SIT aggregates (SIT_wgt*)

        See Also
        --------
        build_sea_ice_satellite_paths
        index_satellite_paths
        dedup_monthly_one_file_per_month
        l3c_paths_to_monthly_grid
        """
        # default roots from module globals if not specified
        if root_esa is None:
            root_esa = globals().get("ROOT_ESA_DEF", None)
        if root_awi is None:
            root_awi = globals().get("ROOT_AWI_DEF", None)
        paths = self.build_sea_ice_satellite_paths(dt0_str      = dt0_str,
                                                   dtN_str      = dtN_str,
                                                   hemisphere   = hemisphere,
                                                   institutions = institutions,
                                                   sensors      = sensors,
                                                   levels       = levels,
                                                   versions     = versions,
                                                   root_esa     = root_esa,
                                                   root_awi     = root_awi)
        # QA & dedup to one per month
        df    = self.index_satellite_paths(paths)
        df_mo = self.dedup_monthly_one_file_per_month(df, prefer=prefer)
        P_mos = df_mo["path"].tolist()
        return self.l3c_paths_to_monthly_grid(paths              = P_mos,
                                              hem                = hemisphere,
                                              sic_thresh_percent = sic_thresh_percent,
                                              mask_strategy      = mask_strategy,
                                              include_snow       = include_snow,
                                              include_flags      = include_flags,
                                              min_obs            = min_obs,
                                              time_midpoint      = time_midpoint)
