import json, os, shutil, sys, time, logging, zarr, re
import xarray as xr
import pandas as pd
import numpy  as np
import xesmf  as xe
sys.path.insert(0, '/home/581/da1339/AFIM/src/AFIM/src')
from sea_ice_science import SeaIceScience
#from grounded_iceberg_processor import GroundedIcebergProcessor
from dask.distributed           import Client, LocalCluster
from collections                import defaultdict
from pathlib                    import Path
from datetime                   import datetime, timedelta
_dask_client = None

#class SeaIceProcessor(SeaIceScience):
class SeaIceProcessor(SeaIceScience):
    def __init__(self, sim_name, **kwargs):
        super().__init__(sim_name, **kwargs)
  

    ############################################################################################
    #                                   NSIDC GO2202_V4
    #-------------------------------------------------------------------------------------------
    def load_local_NSIDC(self, dt0_str=None, dtN_str=None, local_directory=None, monthly_files=False):
        """
        Load NSIDC sea ice concentration data from local NetCDF files over a date range.

        INPUTS:
           dt0_str         : str, optional; start date in 'YYYY-MM-DD' format. Defaults to self.dt0_str.
           dtN_str         : str, optional; end date in 'YYYY-MM-DD' format. Defaults to self.dtN_str.
           local_directory : str or Path, optional; directory containing the NSIDC NetCDF files.
                             If None, uses self.NSIDC_dict["D_original"].
           monthly_files   : bool, optional; iff True, loads monthly files; otherwise loads daily files.

        OUTPUTS
           xarray.Dataset; combined dataset loaded with xarray.open_mfdataset.
        """
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
                                     hem  = self.hemisphere_abbreviation.lower(),
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
        from pyproj import CRS, Transformer
        """
        Convert NSIDC dataset Cartesian coordinates (xgrid, ygrid) into geographic coordinates (longitude, latitude).

        This method transforms the Cartesian coordinates used in NSIDC datasets into
        standard geographic coordinates (longitude, latitude) based on the Polar
        Stereographic projection.

        INPUTS
            ds (xarray.Dataset): The NSIDC dataset containing Cartesian coordinates.

        OUTPUTS
            xarray.Dataset: The dataset with added `lon` and `lat` variables representing
                            geographic coordinates.

        NOTES:
            - The conversion is done using the Proj4 string specified in `self.NSIDC_dict["projection_string"]`.
            - The transformed geographic coordinates are added to the dataset as new variables
              named `lon` and `lat`.
            - The method assumes the dataset contains `xgrid` and `ygrid` variables, which
              represent the Cartesian coordinates in meters.
        """
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

    #-------------------------------------------------------------------------------------------
    #                            RE-ORGANISE MODEL OUTPUT DATA
    #
    # Why this dramatic size reduction makes sense:
    # Zarr uses chunked, compressed storage
    # * original NetCDF files (.nc) are not compressed. Each file is ~238 MB daily.
    # * Zarr stores data in compressed .zip-like chunks per variable, dramatically reducing file size when the data is sparse, smooth, or redundant.
    # * Many sea ice variables contain large areas of constant or near-zero values (e.g., aice, hi, frazil, taubx, etc.). These compress extremely well.
    # * Zarr avoids redundant metadata storage
    # * Each .nc file stores global attributes, dimensions, and coordinate variables repeatedly.
    # * Zarr consolidates this across time within a single file structure.
    # * intelligent chunking: by using chunk={'time': -1, 'nj': 540, 'ni': 1440} you're compressing
    #   monthly time spans while maintaining spatial chunking appropriate for hemisphere-wide analysis.
    #-------------------------------------------------------------------------------------------
    def get_cice_files_between_dates(self, D_iceh, dt0_str, dtN_str):
        dt0 = pd.to_datetime(dt0_str)
        dtN = pd.to_datetime(dtN_str)
        files = []
        for dt in pd.date_range(dt0, dtN, freq="D"):
            fname = f"iceh.{dt.strftime('%Y-%m-%d')}.nc"
            fpath = Path(D_iceh) / fname
            if fpath.exists():
                files.append(fpath)
        if not files:
            self.logger.info(f"*NO* CICE iceh.YYYY-MM-DD.nc files between the dates of {dt0_str} and {dtN_str} ... None being returned")
            return None
        return sorted(files)

    def delete_original_cice(self, P_orgs, P_iceh_zarr, m_str):
        if Path(P_iceh_zarr, ".zgroup").exists():
            self.logger.info(f"🗑️ Deleting original NetCDF files for {m_str}")
            for f in P_orgs:
                try:
                    os.remove(f)
                    self.logger.debug(f"Deleted: {f}")
                except Exception as e:
                    self.logger.warning(f" Could not delete {f}: {e}")
        else:
            self.logger.warning(f"Zarr group {P_iceh_zarr} incomplete — skipping deletion of originals")

    def daily_iceh_to_monthly_zarr(self,
                                   sim_name        = None,
                                   dt0_str         = None,
                                   dtN_str         = None,
                                   D_iceh          = None,
                                   overwrite       = False,
                                   delete_original = False):
        sim_name = sim_name or self.sim_name
        dt0_str  = dt0_str  or self.dt0_str
        dtN_str  = dtN_str  or self.dtN_str
        D_iceh   = D_iceh   or self.D_iceh
        m_grps   = defaultdict(list)
        P_orgs   = self.get_cice_files_between_dates(D_iceh, dt0_str, dtN_str)
        if not P_orgs:
            self.logger.info("No CICE files found. Noting further to do here.")
            return
        for f in P_orgs:
            dt    = datetime.strptime(f.name, "iceh.%Y-%m-%d.nc")
            m_str = dt.strftime("%Y-%m")
            m_grps[m_str].append(f)
        for m_str, P_ in m_grps.items():
            P_iceh_zarr = Path(self.D_zarr, f"iceh_{m_str}.zarr")
            if P_iceh_zarr.exists() and not overwrite:
                self.logger.info(f"Skipping existing {P_iceh_zarr}")
                if delete_original:
                    self.delete_original_cice(P_, P_iceh_zarr, m_str)
                    continue
            else:
                self.logger.info(f"Loading NetCDF files for {m_str} via xarray mfdataset ...")
                CICE_all = xr.open_mfdataset(P_,
                                             engine   = "scipy",
                                             parallel = True,
                                             combine  = "by_coords",
                                             cache    = True,
                                             chunks   = {})
                CICE_all = CICE_all.chunk({'time': -1, 'nj': 540, 'ni': 1440})
                self.logger.info(f"Subtracting one day from original dataset as CICE reports one day ahead for daily-averages")
                CICE_all["time"] = CICE_all["time"] - np.timedelta64(1, "D")
                self.logger.info(f"Writing {m_str} to {P_iceh_zarr}")
                CICE_all.to_zarr(P_iceh_zarr, mode="w", consolidated=True)
                self.get_dir_size(P_iceh_zarr)
                self.count_zarr_files(P_iceh_zarr)
                if delete_original:
                    self.delete_original_cice(P_, P_iceh_zarr, m_str)

    def monthly_zarr_iceh_time_correction(self, P_mnthly_zarr, dry_run=True):
        """
        Fix the time coordinate in a monthly iceh zarr file by shifting all time entries back by 1 day.
        This corrects for CICE's 00:00:00 datestamp which represents daily averages for the *previous* day.
        This method was done after a correction to daily_iceh_to_monthly_zarr() accounts for this timestamp
        discrepancy and hence this method may quickly become outdated/unnecessary. 

        INPUTS:
           P_mnthly_zarr : str or Path; path to the Zarr directory to correct (e.g., iceh_2011-05.zarr).
           dry_run   : bool, optional; if True, do not write any files — only print what would be changed.
        """
        from calendar import monthrange
        P_mnthly_zarr = Path(P_mnthly_zarr)
        if not P_mnthly_zarr.exists():
            self.logger.warning(f"{P_mnthly_zarr} does not exist.")
            return
        try:
            ds = xr.open_zarr(P_mnthly_zarr, consolidated=True)
        except Exception as e:
            self.logger.error(f"Failed to open {P_mnthly_zarr}: {e}")
            return
        if "time" not in ds:
            self.logger.warning(f"No 'time' coordinate found in {P_mnthly_zarr.name}")
            return
        old_time = ds["time"].values
        if len(old_time) == 0:
            self.logger.warning(f"{P_mnthly_zarr.name} has no time entries.")
            return
        # Determine expected month from filename
        try:
            m_str = P_mnthly_zarr.stem.split("_")[1]
            y, m = map(int, m_str.split("-"))
            month_start = np.datetime64(f"{y:04d}-{m:02d}-01")
            month_end = np.datetime64(f"{y:04d}-{m:02d}-{monthrange(y, m)[1]}")
        except Exception as e:
            self.logger.error(f"Could not parse month from {P_mnthly_zarr.name}: {e}")
            return
        # Check if all timestamps fall within the expected month
        time_ok = (old_time[0] >= month_start) and (old_time[-1] <= month_end)
        if time_ok and np.all((old_time >= month_start) & (old_time <= month_end)):
            self.logger.info(f"[GOOD] {P_mnthly_zarr.name} time coordinates already valid — skipping")
            return
        # Otherwise apply correction
        new_time = old_time - np.timedelta64(1, "D")
        self.logger.info(f"Fixing {P_mnthly_zarr.name}: {old_time[0]} → {new_time[0]}")
        ds["time"] = new_time
        if dry_run:
            self.logger.info(f"[dry-run] Would rewrite {P_mnthly_zarr.name} with corrected time index")
            return
        tmp_path = P_mnthly_zarr.with_suffix(".tmp.zarr")
        self.logger.info(f"Writing fixed dataset to temporary path {tmp_path}")
        ds.to_zarr(tmp_path, mode="w", consolidated=True)
        self.logger.info(f"Replacing original {P_mnthly_zarr.name}")
        shutil.rmtree(P_mnthly_zarr)
        tmp_path.rename(P_mnthly_zarr)
        self.logger.info(f"Time fix applied to {P_mnthly_zarr.name}")

    def correct_timestamp_for_all_monthly_zarr_iceh(self, sim_names=None, dry_run=True):
        """
        Loop through all monthly Zarr files in one or more simulation archives and correct the time index.
        See monthly_zarr_iceh_time_correction()

        INPUTS
           sim_name_list : list of str, optional; list of simulation names (e.g., ["FI-heavy", "PI-control"]).
                           If None, defaults to [self.sim_name].
           dry_run       : bool, optional; if True, no files will be changed — actions will be logged only.
        """
        sims = sim_names if sim_names is not None else [self.sim_name]
        for sim in sims:
            D_zarr = Path(self.D_dict['AFIM_out'],sim,'zarr')
            if not D_zarr.exists():
                self.logger.warning(f"Zarr directory not found for {sim}")
                continue
            self.logger.info(f"Scanning Zarr files in {D_zarr} ...")
            zarr_months = sorted(D_zarr.glob("iceh_????-??.zarr"))
            for P_zarr in zarr_months:
                self.monthly_zarr_iceh_time_correction(P_zarr, dry_run=dry_run)

    def get_static_grid(self, var):
        if "time" in var.dims:
            return var.isel(time=0)
        return var
    #-------------------------------------------------------------------------------------------
    #                             COMPUTE FAST ICE CLASSIFICATION
    #-------------------------------------------------------------------------------------------
    def compute_ice_speed_types(self, DS, ispd_type, temporally_average=False, mean_period=None):
        """

        Compute multiple formulations of sea ice speed and internal stress from CICE output.

        This method calculates ice speed (in m/s) and internal ice stress (in N/m²) based on vector
        components of velocity (``uvel``, ``vvel``) and internal stress (``strintx``, ``strinty``) from CICE
        model output. It supports B-grid computation (native grid), and two approaches to interpolate
        to the T-grid: a simple 4-point average (``Ta``) and an xESMF conservative remap (``Tx``).

        Optionally applies a centered temporal rolling average (e.g., 15 days) before returning.

        INPUTS:
           DS                 : xarray.Dataset; Input dataset containing CICE variables on the
                                model grid (`uvel`, `vvel`, `strintx`, `strinty`, `TLAT`, `TLON`, and `time`).
           ispd_type          : str or list of str; Ice speed types to compute. Options are:
                                + 'ispd_B': B-grid speed (√(u² + v²))
                                + 'ispd_Ta': 4-point average of B-grid speed interpolated to T-grid
                                + 'ispd_Tx': xESMF conservative interpolation to T-grid
                                + 'ispd_BT': All of the above (used for composite)
           temporally_average : bool, optional; If True, applies a centered temporal rolling mean over
                               `mean_period` days to each speed/stress field.
           mean_period        : int, optional; Temporal window size (in days) for rolling mean.
                                Defaults to `self.mean_period`.

        OUTPUTS:
           xarray.Dataset; The input dataset `DS` with new variables added for requested speed and stress types:
              + ispd_B, ists_B
              + ispd_Ta, ists_Ta
              + ispd_Tx, ists_Tx
           Each variable includes coordinate information and physical metadata.

        NOTES:
        + The B-grid formulation uses native velocity and stress components.
        + T-grid (`Ta`) is a local bilinear average across adjacent B-grid corners.
        + T-grid (`Tx`) uses regridding weights generated by xESMF for conservative interpolation.
        + Time dimension must exist; rolling averaging requires adequate time resolution.
        + Composite speed `ispd_BT` should be computed separately via `compute_composite_ice_speed`.

        See corresponding method SeaIceProcessor.compute_composite_ice_speed(), which combines the above speed
        types into a composite field.

        """
        mean_period = mean_period if mean_period is not None else self.mean_period
        if "ispd_B" in ispd_type or "ispd_Ta" in ispd_type or "ispd_BT" in ispd_type:
            self.logger.info("Computing B-grid ice speed ('ispd_B') and internal ice stress from strintx and strinty ('ists_B')")
            ispd_B = xr.apply_ufunc(np.hypot, DS['uvel']   , DS['vvel']   , dask="allowed", output_dtypes=[DS['uvel'].dtype])
            ists_B = xr.apply_ufunc(np.hypot, DS['strintx'], DS['strinty'], dask="allowed", output_dtypes=[DS['strintx'].dtype])
            if temporally_average:
                self.logger.info("Temporally-Averaging ispd_B and ists_B")
                ispd_B = ispd_B.rolling(time=mean_period, center=True, min_periods=1).mean()
                ists_B = ists_B.rolling(time=mean_period, center=True, min_periods=1).mean()
            DS['ispd_B'] = ispd_B
            DS['ists_B'] = ists_B
            if "ispd_Ta" in ispd_type or "ispd_BT" in ispd_type:
                self.logger.info("Spatially-Averaging (re-griddding) ispd_B and ists_B to T-grid: CREATING ispd_Ta and ists_Ta")
                ispd_Ta                = xr.full_like(ispd_B, fill_value=np.nan)
                ispd_stack             = xr.concat([ispd_B.isel(nj=slice(None, -1), ni=slice(None, -1)),
                                                    ispd_B.isel(nj=slice(None, -1), ni=slice(1, None)),
                                                    ispd_B.isel(nj=slice(1, None) , ni=slice(None, -1)),
                                                    ispd_B.isel(nj=slice(1, None) , ni=slice(1, None)) ], dim="offset")
                ispd_avg               = ispd_stack.mean(dim="offset", skipna=True)
                ispd_Ta[..., :-1, :-1] = ispd_avg.data
                ispd_Ta[..., :, -1]    = ispd_Ta[..., :, 0]  # basic cyclic wrap
                DS['ispd_Ta']          = xr.DataArray(ispd_Ta,
                                                    dims   = ("time", "nj", "ni"),
                                                    coords = {"time" : (("time"),DS.time.values),
                                                              "nj"   : np.arange(ispd_Ta.sizes["nj"]),
                                                              "ni"   : np.arange(ispd_Ta.sizes["ni"]),
                                                              "TLAT" : (("nj", "ni"), self.get_static_grid(DS.TLAT).values),
                                                              "TLON" : (("nj", "ni"), self.get_static_grid(DS.TLON).values)},
                                                    attrs  = {"long_name" : "T-grid interpolated B-grid ice speed",
                                                              "units"     : ispd_B.attrs.get("units", "m/s")})
                ists_Ta                = xr.full_like(ists_B, fill_value=np.nan)
                ists_stack             = xr.concat([ists_B.isel(nj=slice(None, -1), ni=slice(None, -1)),
                                                    ists_B.isel(nj=slice(None, -1), ni=slice(1, None)),
                                                    ists_B.isel(nj=slice(1, None) , ni=slice(None, -1)),
                                                    ists_B.isel(nj=slice(1, None) , ni=slice(1, None)) ], dim="offset")
                ists_avg               = ists_stack.mean(dim="offset", skipna=True)
                ists_Ta[..., :-1, :-1] = ists_avg.data
                ists_Ta[..., :, -1]    = ists_Ta[..., :, 0]  # basic cyclic wrap
                DS['ists_Ta']          = xr.DataArray(ists_Ta,
                                                    dims   = ("time", "nj", "ni"),
                                                    coords = {"time" : (("time"),DS.time.values),
                                                              "nj"   : np.arange(ists_Ta.sizes["nj"]),
                                                              "ni"   : np.arange(ists_Ta.sizes["ni"]),
                                                              "TLAT" : (("nj", "ni"), self.get_static_grid(DS.TLAT).values),
                                                              "TLON" : (("nj", "ni"), self.get_static_grid(DS.TLON).values)},
                                                    attrs  = {"long_name" : "T-grid interpolated B-grid internal ice stress",
                                                                "units"     : ists_B.attrs.get("units", "N/m^2")})
            if "ispd_Tx" in ispd_type or "ispd_BT" in ispd_type:
                self.logger.info("xESMF regrid to T-grid")
                if not self.reG_weights_defined:
                    self.define_reG_weights()
                DS_reG = self.reG_bgrid_to_tgrid_xesmf(DS)
                self.logger.info("Computing xESMF-based ice speed ('ispd_Tx') and internal ice stress from strintx and strinty ('ists_Tx')")
                ispd_Tx = xr.apply_ufunc(np.hypot, DS_reG['uvel']   , DS_reG['vvel']   , dask="allowed", output_dtypes=[DS_reG['uvel'].dtype])
                ists_Tx = xr.apply_ufunc(np.hypot, DS_reG['strintx'], DS_reG['strinty'], dask="allowed", output_dtypes=[DS_reG['strintx'].dtype])
                if temporally_average:
                    self.logger.info("Temporally-Averaging ispd_Tx and ists_Tx")
                    ispd_Tx = ispd_Tx.rolling(time=mean_period, center=True, min_periods=1).mean()
                    ists_Tx = ists_Tx.rolling(time=mean_period, center=True, min_periods=1).mean()
                DS['ispd_Tx'] = ispd_Tx.copy()
                DS['ists_Tx'] = ists_Tx.copy()
        return DS

    def compute_composite_ice_speed(self, DS):
        """

        Construct a composite sea ice speed field from multiple grid-based definitions.

        This method computes the element-wise average of three distinct ice speed fields:
        + `ispd_B`: speed computed on the native B-grid (model corners)
        + `ispd_Ta`: T-grid approximation via local 4-point B-grid average
        + `ispd_Tx`: T-grid approximation via conservative xESMF regridding

        The resulting variable, `ispd_BT`, provides a more robust and smooth estimate of ice speed
        for use in landfast ice detection and is especially useful in areas with variable grid orientation.

        INPUTS:
           DS : xarray.Dataset; Dataset containing `ispd_B`, `ispd_Ta`, and `ispd_Tx`.

        OUPUTS:
           xarray.DataArray or None; composite sea ice speed field (`ispd_BT`), or None if required
           components are missing. Metadata includes units and a descriptive long_name.

        NOTES:
        + All three speed fields must be precomputed using `compute_ice_speed_types`.
        + NaNs are ignored when averaging, enabling flexible masking near coastlines.
        + The result is useful for consistent fast ice classification across spatial domains.


        See corresponding SeaIceProcessor.compute_ice_speed_types() method, which is used to
        compute required `ispd_*` variables before combining.

        """
        ispd_types_reqd = ["ispd_B", "ispd_Ta", "ispd_Tx"]
        missing = [v for v in ispd_types_reqd if v not in DS]
        if missing:
            self.logger.warning(f"⛔ Cannot compute ispd_BT — missing variables: {missing}")
            return None
        self.logger.info("➕ Computing ispd_BT (mean of ispd_B, ispd_Ta, ispd_Tx)")
        ispd_BT      = xr.concat([DS["ispd_B"], DS["ispd_Ta"], DS["ispd_Tx"]], dim="tmp").mean(dim="tmp", skipna=True)
        ispd_BT.name = "ispd_BT"
        ispd_BT.attrs.update({"long_name"  : "Composite ice speed (average of ispd_B, ispd_Ta, ispd_Tx)",
                              "description": "Element-wise average of ispd_B, ispd_Ta, and ispd_Tx",
                              "units"      : DS["ispd_B"].attrs.get("units", "m/s")})
        return ispd_BT

    def compute_composite_internal_ice_stress(self, DS):
        """
        Function is equivalent to SeaIceProcessor.compute_composite_ice_speed() method and future revisions
        will see an integration of these two methods into a more general method that handles and component-based
        field output from CICE.
        """
        ists_types_reqd = ["ists_B", "ists_Ta", "ists_Tx"]
        missing = [v for v in ists_types_reqd if v not in DS]
        if missing:
            self.logger.warning(f"⛔ Cannot compute ists_BT — missing variables: {missing}")
            return None
        self.logger.info("➕ Computing ists_BT (mean of ists_B, ists_Ta, ists_Tx)")
        ists_BT      = xr.concat([DS["ists_B"], DS["ists_Ta"], DS["ists_Tx"]], dim="tmp").mean(dim="tmp", skipna=True)
        ists_BT.name = "ists_BT"
        ists_BT.attrs.update({"long_name"  : "Composite internal ice stress (average of ists_B, ists_Ta, ists_Tx)",
                              "description": "Element-wise average of ists_B, ists_Ta, and ists_Tx",
                              "units"      : DS["ists_B"].attrs.get("units", "N/m^2")})
        return ists_BT

    def reapply_landmask(self, DS, apply_unmodified=False):
        """

        Apply landmask to all spatial variables in the dataset.

        Uses the modified bathymetry (`kmt_mod`) field from the model grid to mask out land cells.
        Applies this mask to all variables with dimensions `("nj", "ni")`, ensuring that land areas
        are excluded from any subsequent analysis or output.

        INPUTS:
           DS : xarray.Dataset; Dataset containing sea ice fields to be masked.

        OUTPUTS:
           xarray.Dataset; Same dataset with land cells set to NaN for spatial variables.

        """
        if self.use_gi and not apply_unmodified:
            kmt_mask = xr.open_dataset(self.P_KMT_mod)['kmt'] == 1
        else:
            kmt_mask = xr.open_dataset(self.P_KMT_org)['kmt'] == 1
            #kmt_mask = self.GI_proc.G_t['kmt_mod'] == 1  # True for ocean, False for land
        for var in DS.data_vars:
            da = DS[var]
            if {"nj", "ni"}.issubset(da.dims):  # Only apply to spatial fields
                self.logger.debug(f"Masking land for variable: {var}")
                DS[var] = da.where(kmt_mask)
        self.logger.info("Applied landmask to rolled dataset")
        return DS

    def create_fast_ice_mask(self, DS, ispd_type, ispd_thresh):
        """

        Generate fast ice masks based on sea ice concentration and speed. Intentionally requires
        all three inputs. For each specified speed type (`ispd_type`), a mask is created where sea ice concentration
        exceeds `self.icon_thresh` (typically 0.9) and the ice speed is within the specified threshold range.
        This binary mask defines areas considered to be landfast ice for each speed formulation.

        INPUTES:
           DS          : xarray.Dataset; input dataset with sea ice concentration (`aice`) and ice speed fields (`ispd_*`).
           ispd_type   : list of strings; speed types to process (e.g., "ispd_B", "ispd_Ta", "ispd_Tx", "ispd_BT").
           ispd_thresh : float; threshold ice speed (in m/s) below which ice is considered "fast".

        OUTPUTS:
           dict; dictionary of masks keyed by group name (e.g., 'FI_B', 'FI_Ta', etc.).

        """
        masks = {}
        sic_mask = DS['aice'] > self.icon_thresh
        if "ispd_B" in ispd_type:
            masks['FI_B'] = sic_mask & (DS['ispd_B'] > 0) & (DS['ispd_B'] <= ispd_thresh)
        if "ispd_Ta" in ispd_type:
            masks['FI_Ta'] = sic_mask & (DS['ispd_Ta'] > 0) & (DS['ispd_Ta'] <= ispd_thresh)
        if "ispd_Tx" in ispd_type:
            masks['FI_Tx'] = sic_mask & (DS['ispd_Tx'] > 0) & (DS['ispd_Tx'] <= ispd_thresh)
        if "ispd_BT" in ispd_type:
            masks['FI_BT'] = sic_mask & (DS['ispd_BT'] > 0) & (DS['ispd_BT'] <= ispd_thresh)
        return masks

    def groupby_fast_ice_masks(self, DS, masks, m_str):
        """

        Used in both processing methods (daily and rolling) to group masked datasets by fast ice
        type for a given monthly period.

        Applies each fast ice mask (see create_fast_ice_masks) to the dataset and stores the masked result,
        along with the fast ice (`FI_mask`) and complementary pack ice (`PI_mask`) masks.
        Results are grouped in a nested dictionary by month string (`m_str`) and variable group.

        INPUTS:
           DS    : xarray.Dataset; The full input dataset for a single month.
           masks : dict; Dictionary of fast ice masks (from `create_fast_ice_mask`).
           m_str : str; Month string (e.g., "1993-07") used as a dictionary key.

        OUTPUTS:
           dict : Nested dictionary: {m_str: {FI_group: [masked Dataset]}}. Each group corresponds to a
                  different ice speed type (e.g., 'FI_B', 'FI_BT').

        """
        DS_grouped = defaultdict(lambda: {k: [] for k in self.valid_zarr_datasets})
        for group in masks:
            fi_mask          = masks[group]
            ds_fi            = DS.where(fi_mask)
            ds_fi['FI_mask'] = fi_mask
            ds_fi['PI_mask'] = ~fi_mask
            DS_grouped[m_str][group].append(ds_fi)
        return DS_grouped

    def write_to_zarr(self, DS_grouped, P_zarr_root, ispd_thresh, ispd_type, m_str, groups_to_write=None):
        """

         Write monthly fast ice datasets to disk in Zarr format, grouped by ice speed threshold type.

        This method takes preprocessed and masked datasets grouped by fast ice type (e.g., 'FI_B', 'FI_BT'), and writes
        each to a corresponding group in a consolidated Zarr store. It includes logic to optionally overwrite existing
        Zarr groups, annotate metadata, and apply consistent chunking.

        INPUTS:
           DS_grouped      : dict; nested dictionary structure {month_string: {group_name: [list of Datasets]}}
                             as returned by `groupby_fast_ice_masks`.
           P_zarr_root     : pathlib.Path; path to the root of the Zarr store for the current month
                             (e.g., `.../cice_daily_1993-07.zarr`).
           ispd_thresh     : float; ice speed threshold used for fast ice masking (added to dataset metadata).
           ispd_type       : list of str; list of ice speed types used in processing. Controls which groups are
                             eligible for writing.
           m_str           : str; month string (e.g., "1993-07") corresponding to the current batch of data.
           groups_to_write : list of str, optional; if provided, restricts Zarr writing to only those group names.
                             Otherwise, writes all available groups.

        OUTPUTS:
           None

        NOTES:
           + Fast ice groups are matched to corresponding Zarr subgroups: e.g., 'FI_B' → 'ispd_B'.
           + Existing Zarr groups are deleted if `self.overwrite_zarr_group` is True.
           + Each dataset is concatenated along `time`, sorted, and stored using a fixed chunking scheme
             (`time=-1`, `nj=540`, `ni=1440`).
           + Metadata for `ispd_thresh` and `mask_type` is stored in group-level attributes.
           + Groups named 'SO' (Southern Ocean unmasked) are treated as special and written without threshold metadata.

        """
        mask_type_map = {"FI_B"  : "fast ice mask based on thresholding ispd_B (native U-grid)",
                         "FI_Ta" : "fast ice mask based on thresholding ispd_Ta (spatially averaged ispd_B)",
                         "FI_Tx" : "fast ice mask based on thresholding ispd_Tx (xESMF regridded uvel/vvel)",
                         "FI_BT" : "fast ice mask based on thresholding ispd_BT (composite of B, Ta, Tx)"}
        for group, datasets in DS_grouped[m_str].items():
            if not datasets or (groups_to_write and group not in groups_to_write):
                continue
            if group.startswith("FI_"):
                ispd_var = group.replace("FI_", "ispd_")
                if ispd_var not in ispd_type:
                    continue
            group_dir = Path(P_zarr_root,group)
            if group_dir.exists():
                if self.overwrite_zarr_group:
                    self.logger.warning(f"Deleting existing group '{group}' at: {group_dir}")
                    try:
                        shutil.rmtree(group_dir)
                    except Exception as e:
                        self.logger.error(f"Failed to delete group '{group}': {e}")
                        continue
                else:
                    self.logger.info(f"Skipping group '{group}' (already exists and overwrite disabled)")
                    continue
            self.logger.info(f"Writing group '{group}' to Zarr store: {P_zarr_root}")
            try:
                ds_monthly = xr.concat(datasets, dim="time").sortby("time")
                if group != "SO":
                    ds_monthly.attrs["ispd_thresh"] = ispd_thresh
                    ds_monthly.attrs["mask_type"] = mask_type_map.get(group, "unknown")
                ds_monthly = ds_monthly.chunk({"time": -1, "nj": 540, "ni": 1440})
                ds_monthly.to_zarr(P_zarr_root, group=group, mode="a", consolidated=True)
            except Exception as e:
                self.logger.error(f"Failed to write group '{group}': {e}")
        DS_grouped[m_str].clear()

    def drop_unwanted_ispd_vars(self, ds, ispd_type_requested):
        keep_vars = set(ispd_type_requested) | set(ds.data_vars) - set(self.valid_ispd_types)
        drop_vars = [v for v in self.valid_ispd_types if v in ds and v not in keep_vars]
        if drop_vars:
            self.logger.info(f"Dropping unused ice speed variables: {drop_vars}")
            return ds.drop_vars(drop_vars)
        return ds

    def create_monthly_strings(self, dt0_str=None, dtN_str=None):
        dt0_str = dt0_str if dt0_str is not None else self.dt0_str
        dtN_str = dtN_str if dtN_str is not None else self.dtN_str
        dts     = pd.date_range(dt0_str, dtN_str, freq="D")
        return sorted(set(dt.strftime("%Y-%m") for dt in dts))

    def create_empty_valid_DS_dictionary(self, valid_zarr_DS_list=None):
        valid_DS_list = valid_zarr_DS_list if valid_zarr_DS_list is not None else self.valid_zarr_datasets
        return defaultdict(lambda: {k: [] for k in valid_DS_list})

    def load_iceh_zarr(self,
                       sim_name = None,
                       dt0_str  = None,
                       dtN_str  = None,
                       var_list = None):
        """

        Load a time series of monthly sea ice history datasets from Zarr archives.

        This method constructs a list of monthly date strings between the provided start and end dates,
        attempts to open each corresponding Zarr file from the simulation archive, and concatenates them
        into a continuous `xarray.Dataset` along the `time` dimension. Optionally limits loading to a
        subset of variables.

        INPUTS:
           sim_name : str, optional; simulation name (defaults to `self.sim_name`).
           dt0_str  : str, optional; start date in "YYYY-MM-DD" format (defaults to `self.dt0_str`).
           dtN_str  : str, optional; nd date in "YYYY-MM-DD" format (defaults to `self.dtN_str`).
           var_list : list of str, optional; subset of variables to load from each Zarr file.
                      If None, loads all variables.

        OUTPUTS:
           xarray.Dataset or None; concatenated dataset containing selected variables and sorted by
           time, or None if no valid data was found.

        NOTES:
           + Zarr files must follow the naming convention `iceh_YYYY-MM.zarr`.
           + Missing or corrupted Zarr files are skipped with a warning.
           + If no valid files are found, the method logs an error and returns `None`.

        """
        sim_name = sim_name or self.sim_name
        dt0_str  = dt0_str  or self.dt0_str
        dtN_str  = dtN_str  or self.dtN_str
        mo_strs  = self.create_monthly_strings( dt0_str=dt0_str, dtN_str=dtN_str )
        datasets = []
        for m_str in mo_strs:
            P_iceh = Path(self.D_zarr, f"iceh_{m_str}.zarr")
            if not P_iceh.exists():
                self.logger.warning(f"Missing monthly Zarr file: {P_iceh}")
                continue
            self.logger.info(f"Loading monthly Zarr: {P_iceh}")
            try:
                ds = xr.open_zarr(P_iceh, consolidated=True)
                if var_list is not None:
                    ds = ds[var_list]
                datasets.append(ds)
            except Exception as e:
                self.logger.error(f"Failed to load {P_iceh}: {e}")
        if not datasets:
            self.logger.error("No valid Zarr datasets found.")
            return None
        ds_all = xr.concat(datasets, dim="time").sortby("time")
        return ds_all

    def process_daily_cice(self,
                           sim_name             = None,
                           dt0_str              = None,
                           dtN_str              = None,
                           ispd_thresh          = None,
                           ispd_type            = None,
                           overwrite_zarr_group = False):
        """

        Process daily-averaged CICE model output to compute fast ice masks for each day of each month.

        This method loads daily Zarr files for a given simulation and time range, computes the specified types of
        sea ice speed fields, applies the fast ice threshold mask, groups the results by fast ice type, and saves
        the masked outputs to disk as monthly Zarr files.

        INPUTS:
           sim_name             : str, optional; Name of the simulation. Defaults to `self.sim_name`.
           dt0_str              : str, optional; Start date in "YYYY-MM-DD" format. Defaults to `self.dt0_str`.
           dtN_str              : str, optional; End date in "YYYY-MM-DD" format. Defaults to `self.dtN_str`.
           ispd_thresh          : float, optional; Ice speed threshold (in m/s) used to define landfast ice.
                                  If None, uses `self.ispd_thresh`.
           ispd_type            : str or list of str, optional; One or more speed types to use
                                  ("ispd_B", "ispd_Ta", "ispd_Tx", "ispd_BT"). If None, defaults
                                  to all valid types.
           overwrite_zarr_group : bool, optional; If True, overwrite any existing Zarr groups with the same name.
                                  Default is False

        OUTPUTS:
           xarray.Dataset or None : A time-sorted xarray.Dataset of daily fast ice concentration fields for the
                                    requested types, or None if no valid data was processed.

        """
        sim_name                  = sim_name    or self.sim_name
        dt0_str                   = dt0_str     or self.dt0_str
        dtN_str                   = dtN_str     or self.dtN_str
        ispd_type_req             = ispd_type   or self.valid_ispd_types
        ispd_thresh               = float(ispd_thresh) if ispd_thresh is not None else self.ispd_thresh
        ispd_thresh_str           = f"{ispd_thresh:.1e}".replace("e-0", "e-")
        datasets_to_return        = []
        self.overwrite_zarr_group = overwrite_zarr_group
        if isinstance(ispd_type_req, str):
            ispd_type_req = [ispd_type_req]
        assertion_err = f"Invalid requested sea ice speed 'type': {ispd_type_req}. Must be one or more of these valid types: {self.valid_ispd_types}"
        assert all(t in set(self.valid_ispd_types) for t in ispd_type_req), assertion_err
        m_DS     = self.create_empty_valid_DS_dictionary()
        CICE_all = self.load_iceh_zarr(sim_name=sim_name, dt0_str=dt0_str, dtN_str=dtN_str, var_list=self.cice_var_list)
        if CICE_all is None:
            self.logger.error("No valid CICE Zarr data to process.")
            return
        for m_str in sorted(set(pd.to_datetime(CICE_all.time.values).strftime("%Y-%m"))):
            P_zarr              = Path(self.D_zarr, f"ispd_thresh_{ispd_thresh_str}", f"cice_daily_{m_str}.zarr")
            CICE_month          = CICE_all.sel(time=CICE_all.time.dt.strftime("%Y-%m") == m_str)
            CICE_ispd           = self.compute_ice_speed_types(CICE_month, ispd_type_req)
            CICE_reM            = self.reapply_landmask(CICE_ispd)
            CICE_reM['ispd_BT'] = self.compute_composite_ice_speed(CICE_reM)
            CICE_reM['ists_BT'] = self.compute_composite_internal_ice_stress(CICE_reM)
            CICE_reM            = self.drop_unwanted_ispd_vars(CICE_reM, ispd_type_req)
            self.logger.info("Subsetting Ocean into either southern or northern hemisphere (default: southern)")
            CICE_SO = CICE_reM.isel(nj=self.hemisphere_nj_slice)
            self.logger.info("Create fast ice masks")
            masks = self.create_fast_ice_mask(CICE_SO, ispd_type_req, ispd_thresh) #
            self.logger.info("Apply fast ice masks to dataset")
            CICE_grouped = self.groupby_fast_ice_masks(CICE_SO, masks, m_str)
            for key in CICE_grouped[m_str]:
                m_DS[m_str][key].extend(CICE_grouped[m_str][key])
            fast_group = [f"FI_{t.split('_')[-1]}" for t in ispd_type_req]
            for group in fast_group:
                if group in CICE_grouped[m_str]:
                    datasets_to_return.extend(CICE_grouped[m_str][group])
            self.write_to_zarr(m_DS, P_zarr, ispd_thresh, ispd_type_req, m_str, groups_to_write=fast_group) #
            m_DS[m_str].clear()
        if datasets_to_return:
            return xr.concat(datasets_to_return, dim="time").sortby("time")
        else:
            self.logger.warning("⚠️ No fast ice datasets to return.")
            return None

    def process_rolling_cice(self,
                             sim_name             = None,
                             dt0_str              = None,
                             dtN_str              = None,
                             mean_period          = None,
                             ispd_thresh          = None,
                             ispd_type            = None,
                             overwrite_zarr_group = False):
        """

        Compute fast ice masks using a rolling mean on sea ice fields.

        This method first performs a centered temporal rolling average over `mean_period` days on the
        relevant sea ice fields, then applies threshold-based fast ice detection. Outputs are grouped
        by month and fast ice type, and saved to monthly Zarr files.

        INPUTS:
           sim_name             : str, optional; Name of the simulation. Defaults to `self.sim_name`.
           dt0_str              : str, optional; Start date in "YYYY-MM-DD" format. Defaults to `self.dt0_str`.
           dtN_str              : str, optional; End date in "YYYY-MM-DD" format. Defaults to `self.dtN_str`.
           mean_period          : int, optional; Temporal averaging window in days (e.g., 15). Defaults to `self.mean_period`.
           ispd_thresh          : float, optional; Ice speed threshold (in m/s) to define fast ice. Defaults to `self.ispd_thresh`.
           ispd_type            : str or list of str, optional; One or more speed types to use ("ispd_B", "ispd_Ta", "ispd_Tx", "ispd_BT").
                                  If None, defaults to all valid types.
           overwrite_zarr_group : bool, optional; If True, overwrite any existing Zarr groups with the same name.

        OUTPUTS:
           xarray.Dataset or None : A time-sorted xarray.Dataset of daily fast ice concentration fields for the
                                    requested types, or None if no valid data was processed.

        """
        sim_name     = sim_name    or self.sim_name
        dt0_str      = dt0_str     or self.dt0_str
        dtN_str      = dtN_str     or self.dtN_str
        mean_period  = mean_period or self.mean_period
        ispd_type    = ispd_type   or self.valid_ispd_types
        ispd_thresh  = float(ispd_thresh) if ispd_thresh is not None else self.ispd_thresh
        ispd_thresh_str = f"{ispd_thresh:.1e}".replace("e-0", "e-")
        self.overwrite_zarr_group = overwrite_zarr_group
        if isinstance(ispd_type, str):
            ispd_type = [ispd_type]
        valid_types = set(self.valid_ispd_types)
        assert all(t in valid_types for t in ispd_type), f"❌ Invalid ispd_type: {ispd_type}. Must be one or more of {self.valid_ispd_types}"
        self.logger.info(f"Rolling mean first, then fast ice masking between {dt0_str} and {dtN_str}")
        self.logger.info("Loading monthly Zarr datasets")
        CICE_all = self.load_iceh_zarr(sim_name=sim_name, dt0_str=dt0_str, dtN_str=dtN_str, var_list=self.cice_var_list)
        if CICE_all is None:
            self.logger.error("❌ No valid CICE Zarr data to process.")
            return
        datasets_to_return = []
        CICE_all = CICE_all.chunk({'time': mean_period * 2})
        CICE_ispd = self.compute_ice_speed_types(CICE_all, ispd_type, temporally_average=True, mean_period=mean_period)
        self.logger.info("Computing selective rolling means")
        CICE_roll_vars = {}
        for var in CICE_ispd.data_vars:
            da = CICE_ispd[var]
            if var.endswith(("_B", "_Ta", "_Tx")):
                self.logger.debug(f"⏭️ Skipping rolling mean for {var} (already derived)")
                CICE_roll_vars[var] = da
                continue
            cell_meas = da.attrs.get("cell_measures", "")
            if "area: uarea" in cell_meas:
                self.logger.info(f"⏭️ Skipping temporal mean for {var} due to 'cell_measures = {cell_meas}'")
                CICE_roll_vars[var] = da
                continue
            self.logger.info(f"Rolling mean on variable: {var}")
            CICE_roll_vars[var] = da.rolling(time=mean_period, center=True, min_periods=1).mean()
        CICE_roll = xr.Dataset(CICE_roll_vars, coords=CICE_ispd.coords)
        CICE_roll['time'] = CICE_roll['time'] - np.timedelta64(1, 'D')
        CICE_roll = CICE_roll.where(~np.isnan(CICE_roll['aice']), drop=False)
        CICE_roll = CICE_roll.dropna(dim="time", how="all", subset=["aice"])
        CICE_reM = self.reapply_landmask(CICE_roll)
        CICE_reM['ispd_BT'] = self.compute_composite_ice_speed(CICE_reM)
        CICE_reM['ists_BT'] = self.compute_composite_internal_ice_stress(CICE_reM)
        self.logger.info(f"Original time steps: {CICE_all.time.size}")
        self.logger.info(f"Rolled time steps: {CICE_reM.time.size}")
        time_vals = pd.to_datetime(CICE_reM.time.values)
        mo_str_da = xr.DataArray(time_vals.strftime("%Y-%m"), coords={"time": CICE_reM.time}, dims="time")
        mo_uni = np.unique(mo_str_da.values)
        for m_str in mo_uni:
            P_FI_zarr = Path(self.D_zarr, f"ispd_thresh_{ispd_thresh_str}", f"cice_rolling_{m_str}.zarr")
            CICE_month = CICE_reM.sel(time=mo_str_da == m_str)
            if CICE_month.time.size == 0:
                self.logger.warning(f"⚠️ No data for month: {m_str}")
                continue
            CICE_SO = CICE_month.isel(nj=self.hemisphere_nj_slice)
            self.logger.info(f"Rolling monthly group: {m_str} with {CICE_SO.time.size} time steps")
            masks = self.create_fast_ice_mask(CICE_SO, ispd_type, ispd_thresh)
            CICE_grouped = self.groupby_fast_ice_masks(CICE_SO, masks, m_str)
            fast_group = [f"FI_{t.split('_')[-1]}" for t in ispd_type]
            for group in fast_group:
                if group in CICE_grouped[m_str]:
                    self.logger.debug(f"Adding group '{group}' from {m_str} to return list")
                    datasets_to_return.extend(CICE_grouped[m_str][group])
            self.write_to_zarr(CICE_grouped, P_FI_zarr, ispd_thresh, ispd_type, m_str, groups_to_write=fast_group)
        if datasets_to_return:
            return xr.concat(datasets_to_return, dim="time").sortby("time")
        else:
            self.logger.warning("⚠️ No fast ice datasets to return.")
            return None

    #-------------------------------------------------------------------------------------------
    #                             SEA ICE METRIC PREPARATION
    #-------------------------------------------------------------------------------------------
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

    def load_processed_cice(self,
                            sim_name    = None,
                            rolling     = False,
                            ispd_thresh = None,
                            ice_type    = "FI_B",
                            dt0_str     = None,
                            dtN_str     = None,
                            D_zarr      = None,
                            zarr_CICE   = False,
                            chunks      = None,
                            slice_hem   = False):
        """

        Load previously processed fast ice Zarr datasets, and optionally load the original CICE iceh Zarr files.

        This utility is used to retrieve Zarr fast ice datasets from disk based on simulation metadata
        and processing flags (e.g., daily vs. rolling). It can also return the original CICE dataset from
        monthly `iceh_*.zarr` files to allow comparison, remasking, or pack ice reconstruction.

        INPUTS:
           sim_name    : str, optional; simulation name (defaults to `self.sim_name`).
           rolling     : bool, optional; if True, load Zarr files from the rolling-mean directory (`cice_rolling_*.zarr`).
                         Otherwise, load daily outputs (`cice_daily_*.zarr`).
           ispd_thresh : float, optional; threshold value used during fast ice processing, used to construct Zarr path.
           ice_type    : str or list of str; fast ice variable group(s) to load (e.g., 'FI_B', 'FI_BT').
           dt0_str     : str, optional; start date (YYYY-MM-DD).
           dtN_str     : str, optional; end date (YYYY-MM-DD).
           D_zarr      : Path, optional; root directory for the Zarr store. Defaults to the configured output path.
           zarr_CICE   : bool, optional; if True, also load the original CICE `iceh_*.zarr` files over the date range.
           chunks      : dict, optional; dictionary for Dask chunking to apply to the loaded datasets.
           slice_hem   : bool, optional; f True, spatially restricts loaded `iceh` dataset to the configured hemisphere slice.

        OUTPUTS:
           tuple of (xarray.Dataset, xarray.Dataset or None); a tuple: (fast ice dataset, original `iceh` dataset
           or None if `zarr_CICE=False`).

        NOTES
           + `ice_type` may be a single group or list of fast ice groups (e.g., 'FI_BT').
           + Skips Zarr files that do not contain the requested group or are unreadable.
           + Dates are matched to files using a `*_YYYY-MM.zarr` filename convention.

        """
        sim_name        = sim_name    or self.sim_name
        ispd_thresh     = ispd_thresh or self.ispd_thresh
        dt0_str         = dt0_str     or self.dt0_str
        dtN_str         = dtN_str     or self.dtN_str
        ispd_thresh_str = f"{ispd_thresh:.1e}".replace("e-0", "e-")
        D_zarr          = D_zarr or Path(self.config['D_dict']['AFIM_out'],sim_name,"zarr")
        if isinstance(ice_type, str) and "," in ice_type:
            ice_type = ice_type.split(",")
        if isinstance(ice_type, list):
            for it in ice_type:
                assert it in self.valid_zarr_datasets, f"Invalid ice_type: {it}"
        else:
            assert ice_type in self.valid_zarr_datasets, f"Invalid ice_type: {ice_type}"
        F_      = "cice_rolling*.zarr" if rolling else "cice_daily*.zarr"
        P_zarrs = sorted(Path(D_zarr,f"ispd_thresh_{ispd_thresh_str}").glob(F_))
        if dt0_str and dtN_str:
            dt0 = datetime.strptime(dt0_str, "%Y-%m-%d")
            dtN = datetime.strptime(dtN_str, "%Y-%m-%d")
            def file_in_range(path):
                try:
                    dt_file = datetime.strptime(path.name.split("_")[-1].split(".zarr")[0], "%Y-%m")
                    return dt0 <= dt_file <= dtN
                except Exception:
                    return False
            P_zarrs = [p for p in P_zarrs if file_in_range(p)]
        self.logger.info(f"📁 Found {len(P_zarrs)} zarr files: {[p.name for p in P_zarrs]}")
        if not P_zarrs:
            self.logger.warning(f"⚠️ No Zarr datasets found in {D_zarr}")
            return None
        DS_list = []
        for P_zarr in P_zarrs:
            self.logger.debug(f"attempting to load: {P_zarr}")
            try:
                ds = xr.open_zarr(P_zarr, group=ice_type, consolidated=True)
                DS_list.append(ds)
            except (OSError, KeyError) as e:
                self.logger.warning(f"⚠️ Skipping {P_zarr} ({ice_type}): {e}")
        if not DS_list:
            self.logger.warning(f"⚠️ No {ice_type} datasets found in any Zarr group")
            return None
        DS_FI = xr.concat(DS_list, dim="time")
        DS_FI = DS_FI.chunk(chunks)
        self.logger.info(f"✅ Loaded {ice_type}: {len(DS_FI.time)} time steps from {len(DS_list)} files")
        if zarr_CICE:
            self.logger.info(f"📦 Load monthly iceh_*.zarr files between {dt0_str} and {dtN_str}")
            P_monthly_zarrs = []
            dt0 = datetime.strptime(dt0_str, "%Y-%m-%d")
            dtN = datetime.strptime(dtN_str, "%Y-%m-%d")
            for m in pd.date_range(dt0, dtN, freq="MS"):
                m_str = m.strftime("%Y-%m")
                P_zarr = Path(self.D_zarr, f"iceh_{m_str}.zarr")
                if P_zarr.exists():
                    P_monthly_zarrs.append(P_zarr)
            if not P_monthly_zarrs:
                raise FileNotFoundError(f"No Zarr files found between {dt0_str} and {dtN_str}")
            self.logger.info(f"📁 Found {len(P_monthly_zarrs)} zarr files: {[p.name for p in P_monthly_zarrs]}")
            CICE = xr.open_mfdataset(P_monthly_zarrs,
                                         engine     = "zarr",
                                         concat_dim = "time",
                                         combine    = "nested",
                                         parallel   = True,
                                         chunks     = chunks or {})
            if slice_hem:
                CICE = CICE.isel(nj=self.hemisphere_nj_slice)
        else:
            CICE = None
        return DS_FI, CICE

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
        FI_bool      = (FI_roll_mask >= min_count).persist()
        return FI_bool

    #-------------------------------------------------------------------------------------------
    #                                      SEA ICE METRICS
    #-------------------------------------------------------------------------------------------
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
        from sea_ice_plotter import SeaIcePlotter
        sim_name        = sim_name    or self.sim_name
        ispd_thresh     = ispd_thresh or self.ispd_thresh
        ispd_thresh_str = f"{ispd_thresh:.1e}".replace("e-0", "e-")
        D_out           = D_out if D_out is not None else Path(self.D_zarr, f"ispd_thresh_{ispd_thresh_str}", "metrics")
        D_out.mkdir(parents=True, exist_ok=True)
        SI_plot      = SeaIcePlotter(sim_name  = sim_name)
        ktens        = self.sim_config.get("Ktens", "?")
        elps         = self.sim_config.get("e_f", "?")
        GI_thin      = 1 - float(self.sim_config.get("GI_thin_fact", 0.0))
        dt_range_str = f"{self.dt0_str[:4]}-{self.dtN_str[:4]}"
        af2020_df    = pd.read_csv(self.sea_ice_dict['P_AF2020_cli_csv'])
        obs_clim     = self.interpolate_obs_fia(af2020_df)
        FIA_comp     = {}
        ice_types    = [ice_type, f"{ice_type}_roll", f"{ice_type}_bool"]
        for i_type in ice_types:
            P_METS = Path(D_out, f"{i_type}_mets.zarr")
            P_sum  = Path(D_out, f"{i_type}_summary.csv")
            if P_METS.exists() and not overwrite_zarr:
                self.logger.info(f"{P_METS} exists and not overwriting--loading")
                METS             = xr.open_zarr(P_METS)
                FIA_comp[i_type] = METS['FIA']
            else:
                self.logger.info(f"{P_METS} does NOT exists and/or overwriting--computing")
                DS, CICE_SO = self.load_processed_cice(ispd_thresh = ispd_thresh,
                                                       ice_type    = ice_type,
                                                       zarr_CICE   = True,
                                                       rolling     = True if i_type==f"{ice_type}_roll" else False,
                                                       slice_hem   = True)
                if i_type==f"{ice_type}_bool":
                    bool_mask          = self.boolean_fast_ice(DS['FI_mask'], dim="time", window=7, min_count=6)
                    DS_bool            = CICE_SO.where(bool_mask)
                    DS_bool["FI_mask"] = DS["FI_mask"]
                    DS                 = DS_bool
                METS = self.compute_sea_ice_metrics(DS, sim_name, i_type, ispd_thresh_str, P_METS, P_sum, obs_clim)
            FIA_comp[i_type] = METS['FIA']
            tit_str = f"{sim_name} {i_type} ispd_thresh={ispd_thresh_str}: ktens={ktens}, elps={elps}, GI-thin={GI_thin:.2f}"
            SI_plot.plot_persistence_map(METS['FIP'],
                                         tit_str  = tit_str,
                                         ispd_str = ispd_thresh_str,
                                         ice_type = i_type,
                                         sim_name = sim_name,
                                         regional = True,
                                         plot_GI  = self.use_gi,
                                         dt_range_str  = f"{self.dt0_str[:4]}-{self.dtN_str[:4]}",
                                         overwrite_png = overwrite_png)
            SI_plot.plot_persistence_map(METS['FIP'],
                                         tit_str  = tit_str,
                                         ispd_str = ispd_thresh_str,
                                         ice_type = i_type,
                                         sim_name = sim_name,
                                         regional = False,
                                         plot_GI  = self.use_gi,
                                         dt_range_str  = f"{self.dt0_str[:4]}-{self.dtN_str[:4]}",
                                         overwrite_png = overwrite_png)
        tit_str = f"{sim_name} ispd_thresh={ispd_thresh_str}: ktens={ktens}, elps={elps}, GI-thin={GI_thin:.2f}"
        P_png   = Path(SI_plot.D_graph, "timeseries", f"FIA_{sim_name}_{ispd_thresh_str}_smoothed_{dt_range_str}.png")
        SI_plot.plot_ice_area(FIA_comp, tit_str=tit_str, P_png=P_png, roll_days=smooth_FIA_days, obs_clim=obs_clim, keys_to_plot=ice_types)

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
            from grounded_iceberg_processor import GroundedIcebergProcessor
            GI_proc       = GroundedIcebergProcessor(sim_name=self.sim_name)
            GI_total_area = GI_proc.compute_grounded_iceberg_area()
        else:
            GI_total_area = 0
        self.logger.info(f"{GI_total_area:0.2f} m^2 total circumpolar grounded iceberg area for {self.sim_name}")
        ice_area_scale = ice_area_scale if ice_area_scale is not None else self.FIC_scale
        self.logger.info(f"🧮 Spatially-integrating the product of sea ice concentrations and grid cell areas")
        IA = ((SIC * GC_area).sum(dim=("nj", "ni"))).persist()
        IA = IA + GI_total_area
        IA = IA/ice_area_scale
        return IA

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
        from scipy.ndimage import binary_dilation, distance_transform_edt
        if self.use_gi:
            P_kmt = self.P_KMT_mod
        else:
            P_kmt = self.P_KMT_org
        kmt             = xr.open_dataset(P_kmt)['kmt']
        kmt             = kmt.isel(nj=self.hemisphere_nj_slice).values
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

    #-------------------------------------------------------------------------------------------
    #                                 OBSERVATIONAL DATA
    #-------------------------------------------------------------------------------------------
    def load_AF2020_FI_area_CSV(self, doy_start):
        csv_path    = self.sea_ice_dict['P_AF2020_cli_csv']
        df          = pd.read_csv(csv_path)
        closest_doy = min(self.doy_vals, key=lambda x: abs(x - doy_start))
        row         = df[df['DOY_start'] == closest_doy].copy()
        row         = row.rename(columns={'Circumpolar': 'circumpolar',
                                          'IOsector'   : 'IOsector',
                                          'WPOsector'  : 'WPOsector',
                                          'RSsector'   : 'RSsector',
                                          'BASsector'  : 'BASsector',
                                          'WSsector'   : 'WSsector'})
        sectors     = ['Circumpolar', 'IOsector', 'WPOsector', 'RSsector', 'BASsector', 'WSsector']
        self.logger.info(f"loaded: {csv_path}")
        return row[sectors].reset_index(drop=True)

    def interpolate_obs_fia(self, csv_df, full_doy=np.arange(1, 366)):
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

    def define_AF2020_reG_weights(self, FI_obs_da):
        from pyproj import CRS, Transformer
        self.logger.info("define model grid")
        G_t = self.define_cice_grid( grid_type='t' , mask=False , build_grid_corners=False )
        G_t.to_netcdf("/g/data/gv90/da1339/grids/CICE_Tgrid_for_AF2020db_regrid.nc")
        self.logger.info("defining AF2020 regridder weights")
        F_weights     = "/g/data/gv90/da1339/grids/weights/AF_FI_2020db_to_AOM2-0p25_Tgrid_bil_meth2.nc" #self.sea_ice_dict["AF_reG_weights"]
        weights_exist = os.path.exists(F_weights)
        self.logger.info("convert 'AF_FI_OBS_2020db' Cartesian coordinates to spherical coordindates")
        crs_obs          = CRS.from_epsg(self.sea_ice_dict["projection_FI_obs"]) #unique to observations
        crs_spherical    = CRS.from_epsg(self.sea_ice_dict["projection_wgs84"])  #spherical coordinates
        transformer      = Transformer.from_crs(crs_obs, crs_spherical, always_xy=True)
        X, Y             = np.meshgrid(FI_obs_native['x'].isel(time=0).values, FI_obs_native['y'].isel(time=0).values)
        lon_obs, lat_obs = transformer.transform(X,Y)
        da_obs           = xr.DataArray(data   = FI_obs_da,
                                        dims   = ["nj", "ni"],
                                        coords = dict(lon = (["nj", "ni"], lon_obs),
                                                      lat = (["nj", "ni"], lat_obs)))
        self.logger.info(f"Model lon: {G_t['lon'].values.min()} to {G_t['lon'].values.max()}")
        self.logger.info(f"Obs lon:   {G_obs['lon'].min()} to {G_obs['lon'].max()}")
        self.logger.info(f"{'🔁 Reusing' if weights_exist else '⚙️ Creating'} regrid weights: {F_weights}")
        self.reG_AF2020 = xe.Regridder(G_obs, G_t,
                                       method            = "bilinear", #self.sea_ice_dict["AF_reG_weights_method"],
                                       periodic          = False,
                                       ignore_degenerate = True,
                                       #extrap_method     = "nearest_s2d",
                                       reuse_weights     = weights_exist,
                                       filename          = F_weights)

    def load_AF2020_FI_org_netcdf(self, P_orgs):
        from pyproj import CRS, Transformer
        F_weights        = self.sea_ice_dict["AF_reG_weights"]
        reuse_weights    = os.path.exists(F_weights)
        self.logger.info("loading FI observations with xarray mfdataset")
        self.logger.info(f"loading these files:\n{P_orgs}")
        FI_obs           = xr.open_mfdataset(P_orgs, engine='netcdf4')
        self.logger.info("masking FI observations for values greater than 4")
        mask             = (FI_obs['Fast_Ice_Time_series'] >= 4).compute()
        #FI_OBS           = FI_obs['Fast_Ice_Time_series'].where(mask)
        FI_OBS           = xr.where( FI_obs['Fast_Ice_Time_series'] >= 4, 1.0, np.nan)
        self.logger.info("define model grid")
        G_t              = self.define_cice_grid( grid_type='t' , mask=False , build_grid_corners=False )
        self.logger.info("convert 'AF_FI_OBS_2020db' Cartesian coordinates to spherical coordindates")
        crs_obs          = CRS.from_epsg(self.sea_ice_dict["projection_FI_obs"]) #unique to observations
        crs_spherical    = CRS.from_epsg(self.sea_ice_dict["projection_wgs84"])  #spherical coordinates
        transformer      = Transformer.from_crs(crs_obs, crs_spherical, always_xy=True)
        X, Y             = np.meshgrid(FI_obs_native['x'].isel(time=0).values, FI_obs_native['y'].isel(time=0).values)
        lon_obs, lat_obs = transformer.transform(X,Y)
        da_obs           = xr.DataArray(data   = FI_OBS.isel(time=0),
                                        dims   = ["nj", "ni"],
                                        coords = dict(lon = (["nj", "ni"], lon_obs),
                                                      lat = (["nj", "ni"], lat_obs)))
        self.logger.info("*** Regridding 'AF_FI_OBS_2020db' to CICE T-grid *** ")
        self.logger.info(f"\tDefine regridder; reusing weights : {reuse_weights}")
        reG_obs          = xe.Regridder( da_obs, G_t,
                                         method            = "bilinear",
                                         periodic          = False,
                                         ignore_degenerate = True,
                                         reuse_weights     = reuse_weights,
                                         filename          = P_weights)
        self.logger.info(f"\tRegridding masked observational dataset")
        FI_OBS_reG        = reG_obs( FI_OBS )
        #FI_obs_reG_dat    = self.reG_AF2020(FI_obs["Fast_Ice_Time_series"])
        #mask              = (FI_obs_reG_dat >= 4).compute()
        #FI_obs_binary     =         xr.where(FI_obs_reG_dat )
        FI                = (('t_FI_obs', 'nj', 'ni'),
                             FI_OBS_reG.values,
                             {'long_name': FI_obs['Fast_Ice_Time_series'].attrs['long_name']})
        t_alt             = (('t_FI_obs'),
                             FI_obs.date_alt.values,
                             {'long_name': FI_obs.date_alt.attrs['long_name'],
                              'description': FI_obs.date_alt.attrs['description']})
        t_coords          = (('t_FI_obs'),
                             FI_obs.time.values,
                             {'description': "Start date of 15- or 20-day image mosaic window.",
                              'units': "days since 2000-1-1 0:0:0"})
        x_coords          = (('nj','ni'),
                             G_t['lon'].values,
                             {'long_name': 'longitude', 'units': 'degrees_east'})
        y_coords          = (('nj','ni'),
                             G_t['lat'].values,
                             {'long_name': 'latitude', 'units': 'degrees_north'})
        self.logger.info("converted AF2020 database for use with SeaIceProcessor")
        return xr.Dataset({'FI'       : FI,
                           'FI_t_alt' : t_alt },
                          coords = {'t_FI_obs' : t_coords,
                                    'lon'      : x_coords,
                                    'lat'      : y_coords})

    def filter_AF2020_FI_by_date(self, dt0_str, dtN_str):
        dt0       = datetime.strptime(dt0_str, "%Y-%m-%d")
        dtN       = datetime.strptime(dtN_str, "%Y-%m-%d")
        D_obs     = Path(self.sea_ice_dict['D_AF2020_db_org'])
        yrs_reqd  = list(range(dt0.year, dtN.year + 1))
        P_orgs    = [D_obs / f"FastIce_70_{yr}.nc" for yr in yrs_reqd]
        ds        = self.load_AF2020_FI_org_netcdf(P_orgs)
        alt_dates = pd.to_datetime(ds['FI_t_alt'].values.astype(str), format='%Y%m%d')
        ds        = ds.assign_coords(obs_date=("t_FI_obs", alt_dates))
        ds        = ds.set_index(t_FI_obs="obs_date")
        matched   = ds.sel(t_FI_obs=slice(dt0, dtN))
        if matched.dims['t_FI_obs'] == 0:
            self.logger.warning(f"No matching observational dates found between {dt0_str} and {dtN_str}")
        return matched

    def create_AF2020_FI_zarr(self):
        P_zarr = Path(self.sea_ice_dict["P_AF_2020db_avg"])
        if P_zarr.exists():
            self.logger.info(f"Averaged observational gridded climatology already exists\n{P_zarr}")
            return xr.open_zarr(P_zarr)
        elif P_zarr.exists() and overwrite:
            self.logger.info(f"Averaged observational gridded climatology exists *but* over-writing has been requested\n\tOVER-WRITING: {P_zarr}")
            shutil.rmtree(P_zarr)
        else:
            self.logger.info(f"Averaged observational gridded climatology does *NOT* exist\n\tCREATING NEW: {P_zarr}")
        D_obs  = Path(self.sea_ice_dict['D_AF2020_db_org'])
        P_orgs = sorted(D_obs.glob("FastIce_70_*.nc"))
        self.logger.info(f"loading all observational fast ice netcdf files:\n{P_orgs}")
        ds_all    = self.load_AF2020_FI_org_netcdf(P_orgs)
        alt_dates = pd.to_datetime(ds_all['FI_t_alt'].astype(str), format='%Y%m%d')
        doy_vals  = alt_dates.dayofyear
        ds_all    = ds_all.assign_coords(doy=("t_FI_obs", doy_vals))
        ds_all    = ds_all.where(ds_all.doy.isin(self.doy_vals), drop=True)
        grouped   = ds_all.groupby("doy").mean(dim="t_FI_obs", skipna=True)
        grouped   = grouped.rename_dims({"doy": "t_doy"})
        grouped   = grouped.rename_vars({"FI": "FI_OBS_GRD"})
        grouped   = grouped.assign_attrs(doy_vals=self.doy_vals).persist()
        self.logger.info(f"gridded climatology now looks like\n{grouped}")
        self.logger.info(f"writing this to {P_zarr}")
        t1 = time.time()
        grouped.to_zarr(P_zarr)
        self.logger.info(f"\tzarr written in {time.time()-t1:0.2f} seconds")
        return grouped

    def create_AF2020_daily_zarr(self, dt0_str="1993-01-01", dtN_str="2023-12-31"):
        """
        Generate a daily, regridded observational dataset of fast ice from AF2020:
        - 14-day mosaics interpolated to daily (2000–2018)
        - Climatological daily average used outside that range
        - Output written to monthly Zarr files
        """
        G_t = self.define_cice_grid( grid_type='t' , mask=False , build_grid_corners=False )
        dt0 = pd.to_datetime(dt0_str)
        dtN = pd.to_datetime(dtN_str)
        all_dates = pd.date_range(dt0, dtN, freq="D")
        years_all = np.unique(all_dates.year)
        years_obs = np.arange(2000, 2019)
        years_clim = [y for y in years_all if y not in years_obs]
        self.logger.info("Interpolating observational fast ice (2000–2018) to daily")
        ds_obs_14d = self.filter_AF2020_FI_by_date("2000-01-01", "2018-12-31")
        ds_obs_daily = ds_obs_14d["FI"].interp( t_FI_obs=pd.date_range("2000-01-01", "2018-12-31", freq="D") ).rename({"t_FI_obs": "time"})
        self.logger.info("Interpolating climatological fast ice to daily DOY")
        ds_clim = self.create_AF2020_FI_zarr()
        ds_clim_daily = ds_clim["FI_OBS_GRD"].interp( t_doy=np.arange(1, 367) )
        self.logger.info(f"Expanding climatology to years: {years_clim}")
        ds_clim_list = []
        for year in years_clim:
            dates = pd.date_range(f"{year}-01-01", f"{year}-12-31", freq="D")
            doy_vals = dates.dayofyear
            da = ds_clim_daily.sel(t_doy=xr.DataArray(doy_vals, dims="time"))
            da = da.assign_coords(time=("time", dates))
            ds_clim_list.append(da)
        ds_clim_expanded = xr.concat(ds_clim_list, dim="time") if ds_clim_list else None
        self.logger.info("Combining interpolated observations and climatology")
        if ds_clim_expanded is not None:
            ds_full = xr.concat([ds_obs_daily, ds_clim_expanded], dim="time").sortby("time")
        else:
            ds_full = ds_obs_daily
        ds_full = ds_full.sel(time=slice(dt0, dtN))
        ds_full = ds_full.expand_dims("time") if "time" not in ds_full.dims else ds_full
        ds_full = ds_full.assign_coords(TLAT = (("nj", "ni"), G_t["lat"]),
                                        TLON = (("nj", "ni"), G_t["lon"]))
        out_dir = Path(self.sea_ice_dict["D_AF2020_daily_zarr"])
        self.logger.info(f"Writing daily observational Zarrs to {out_dir}")
        for year in np.unique(ds_full.time.dt.year.values):
            for month in np.unique(ds_full.time.sel(time=ds_full.time.dt.year == year).dt.month):
                ds_month = ds_full.sel(time=ds_full.time.dt.year == year).sel(time=ds_full.time.dt.month == month)
                if ds_month.time.size == 0:
                    continue
                out_path = out_dir / f"{year}/{month:02d}.zarr"
                out_path.parent.mkdir(parents=True, exist_ok=True)
                self.logger.info(f"Writing: {out_path}")
                ds_month.to_dataset(name="FI_OBS_DAILY").to_zarr(out_path, mode="w")
        self.logger.info("Completed writing observational daily Zarr dataset.")

    def load_AF2020_daily(self, dt0_str: str, dtN_str: str) -> xr.Dataset:
        """
        Load regridded daily observational landfast sea ice from AF2020, written as monthly Zarr files (via create_AF2020_daily_zarr).

        INPUTS:
           dt0_str : str; start date (e.g., "1993-01-01")
           dtN_str : str; end date (e.g., "2023-12-31")

        OUTPUTS:
           xr.Dataset : daily observational fast ice on model T-grid.
                        Contains:
                        + FI_OBS_DAILY (time, nj, ni)
                        + TLAT, TLON (nj, ni)
        """
        dt0 = pd.to_datetime(dt0_str)
        dtN = pd.to_datetime(dtN_str)
        all_dates = pd.date_range(dt0, dtN, freq="D")
        base_dir = Path(self.sea_ice_dict["D_AF2020_daily_zarr"])
        zarr_paths = []
        for date in all_dates:
            path = base_dir / f"{date.year}/{date.month:02d}.zarr"
            if path.exists():
                zarr_paths.append(str(path))
        if not zarr_paths:
            self.logger.warning(f"No AF2020 daily Zarr files found between {dt0_str} and {dtN_str}")
            return None
        self.logger.info(f"Loading {len(zarr_paths)} monthly Zarrs for AF2020 observations")
        ds = xr.open_mfdataset(zarr_paths, engine="zarr", combine="by_coords")
        ds = ds.sel(time=slice(dt0, dtN))
        if "FI_OBS_DAILY" not in ds.data_vars:
            self.logger.error("Variable 'FI_OBS_DAILY' not found in loaded dataset.")
            return None
        return ds

