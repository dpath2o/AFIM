# import json, os, shutil, sys, time, logging, zarr, re
# import xarray as xr
# import pandas as pd
# import numpy  as np
# import xesmf  as xe
# from dask.distributed           import Client, LocalCluster
# from collections                import defaultdict
# from pathlib                    import Path
# from datetime                   import datetime, timedelta
# _dask_client = None
import json, os, shutil, logging
import xarray    as xr
import pandas    as pd
import numpy     as np
import xesmf     as xe
from collections import defaultdict
from pathlib     import Path
from datetime    import datetime, timedelta

class SeaIceModels:
    def __init__(self, sim_name=None, **kwargs):
        return
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
                                   overwrite       = None,
                                   delete_original = None):
        sim_name = sim_name if sim_name is not None else self.sim_name
        dt0_str  = dt0_str  if dt0_str  is not None else self.dt0_str
        dtN_str  = dtN_str  or self.dtN_str
        D_iceh   = D_iceh   or self.D_iceh
        overwrite = overwrite if overwrite is not None else self.overwrite_zarr
        delete_nc = delete_original if delete_original is not None else self.delete_original_cice_iceh_nc
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
                if delete_nc:
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
                if delete_nc:
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
        DS_grouped = defaultdict(lambda: {k: [] for k in self.valid_ice_types})
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
        valid_DS_list = valid_zarr_DS_list if valid_zarr_DS_list is not None else self.valid_ice_types
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
            CICE_SO = CICE_reM.isel(nj=self.hemisphere_dict['nj_slice'])
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
            CICE_SO = CICE_month.isel(nj=self.hemisphere_dict['nj_slice'])
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

    def load_processed_cice(self,
                            sim_name    = None,
                            rolling     = False,
                            ispd_thresh = None,
                            ice_type    = None,
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
                assert it in self.valid_ice_types, f"Invalid ice_type: {it}"
        else:
            assert ice_type in self.valid_ice_types, f"Invalid ice_type: {ice_type}"
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
                CICE = CICE.isel(nj=self.hemisphere_dict['nj_slice'])
        else:
            CICE = None
        return DS_FI, CICE
