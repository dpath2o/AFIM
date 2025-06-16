import json, os, shutil, logging, re
import xarray    as xr
import pandas    as pd
import numpy     as np
import xesmf     as xe
from collections import defaultdict
from pathlib     import Path
from datetime    import datetime, timedelta
import concurrent.futures
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

class SeaIceClassification:
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
    @staticmethod
    def get_month_range(year_month):
        dt0 = datetime.strptime(year_month + "-01", "%Y-%m-%d")
        dtN = (dt0.replace(day=28) + timedelta(days=4)).replace(day=1) - timedelta(days=1)
        return [dt0 + timedelta(days=i) for i in range((dtN - dt0).days + 1)]

    @staticmethod
    def verify_month(args_tuple):
        zarr_path, nc_dir, done_marker, dry_run = args_tuple
        year_month = zarr_path.stem.split("_")[1]
        if done_marker.exists():
            return f"[SKIP] {year_month}: already verified (.done exists)"
        dt_list = SeaIceModels.get_month_range(year_month)
        nc_files = [nc_dir / f"iceh.{dt.strftime('%Y-%m-%d')}.nc" for dt in dt_list]
        existing_nc_files = [f for f in nc_files if f.exists()]
        if not existing_nc_files:
            return f"[SKIP] {year_month}: no NetCDFs remain"
        try:
            zarr = xr.open_zarr(zarr_path)
        except Exception as e:
            return f"[FAIL] {year_month}: cannot open Zarr: {e}"
        zarr_dates = pd.to_datetime(zarr["time"].values).normalize()
        expected_dates = pd.to_datetime([dt.date() for dt in dt_list])
        all_dates_present = set(expected_dates).issubset(set(zarr_dates))
        variables_expected = set()
        for f in existing_nc_files:
            try:
                ds = xr.open_dataset(f)
                variables_expected.update(ds.data_vars)
            except Exception as e:
                return f"[FAIL] {year_month}: cannot open NetCDF {f.name}: {e}"
        variables_missing = [v for v in variables_expected if v not in zarr.data_vars]
        if not all_dates_present:
            return f"[FAIL] {year_month}: missing time steps in Zarr"
        if variables_missing:
            return f"[FAIL] {year_month}: Zarr missing variables: {variables_missing}"
        if not dry_run:
            done_marker.touch()
        return f"[OK]   {year_month}: verified ({len(existing_nc_files)} files)"

    def verify_zarr_and_cleanup_netcdf(self,
                                       dry_run     = True,
                                       delete      = False,
                                       max_workers = 4):
        """
        Verify that each monthly Zarr directory under a given simulation archive
        matches all expected daily NetCDF files in date coverage and variable presence.
        Optionally delete NetCDF files after verification.
        """
        self.D_iceh_nc = Path(self.D_sim,"history","daily")
        P_clean_log    = Path(self.D_sim,"cleanup.log")
        zarr_months    = sorted([p for p in self.D_zarr.glob("iceh_????-??.zarr") if p.is_dir()])
        tasks = []
        for zarr_path in zarr_months:
            ym          = zarr_path.stem.split("_")[1]
            done_marker = self.D_zarr / f".done_{ym}"
            tasks.append((zarr_path, self.D_iceh_nc, done_marker, dry_run))
        with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
            results = list(executor.map(SeaIceModels.verify_month, tasks))
        for res in results:
            print(res)
        with open(P_clean_log, "a") as logf:
            logf.write("\n# Last updated: " + datetime.now().strftime("%Y-%m-%d %H:%M:%S") + "\n")
            for res in results:
                logf.write(res + "\n")
        if delete:
            print("\nDeletion mode active.")
            verified_months = []
            total_files = []
            log_entries = []
            for zarr_path in zarr_months:
                ym = zarr_path.stem.split("_")[1]
                done_marker = self.D_zarr / f".done_{ym}"
                if not done_marker.exists():
                    continue  # skip unverified
                nc_files = list(self.D_iceh_nc.glob(f"iceh.{ym}-??.nc"))
                if nc_files:
                    verified_months.append((ym, nc_files))
                    total_files.extend(nc_files)
            if not total_files:
                print("No deletable NetCDF files found.")
            else:
                print(f"\n🔍 {len(total_files)} NetCDF files across {len(verified_months)} verified months are eligible for deletion.")
                confirm = input("Confirm delete all these files? [y/N] ").strip().lower()
                if confirm == "y":
                    for ym, files in verified_months:
                        for f in files:
                            try:
                                f.unlink()
                                print(f"[DELETED] {f.name}")
                                log_entries.append(f"[DELETED] {f}")
                            except Exception as e:
                                print(f"[ERROR] Could not delete {f.name}: {e}")
                                log_entries.append(f"[ERROR] Failed to delete {f}: {e}")
                    log_entries.append(f"# Deletion complete: {len(total_files)} files removed")
                else:
                    print("Deletion cancelled.")
                    log_entries.append("# Deletion prompt declined — no files deleted")
            with open(P_clean_log, "a") as logf:
                for entry in log_entries:
                    logf.write(entry + "\n")

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
        overwrite = overwrite if overwrite is not None else self.overwrite_zarr_group
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
    def simple_spatial_averaging_bgrid_to_tgrid(self, var):
        """
        Dask-safe 4-point unweighted average from B-grid to T-grid.

        Uses efficient array shifting and avoids costly concatenation over new dimensions.

        Parameters
        ----------
        var : xr.DataArray
            2D or 3D (time, nj, ni) array on B-grid.

        Returns
        -------
        xr.DataArray
            Averaged field on T-grid with shape (time, nj, ni).
        """
        self.logger.info("  ↪ Slicing corner points for averaging...")
        v00          = var.isel(nj=slice(None, -1)  , ni=slice(None, -1  ))
        v01          = var.isel(nj=slice(None, -1)  , ni=slice(1   , None))
        v10          = var.isel(nj=slice(1   , None), ni=slice(None, -1  ))
        v11          = var.isel(nj=slice(1   , None), ni=slice(1   , None))
        self.logger.info("  ↪ Computing mean of four corners...")
        avg          = (v00 + v01 + v10 + v11) / 4.0
        self.logger.info("  ↪ Padding with NaNs to restore original grid size...")
        avg          = avg.pad(nj=(0,1), ni=(0,1), constant_values=np.nan)
        self.logger.info("  ↪ Applying cyclic wrap for last column...")
        avg[..., -1] = avg[..., 0]
        if "time" in var.coords:
            avg = avg.assign_coords(time=var.time)
            self.logger.info("  ↪ Time coordinate restored.")
        return avg

    def compute_ice_magnitude_from_ice_components_on_Bgrid(self, DS,
                                                           ivec_type                  = None,
                                                           ice_vector_component_names = None,
                                                           temporally_average         = False,
                                                           mean_period                = None):
        """
        Compute the magnitude from two vector components(default is sea ice speed; e.g. components are: 'uvel', 'vvel') 
        across multiple grid definitions: B, Ta, Tx, and BT. To be more clear:
        + B  : represents magnitudes that are returned on the native B-grid
        + Ta : magnitudes computed on the B-grid then a four-point spatial average is applied 
               which forms a crude (un-weighted) re-gridding
        + Tx : magnitudes computed on the B-grid then a weighted spatial average is applied using
               xesmf module
        + BT : uses the mean from the above three methods; the reason for this is that it adequately addresses 
               a scientifically sound spatially consistent sea ice velocity field that accounts for 
               the no-slip boundary condition that is a result of the B-grid momentum equations

        INPUTS
        DS : xr.Dataset; input dataset with where vector components reside on a B-grid and will either be

        ivec_type : tuple of str; this name is a legacy from the first method used to do this 
            Types of speed to compute. Options include:
            - 'B' : native B-grid speed
            - 'Ta': 4-pt local average of B-grid to T-grid
            - 'Tx': xESMF regridded speed to T-grid
            - 'BT': composite of the above three (computed last)

        temporally_average : bool
            Whether to apply a rolling temporal mean to each speed field.

        mean_period : int or None
            Length of rolling window. Uses default if None.

        Returns
        -------
        xr.Dataset
            The original dataset with additional `ispd_*` fields.
        """
        mean_period = mean_period if mean_period is not None else self.mean_period
        ivec_type   = ivec_type                  if ivec_type                  is not None else self.ivec_type
        ivec_names  = ice_vector_component_names if ice_vector_component_names is not None else ["uvel","vvel"]
        xvec        = DS[ivec_names[0]]
        yvec        = DS[ivec_names[1]]
        if set(ivec_names) == {"uvel", "vvel"}:
            mag_var_name_prefix = "ispd"
        else:
            match               = re.match(r"^(.*?)(x|y)?$", ivec_names[0])
            mag_var_name_prefix = match.group(1)
        mag_var_name_prefix = mag_var_name_prefix.rstrip("_")
        mag_var_names       = {suffix: f"{mag_var_name_prefix}_{suffix}" for suffix in self.valid_ivec_types}
        B_name              = mag_var_names['B']
        Ta_name             = mag_var_names['Ta']
        Tx_name             = mag_var_names['Tx']
        BT_name             = mag_var_names['BT']
        self.logger.info("Computing native B-grid sea ice vector magnitudes")
        B = xr.apply_ufunc(np.hypot, xvec, yvec, dask="allowed")
        if temporally_average:
            self.logger.info("Averaging B over time")
            B = B.rolling(time=mean_period, center=True, min_periods=1).mean()
        DS[B_name] = B   
        if "Ta" in ivec_type or "BT" in ivec_type:
            self.logger.info("Interpolating B-grid to T-grid using 4-pt average method: simple_spatial_averaging_bgrid_to_tgrid()")
            B = B.chunk({"time": 30, "nj": 540, "ni": 720})
            Ta = self.simple_spatial_averaging_bgrid_to_tgrid(B)
            if temporally_average:
                self.logger.info("Averaging Ta over time")
                Ta = Ta.rolling(time=mean_period, center=True, min_periods=1).mean()
            DS[Ta_name] = Ta.assign_attrs({"long_name": "T-grid interpolated B-grid ice speed",
                                           "units"    : B.attrs.get("units", "m/s")})
        if "Tx" in ivec_type or "BT" in ivec_type:
            self.logger.info("Interpolating to T-grid using xESMF")
            if not self.reG_weights_defined:
                self.define_reG_weights()
            Tx = self.reG_bgrid_to_tgrid_xesmf(B)
            if Tx is None:
                raise RuntimeError("Regridding failed — check coordinate alignment and weights.")
            Tx               = Tx.assign_coords(TLAT=DS["TLAT"], TLON=DS["TLON"])
            Tx["TLAT"].attrs = DS["TLAT"].attrs
            Tx["TLON"].attrs = DS["TLON"].attrs
            if temporally_average:
                self.logger.info("Averaging Tx over time")
                Tx = Tx.rolling(time=mean_period, center=True, min_periods=1).mean()
            DS[Tx_name] = Tx
        if "BT" in ivec_type:
            required = [k for k in [B_name, Ta_name, Tx_name] if k in DS]
            if len(required) < 3:
                self.logger.warning(f"Cannot compute BT — missing one or more of {required}")
            else:
                self.logger.info("Computing composite BT method of re-gridding")
                DS          = self.reapply_landmask(DS)
                BT          = xr.concat([DS[B_name],
                                         DS[Ta_name],
                                         DS[Tx_name]], dim="tmp").mean(dim="tmp", skipna=True)
                DS[BT_name] = BT.assign_attrs({"long_name"  : f"Composite ice speed (average of {B_name}, {Ta_name}, {Tx_name})",
                                               "description": "Element-wise average of grids B, Ta, and Tx",
                                               "units"      : B.attrs.get("units", "m/s")})
        return DS

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

    def create_fast_ice_mask(self, DS, ivec_type, ispd_thresh):
        """

        Generate fast ice masks based on sea ice concentration and speed. Intentionally requires
        all three inputs. For each specified speed type (`ivec_type`), a mask is created where sea ice concentration
        exceeds `self.icon_thresh` (typically 0.9) and the ice speed is within the specified threshold range.
        This binary mask defines areas considered to be landfast ice for each speed formulation.

        INPUTES:
           DS          : xarray.Dataset; input dataset with sea ice concentration (`aice`) and ice speed fields (`ispd_*`).
           ivec_type   : list of strings; speed types to process (e.g., "ispd_B", "ispd_Ta", "ispd_Tx", "ispd_BT").
           ispd_thresh : float; threshold ice speed (in m/s) below which ice is considered "fast".

        OUTPUTS:
           dict; dictionary of masks keyed by group name (e.g., 'FI_B', 'FI_Ta', etc.).

        """
        masks = {}
        sic_mask = DS['aice'] > self.icon_thresh
        if "B" in ivec_type:
            masks['FI_B'] = sic_mask & (DS['ispd_B'] > 0) & (DS['ispd_B'] <= ispd_thresh)
        if "Ta" in ivec_type:
            masks['FI_Ta'] = sic_mask & (DS['ispd_Ta'] > 0) & (DS['ispd_Ta'] <= ispd_thresh)
        if "Tx" in ivec_type:
            masks['FI_Tx'] = sic_mask & (DS['ispd_Tx'] > 0) & (DS['ispd_Tx'] <= ispd_thresh)
        if "BT" in ivec_type:
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

    def write_to_zarr(self, DS_grouped, P_zarr_root, ispd_thresh, ivec_type, m_str, groups_to_write=None):
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
           ivec_type       : list of str; list of ice speed types used in processing. Controls which groups are
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
            self.logger.info(f"Inspecting group: {group} with {len(datasets)} datasets")
            if not datasets or (groups_to_write and group not in groups_to_write):
                self.logger.info(f"Skipping group: {group} (empty or not in groups_to_write)")
                continue
            if group.startswith("FI_"):
                ispd_var = group.replace("FI_", "")  # e.g. 'BT'
                if ispd_var not in ivec_type:
                    self.logger.info(f"Skipping '{group}' as '{ispd_var}' not in ivec_type: {ivec_type}")
                    continue
            group_dir = Path(P_zarr_root,group)
            print(group_dir)
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

    def drop_unwanted_ispd_vars(self, ds, ivec_type_requested):
        keep_vars = set(ivec_type_requested) | set(ds.data_vars) - set(self.valid_ivec_types)
        drop_vars = [v for v in self.valid_ivec_types if v in ds and v not in keep_vars]
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
        N_loaded = 0
        self.logger.info(f"Loading monthly Zarr files: {self.D_zarr}")
        for m_str in mo_strs:
            P_iceh = Path(self.D_zarr, f"iceh_{m_str}.zarr")
            if not P_iceh.exists():
                self.logger.warning(f"Missing monthly Zarr file: {P_iceh}")
                continue
            self.logger.debug(f"Loading monthly Zarr: {P_iceh}")
            try:
                ds = xr.open_zarr(P_iceh, consolidated=True)
                if var_list is not None:
                    ds = ds[var_list]
                datasets.append(ds)
                N_loaded += 1
            except Exception as e:
                self.logger.error(f"Failed to load {P_iceh}: {e}")
        if not datasets:
            self.logger.error("No valid Zarr datasets found.")
            return None
        ds_all = xr.concat(datasets, dim="time").sortby("time")
        n_time = ds_all.sizes.get("time", 0)
        self.logger.info(f"Loaded {N_loaded}-zarr files covering {n_time} time steps from {dt0_str} to {dtN_str}")        
        return ds_all

    def process_daily_cice(self,
                           sim_name             = None,
                           dt0_str              = None,
                           dtN_str              = None,
                           ispd_thresh          = None,
                           ivec_type            = None,
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
           ivec_type            : str or list of str, optional; One or more speed types to use
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
        ivec_type_req             = ivec_type   or self.valid_ivec_types
        ispd_thresh               = float(ispd_thresh) if ispd_thresh is not None else self.ispd_thresh
        ispd_thresh_str           = f"{ispd_thresh:.1e}".replace("e-0", "e-")
        datasets_to_return        = []
        self.overwrite_zarr_group = overwrite_zarr_group
        if isinstance(ivec_type_req, str):
            ivec_type_req = [ivec_type_req]
        assertion_err = f"Invalid requested sea ice speed 'type': {ivec_type_req}. Must be one or more of these valid types: {self.valid_ivec_types}"
        assert all(t in set(self.valid_ivec_types) for t in ivec_type_req), assertion_err
        m_DS     = self.create_empty_valid_DS_dictionary()
        CICE_all = self.load_iceh_zarr(sim_name=sim_name, dt0_str=dt0_str, dtN_str=dtN_str, var_list=self.cice_var_list)
        if CICE_all is None:
            self.logger.error("No valid CICE Zarr data to process.")
            return
        for m_str in sorted(set(pd.to_datetime(CICE_all.time.values).strftime("%Y-%m"))):
            P_zarr       = Path(self.D_zarr, f"ispd_thresh_{ispd_thresh_str}", f"cice_daily_{m_str}.zarr")
            CICE_month   = CICE_all.sel(time=CICE_all.time.dt.strftime("%Y-%m") == m_str)
            CICE_ispd    = self.compute_ice_magnitude_from_ice_components_on_Bgrid(CICE_month)
            CICE_reM     = self.drop_unwanted_ispd_vars(CICE_ispd, ivec_type_req)
            self.logger.info("Subsetting Ocean into either southern or northern hemisphere (default: southern)")
            CICE_SO      = CICE_reM.isel(nj=self.hemisphere_dict['nj_slice'])
            self.logger.info("Create fast ice masks")
            masks        = self.create_fast_ice_mask(CICE_SO, ivec_type_req, ispd_thresh) #
            self.logger.info("Apply fast ice masks to dataset")
            CICE_grouped = self.groupby_fast_ice_masks(CICE_SO, masks, m_str)
            for key in CICE_grouped[m_str]:
                m_DS[m_str][key].extend(CICE_grouped[m_str][key])
            fast_group = [f"FI_{t.split('_')[-1]}" for t in ivec_type_req]
            for group in fast_group:
                if group in CICE_grouped[m_str]:
                    datasets_to_return.extend(CICE_grouped[m_str][group])
            self.write_to_zarr(m_DS, P_zarr, ispd_thresh, ivec_type_req, m_str, groups_to_write=fast_group) #
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
                             ivec_type            = None,
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
           ivec_type            : str or list of str, optional; One or more speed types to use ("ispd_B", "ispd_Ta", "ispd_Tx", "ispd_BT").
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
        ivec_type    = ivec_type   or self.valid_ivec_types
        ispd_thresh  = float(ispd_thresh) if ispd_thresh is not None else self.ispd_thresh
        ispd_thresh_str = f"{ispd_thresh:.1e}".replace("e-0", "e-")
        self.overwrite_zarr_group = overwrite_zarr_group
        if isinstance(ivec_type, str):
            ivec_type = [ivec_type]
        valid_types = set(self.valid_ivec_types)
        assert all(t in valid_types for t in ivec_type), f"❌ Invalid ivec_type: {ivec_type}. Must be one or more of {self.valid_ivec_types}"
        self.logger.info(f"Rolling mean first, then fast ice masking between {dt0_str} and {dtN_str}")
        self.logger.info("Loading monthly Zarr datasets")
        CICE_all = self.load_iceh_zarr(sim_name=sim_name, dt0_str=dt0_str, dtN_str=dtN_str, var_list=self.cice_var_list)
        if CICE_all is None:
            self.logger.error("❌ No valid CICE Zarr data to process.")
            return
        datasets_to_return = []
        CICE_all  = CICE_all.chunk({'time': mean_period * 2})
        CICE_ispd = self.compute_ice_magnitude_from_ice_components_on_Bgrid(CICE_all, ivec_type, temporally_average=True, mean_period=mean_period)
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
        CICE_roll         = xr.Dataset(CICE_roll_vars, coords=CICE_ispd.coords)
        CICE_roll['time'] = CICE_roll['time'] - np.timedelta64(1, 'D')
        CICE_roll         = CICE_roll.where(~np.isnan(CICE_roll['aice']), drop=False)
        CICE_reM          = CICE_roll.dropna(dim="time", how="all", subset=["aice"])
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
            masks = self.create_fast_ice_mask(CICE_SO, ivec_type, ispd_thresh)
            CICE_grouped = self.groupby_fast_ice_masks(CICE_SO, masks, m_str)
            fast_group = [f"FI_{t.split('_')[-1]}" for t in ivec_type]
            for group in fast_group:
                if group in CICE_grouped[m_str]:
                    self.logger.debug(f"Adding group '{group}' from {m_str} to return list")
                    datasets_to_return.extend(CICE_grouped[m_str][group])
            self.write_to_zarr(CICE_grouped, P_FI_zarr, ispd_thresh, ivec_type, m_str, groups_to_write=fast_group)
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
        ice_type        = ice_type    or self.ice_type
        dt0_str         = dt0_str     or self.dt0_str
        dtN_str         = dtN_str     or self.dtN_str
        chunks          = chunks      or {'time': 31}
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
        self.logger.info(f"Found {len(P_zarrs)} zarr files")
        self.logger.debug(f"{[p.name for p in P_zarrs]}")
        if not P_zarrs:
            self.logger.warning(f"No Zarr datasets found in {D_zarr}")
            return None
        DS_list = []
        for P_zarr in P_zarrs:
            self.logger.debug(f"attempting to load: {P_zarr}")
            try:
                ds = xr.open_zarr(P_zarr, group=ice_type, consolidated=True).chunk(chunks)
                DS_list.append(ds)
            except (OSError, KeyError) as e:
                self.logger.warning(f"Skipping {P_zarr} ({ice_type}): {e}")
        if not DS_list:
            self.logger.warning(f"No {ice_type} datasets found in any Zarr group")
            return None
        DS_FI = xr.concat(DS_list, dim="time", coords='minimal')
        self.logger.info(f"Loaded {ice_type}: {len(DS_FI.time)} time steps from {len(DS_list)} files")
        if zarr_CICE:
            self.logger.info(f"Load monthly iceh_*.zarr files between {dt0_str} and {dtN_str}")
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
            self.logger.info(f"Found {len(P_monthly_zarrs)} zarr files")
            self.logger.debug(f"{[p.name for p in P_monthly_zarrs]}")
            CICE = xr.open_mfdataset(P_monthly_zarrs,
                                    engine     = "zarr",
                                    concat_dim = "time",
                                    combine    = "nested",
                                    parallel   = True,
                                    chunks     = chunks)
            if slice_hem:
                CICE = CICE.isel(nj=self.hemisphere_dict['nj_slice'])
        else:
            CICE = None
        return DS_FI, CICE

    def coarsen_and_align_simulated_FI_to_observed_FI(self, sim_ds, obs_ds, doy_vals=None, method="mean"):
        """
        Coarsen daily sim data into windows defined by AF2020 observation periods.

        INPUTS:
            sim_ds : xr.Dataset; daily CICE model output with time dimension.
            obs_ds : xr.Dataset; oservational fast ice dataset with `t_FI_obs` coordinate.
            doy_vals : list of int; DOY start values for each observation period. Defaults to AF2020 standard (24 bins).
            method : str; Aggregation method ("mean" or "median").

        OUTPUTS:
            Dataset with same shape as obs_ds[t_FI_obs] and matched time resolution.
        """
        if doy_vals is None:
            doy_vals = self.AF_FI_dict["DOY_vals"]
        sim_time = sim_ds["time"].values
        obs_times = pd.to_datetime(obs_ds["t_FI_obs"].values)
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
                ds_agg = ds_agg.expand_dims({"t_FI_obs": [dt_start]}).persist()
                grouped.append(ds_agg)
        if not grouped:
            raise ValueError("❌ No observation periods matched simulation data")
        ds_aligned = xr.concat(grouped, dim="t_FI_obs")
        self.logger.info(f"✅ Aligned model output to {len(ds_aligned.t_FI_obs)} obs windows")
        return ds_aligned


