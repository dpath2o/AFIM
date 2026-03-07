from __future__ import annotations
import xarray as xr
import numpy  as np
import pandas as pd
from pathlib  import Path
import os
__all__ = ["SeaIceCICE"]
class SeaIceCICE:

    def __init__():
        return

    def get_month_range(self, year_month):
        from datetime import datetime, timedelta
        dt0 = datetime.strptime(year_month + "-01", "%Y-%m-%d")
        dtN = (dt0.replace(day=28) + timedelta(days=4)).replace(day=1) - timedelta(days=1)
        return [dt0 + timedelta(days=i) for i in range((dtN - dt0).days + 1)]

    def verify_month(self,args_tuple):
        zarr_path, nc_dir, done_marker, dry_run = args_tuple
        year_month = zarr_path.stem.split("_")[1]
        if done_marker.exists():
            return f"[SKIP] {year_month}: already verified (.done exists)"
        dt_list = self.get_month_range(year_month)
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
        Verify monthly Zarr groups against expected daily NetCDFs and optionally delete NetCDFs.

        The procedure scans ``{self.D_zarr}/iceh_YYYY-MM.zarr`` directories, and for
        each month it checks:
        1) that the Zarr `time` coordinate covers all expected daily dates;
        2) that Zarr contains the union of variables present in the month’s
            remaining NetCDFs.

        A log of results is appended to ``{self.D_sim}/cleanup.log``. When
        `delete=True`, an interactive confirmation is prompted and, if accepted,
        the verified month’s NetCDFs are removed.

        Parameters
        ----------
        dry_run : bool, default True
            If False, write a ``.done_YYYY-MM`` marker in the Zarr directory after
            successful verification. When True, do not create markers.
        delete : bool, default False
            If True, prompt to delete verified NetCDFs after verification.
        max_workers : int, default 4
            Number of processes used for parallel month verification.

        Returns
        -------
        None
            Writes to log and may delete files on user confirmation.

        Notes
        -----
        - Verification is parallelised with ``ProcessPoolExecutor``.
        - Month logic expects Zarr directories named ``iceh_YYYY-MM.zarr`` and daily
        files at ``{self.D_sim}/history/daily/iceh.YYYY-MM-DD.nc``.
        - The verification routines are implemented as static helpers (`get_month_range`,
        `verify_month`).
        """
        from datetime import datetime
        self.D_iceh_nc = Path(self.D_sim,"history","daily")
        P_clean_log    = Path(self.D_sim,"cleanup.log")
        zarr_months    = sorted([p for p in self.D_zarr.glob("iceh_????-??.zarr") if p.is_dir()])
        tasks = []
        for zarr_path in zarr_months:
            ym          = zarr_path.stem.split("_")[1]
            done_marker = self.D_zarr / f".done_{ym}"
            tasks.append((zarr_path, self.D_iceh_nc, done_marker, dry_run))
        with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
            results = list(executor.map(self.verify_month, tasks))
        for res in results:
            self.logger.info(res)
        with open(P_clean_log, "a") as logf:
            logf.write("\n# Last updated: " + datetime.now().strftime("%Y-%m-%d %H:%M:%S") + "\n")
            for res in results:
                logf.write(res + "\n")
        if delete:
            self.logger.info("\nDeletion mode active.")
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
                self.logger.info("No deletable NetCDF files found.")
            else:
                self.logger.info(f"\n🔍 {len(total_files)} NetCDF files across {len(verified_months)} verified months are eligible for deletion.")
                confirm = input("Confirm delete all these files? [y/N] ").strip().lower()
                if confirm == "y":
                    for ym, files in verified_months:
                        for f in files:
                            try:
                                f.unlink()
                                self.logger.info(f"[DELETED] {f.name}")
                                log_entries.append(f"[DELETED] {f}")
                            except Exception as e:
                                self.logger.info(f"[ERROR] Could not delete {f.name}: {e}")
                                log_entries.append(f"[ERROR] Failed to delete {f}: {e}")
                    log_entries.append(f"# Deletion complete: {len(total_files)} files removed")
                else:
                    self.logger.info("Deletion cancelled.")
                    log_entries.append("# Deletion prompt declined — no files deleted")
            with open(P_clean_log, "a") as logf:
                for entry in log_entries:
                    logf.write(entry + "\n")

    def get_cice_files_between_dates(self, D_iceh, dt0_str, dtN_str):
        """
        List existing daily CICE iceh NetCDF files for a date range.

        Parameters
        ----------
        D_iceh : str or pathlib.Path
            Directory containing daily files named ``iceh.YYYY-MM-DD.nc``.
        dt0_str, dtN_str : str
            Inclusive date range in ``YYYY-MM-DD``.

        Returns
        -------
        list[pathlib.Path] or None
            Sorted list of existing files in the range. Returns ``None`` and logs
            an info message if none are found.
        """
        import re
        date_re = re.compile(r"(\d{4}-\d{2}-\d{2})\.nc$")
        D_iceh  = Path(D_iceh)
        dt0     = pd.to_datetime(dt0_str).date()
        dtN     = pd.to_datetime(dtN_str).date()
        files   = []
        for fpath in D_iceh.glob("*.nc"):
            m = date_re.search(fpath.name)
            if not m:
                continue
            fdate = pd.to_datetime(m.group(1)).date()
            if dt0 <= fdate <= dtN:
                files.append(fpath)
        if not files:
            self.logger.info(f"*NO* CICE daily NetCDF files ending in .YYYY-MM-DD.nc between {dt0_str} and {dtN_str}")
            return None
        return sorted(files)

    def delete_original_cice(self, P_orgs, P_iceh_zarr, m_str):
        """
        Delete the original daily NetCDF files for a given month if the monthly Zarr exists.

        Parameters
        ----------
        P_orgs : list[pathlib.Path]
            Paths to the month’s daily ``iceh.YYYY-MM-DD.nc`` files.
        P_iceh_zarr : str or pathlib.Path
            Path to the Zarr root (e.g., ``.../iceh_daily.zarr``).
        m_str : str
            Month string in ``YYYY-MM`` used as the Zarr group name.

        Behavior
        --------
        If ``{P_iceh_zarr}/.zgroup`` exists, the method attempts to delete each
        file in `P_orgs`. Failures are logged with a warning. If the Zarr store is
        incomplete (no ``.zgroup``), deletion is skipped with a warning.

        Returns
        -------
        None
        """
        grp = Path(P_iceh_zarr, m_str)
        # 1) require the month group to exist
        if not (grp / ".zgroup").exists():
            self.logger.warning(f"Zarr group {grp} missing/incomplete — skipping deletion")
            return
        # 2) verify open/read works (basic sanity)
        try:
            ds = xr.open_zarr(P_iceh_zarr, group=m_str, consolidated=True)
            ntime = ds.sizes.get("time", 0)
            ds.close()
        except Exception as e:
            self.logger.warning(f"Could not open Zarr group {m_str}; skipping deletion: {e}")
            return
        # 3) only then delete
        self.logger.info(f"Deleting original NetCDF files for {m_str} (ntime={ntime})")
        for f in P_orgs:
            try:
                Path(f).unlink()
            except Exception as e:
                self.logger.warning(f" Could not delete {f}: {e}")

    def daily_iceh_to_monthly_zarr(self,
                                   sim_name        = None,
                                   dt0_str         = None,
                                   dtN_str         = None,
                                   D_iceh          = None,
                                   netcdf_engine   = "scipy",
                                   overwrite       = None,
                                   delete_original = None):
        """
        Convert daily CICE ``iceh.YYYY-MM-DD.nc`` files into a monthly Zarr store.

        For each month intersecting the requested date window, this method:
        1) groups existing daily NetCDFs by month,
        2) opens them with ``xarray.open_mfdataset(..., engine='scipy')``,
        3) shifts the time coordinate **back by 1 day** (CICE writes 00:00:00 of
        the following day for daily averages),
        4) writes to a grouped Zarr store at ``{self.D_zarr}/iceh_daily.zarr``
        under group ``YYYY-MM`` (consolidated metadata),
        5) optionally deletes the original NetCDF files if the Zarr group exists.

        Parameters
        ----------
        sim_name : str, optional
            Simulation name. Defaults to ``self.sim_name``.
        dt0_str, dtN_str : str, optional
            Date window (inclusive) in ``YYYY-MM-DD``. Defaults to ``self.dt0_str``
            and ``self.dtN_str`` when omitted.
        D_iceh : str or pathlib.Path, optional
            Directory containing the per-day ``iceh.YYYY-MM-DD.nc`` files.
            Defaults to ``self.D_iceh``.
        overwrite : bool, optional
            If False and a monthly group already exists, the group is skipped (but
            originals may still be deleted if `delete_original` is True).
            Defaults to ``self.overwrite_zarr_group``.
        delete_original : bool, optional
            If True, remove the month’s original NetCDF files **after** successful
            Zarr write. Defaults to ``self.del_org_cice_iceh_nc``.

        Returns
        -------
        None
            Writes Zarr groups and logs progress; may delete input NetCDFs.

        Notes
        -----
        - The time shift of ``-1 day`` aligns Zarr entries with the *date of the
        average*, not the CICE timestamp.
        - Output chunks are set via ``self.CICE_dict['FI_chunks']`` after open.
        - The Zarr root is ``{self.D_zarr}/iceh_daily.zarr`` with one subgroup per month.
        """
        import re
        date_re = re.compile(r"(\d{4}-\d{2}-\d{2})\.nc$")
        from datetime    import datetime
        from collections import defaultdict
        sim_name = sim_name if sim_name is not None else self.sim_name
        dt0_str  = dt0_str  if dt0_str  is not None else self.dt0_str
        dtN_str  = dtN_str  or self.dtN_str
        D_iceh   = D_iceh   or self.D_iceh
        overwrite = overwrite if overwrite is not None else self.overwrite_zarr_group
        delete_nc = delete_original if delete_original is not None else self.del_org_cice_iceh_nc
        m_grps   = defaultdict(list)
        P_orgs   = self.get_cice_files_between_dates(D_iceh, dt0_str, dtN_str)
        if not P_orgs:
            self.logger.info("No CICE files found. Noting further to do here.")
            return
        for f in P_orgs:
           m = date_re.search(f.name)
           if not m:
               self.logger.warning(f"Skipping unrecognized filename: {f.name}")
               continue
           dt = datetime.strptime(m.group(1), "%Y-%m-%d")
           m_str = dt.strftime("%Y-%m")
           m_grps[m_str].append(f)
        P_iceh_zarr = Path(self.D_zarr, "iceh_daily.zarr")
        for m_str, P_ in m_grps.items():
            P_iceh_zarr_group = Path(P_iceh_zarr, m_str)
            if P_iceh_zarr_group.exists() and not overwrite:
                self.logger.info(f"Skipping existing {P_iceh_zarr_group}")
                if delete_nc:
                    self.delete_original_cice(P_, P_iceh_zarr, m_str)
                    continue
            else:
                self.logger.info(f"Loading NetCDF files for {m_str} via xarray mfdataset ...")
                CICE_all = xr.open_mfdataset(P_,
                                             engine     = netcdf_engine,
                                             parallel   = True,
                                             combine    = "nested",
                                             concat_dim = "time",
                                             coords     = "minimal",
                                             data_vars  = "minimal",
                                             compat     = "override",
                                             join       = "override",
                                             cache      = False,
                                             chunks     = self.CICE_dict["FI_chunks"])
                CICE_all = CICE_all.chunk(self.CICE_dict['FI_chunks'])
                self.logger.info(f"Subtracting one day from original dataset as CICE reports one day ahead for daily-averages")
                CICE_all["time"] = CICE_all["time"] - np.timedelta64(1, "D")
                self.logger.info(f"Writing {P_iceh_zarr} and group ('YYYY-MM'): {m_str}")
                CICE_all.to_zarr(P_iceh_zarr, group=m_str, mode="w", consolidated=True, zarr_format=2)
                self.get_dir_size(P_iceh_zarr_group)
                self.count_zarr_files(P_iceh_zarr_group)
                if delete_nc:
                    self.delete_original_cice(P_, P_iceh_zarr, m_str)

    def monthly_zarr_iceh_time_correction(self, P_mnthly_zarr, dry_run=True):
        """
        Correct the time coordinate in a monthly iceh Zarr store by shifting timestamps back by one day.

        CICE ice history files may encode timestamps at 00:00 that represent averages
        for the *previous* day. This utility rewrites a monthly Zarr store such that all
        time entries are shifted by -1 day.

        Parameters
        ----------
        P_mnthly_zarr : str or pathlib.Path
            Path to the monthly Zarr directory to correct (e.g., "iceh_2011-05.zarr").
        dry_run : bool, default True
            If True, do not write changes; only log what would be done.

        Notes
        -----
        - This method rewrites the store by writing a temporary Zarr and then replacing
        the original directory. Ensure you have permissions and sufficient disk space.
        - If timestamps already lie within the expected month bounds inferred from the
        filename, the store is skipped.

        Warnings
        --------
        This method may become obsolete if upstream conversion routines already perform
        the correct timestamp handling.
        """
        from calendar import monthrange
        import shutil
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
        ds.to_zarr(tmp_path, mode="w", consolidated=True, zarr_format=2)
        self.logger.info(f"Replacing original {P_mnthly_zarr.name}")
        shutil.rmtree(P_mnthly_zarr)
        tmp_path.rename(P_mnthly_zarr)
        self.logger.info(f"Time fix applied to {P_mnthly_zarr.name}")

    def correct_timestamp_for_all_monthly_zarr_iceh(self, sim_names=None, dry_run=True):
        """
        Apply `monthly_zarr_iceh_time_correction` across all monthly iceh Zarr stores.

        Parameters
        ----------
        sim_names : list[str], optional
            Simulation names to scan under the AFIM output root. If None, defaults to
            `[self.sim_name]`.
        dry_run : bool, default True
            If True, do not write changes; only log what would be done.

        Notes
        -----
        Monthly stores are identified using the glob pattern "iceh_????-??.zarr".
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

    def _get_zarr_group_index(self, zarr_root):
        """
        Index available Zarr groups under a root store and cache the results.

        This helper scans for groups of the form "YYYY-MM", determines overall time
        bounds by opening the first and last group, and detects whether consolidated
        metadata (`.zmetadata`) exists.

        Parameters
        ----------
        zarr_root : str or pathlib.Path
            Root path to a grouped Zarr store (expects "YYYY-MM" groups underneath).

        Returns
        -------
        tuple
            (groups, available_dt0, available_dtN, consolidated) where:
            - groups : list[str]
                Sorted list of available group names.
            - available_dt0 : pandas.Timestamp
                Earliest available timestamp in the first group.
            - available_dtN : pandas.Timestamp
                Latest available timestamp in the last group.
            - consolidated : bool
                True if `.zmetadata` exists at the root.

        Raises
        ------
        FileNotFoundError
            If no valid Zarr groups are found.

        Notes
        -----
        Results are cached per `zarr_root` string key in `self._zarr_group_index` to
        avoid repeated filesystem scans and metadata opens.
        """
        key = str(zarr_root)
        if not hasattr(self, "_zarr_group_index"):
            self._zarr_group_index = {}
        if key in self._zarr_group_index:
            return self._zarr_group_index[key]
        groups = sorted([p.name for p in zarr_root.glob("????-??")
                         if (zarr_root / p.name / ".zgroup").exists()])
        if not groups:
            raise FileNotFoundError(f"No Zarr groups found in {zarr_root}.")
        # Use consolidated metadata if present; else fall back
        consolidated = (zarr_root / ".zmetadata").exists()
        ds0 = xr.open_zarr(zarr_root, group=groups[0], consolidated=consolidated)
        dsN = xr.open_zarr(zarr_root, group=groups[-1], consolidated=consolidated)
        available_dt0 = pd.to_datetime(ds0.time.values[0])
        available_dtN = pd.to_datetime(dsN.time.values[-1])
        self._zarr_group_index[key] = (groups, available_dt0, available_dtN, consolidated)
        self.logger.info(f"Indexed Zarr root {zarr_root} (groups={len(groups)}, "
                         f"bounds={available_dt0.date()}..{available_dtN.date()}, consolidated={consolidated})")
        return self._zarr_group_index[key]

    def load_cice_zarr(self,
                        sim_name  = None,
                        dt0_str   = None,
                        dtN_str   = None,
                        D_alt     = None,
                        variables = None,
                        slice_hem = False):
        """
        Open and concatenate monthly-grouped CICE ice history Zarr data over a date window.

        This method loads only those monthly Zarr groups that intersect the requested
        date window, optionally subselects variables, crops each group to the requested
        time range, and concatenates along time.

        Parameters
        ----------
        sim_name : str, optional
            Simulation name. Defaults to `self.sim_name`.
        dt0_str, dtN_str : str, optional
            Requested analysis window in ``YYYY-MM-DD``. Defaults to `self.dt0_str` and
            `self.dtN_str`. Requests are clamped to available data.
        D_alt : str or pathlib.Path, optional
            Alternative simulation directory to use instead of `self.D_sim`.
        variables : list[str], optional
            If provided, select only these variables from each monthly group. Groups
            missing any requested variable are skipped with a warning.
        slice_hem : bool, default False
            If True, apply `self.slice_hemisphere(...)` after concatenation.

        Returns
        -------
        xarray.Dataset
            Concatenated dataset cropped to the intersection of requested and available
            times.

        Raises
        ------
        ValueError
            If no datasets remain after filtering/cropping (e.g., all groups skipped).

        Notes
        -----
        - Expects a grouped Zarr store at `self.D_iceh_zarr` with groups named "YYYY-MM".
        - Uses dask configuration overrides to control slicing and fusion behaviour.
        """
        import dask
        sim_name = sim_name or self.sim_name
        dt0_str  = dt0_str  or self.dt0_str
        dtN_str  = dtN_str  or self.dtN_str
        D_alt    = D_alt    or self.D_sim
        self.define_iceh_dirs(D_alt)
        zarr_root = self.D_iceh_zarr
        available_groups, available_dt0, available_dtN, consolidated = self._get_zarr_group_index(zarr_root)
        # Clamp user request to data availability
        user_dt0 = max(pd.to_datetime(dt0_str), available_dt0)
        user_dtN = min(pd.to_datetime(dtN_str), available_dtN)
        required_groups = [g for g in available_groups
                           if (pd.to_datetime(f"{g}-01") <= user_dtN) and 
                           (pd.to_datetime(f"{g}-01") + pd.offsets.MonthEnd(1) >= user_dt0)]
        self.logger.info(f"Loading Zarr groups between {user_dt0.date()} and {user_dtN.date()} ({len(required_groups)} groups)")
        with dask.config.set({
            "array.slicing.split_large_chunks": True,
            "array.chunk-size": "256MiB",
            # Consider turning fuse back on unless you have a specific reason to disable it:
            # "optimization.fuse.active": True,
            "optimization.fuse.active": False}):
            ds_list = []
            missing_anywhere = set()
            selected_anywhere = set()
            for g in required_groups:
                self.logger.debug(f"  - opening group {g}")
                ds = xr.open_zarr(zarr_root, group=g, consolidated=consolidated)
                if variables:
                    present = [v for v in variables if v in ds.data_vars]  # use data_vars, not coords
                    missing = [v for v in variables if v not in ds.data_vars]
                    if missing:
                        missing_anywhere.update(missing)
                        self.logger.warning(f"  [{g}] missing requested vars: {missing}")
                    if not present:
                        self.logger.warning(f"  > Skipping {g}: none of the requested variables are present.")
                        continue
                    selected_anywhere.update(present)
                    self.logger.info(f"  [{g}] selecting variables present: {present}")
                    ds = ds[present]
                ds = ds.sel(time=slice(user_dt0, user_dtN))
                ds_list.append(ds)
            if not ds_list:
                raise ValueError("No datasets to concatenate after filtering/cropping.")
            ds_all = xr.concat(ds_list, dim="time", coords="minimal", compat="override")
            # Warn once about variables that never appeared in any group
            if variables:
                never_found = [v for v in variables if v not in selected_anywhere]
                if never_found:
                    self.logger.warning(f"Requested variables not found in ANY group: {never_found}")
        if slice_hem:
            self.logger.info("  slicing hemisphere")
            ds_all = self.slice_hemisphere(ds_all)
        return ds_all
