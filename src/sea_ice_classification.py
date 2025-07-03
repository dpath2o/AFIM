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

    def boolean_fast_ice(self, FI_mask, dim="time", window=7, min_count=6, min_periods=1):
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
        self.logger.info(f"Rolling boolean presence: window = {window}, min_count = {min_count}")
        FI_roll_mask = (FI_mask.astype(bool).astype(int)
                               .rolling({dim: window}, center=True, min_periods=1)
                               .construct(f"{dim}_window")
                               .sum(dim=f"{dim}_window"))
        FI_bool = (FI_roll_mask >= min_count)
        return FI_bool

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


