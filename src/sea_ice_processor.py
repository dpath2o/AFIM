import json, os, shutil, sys, time, logging
import xarray as xr
import pandas as pd
import numpy  as np
import xesmf  as xe
from pathlib      import Path
from datetime     import datetime, timedelta
from collections  import defaultdict
sys.path.insert(0, '/home/581/da1339/AFIM/src/AFIM/src')
from grounded_iceberg_processor import GroundedIcebergProcessor
from dask.distributed import Client, LocalCluster
_dask_client = None  # <-- add this at the module level (outside the class)

class SeaIceProcessor:
    """
    Class to compute sea ice metrics from sea ice model output.

    Default mode of usage is for processing landfast sea ice from CICE model output using a rolling
    window climatology aligned with Fraser et al. (2020).

    Optionally can compute 'pack ice' that is either inclusive/exclusive of fast ice.

    This class is tightly coupled with:
    - A JSON configuration file providing paths, thresholds, and parameters:
      https://github.com/dpath2o/AFIM/blob/main/src/AFIM/src/JSONs/afim_cice_analysis.json
    - GroundedIcebergProcessor for land and iceberg masking
    - Observational climatology from Fraser et al. 2020 (DOI: 10.5194/essd-12-2987-2020):
      https://data.aad.gov.au/metadata/AAS_4116_Fraser_fastice_circumantarctic
    - Optional NSIDC pack ice data: https://nsidc.org/data/g02202/versions/4

    Supports both interactive in-memory workflows and Zarr output for batch processing.

    For more, see the full project repository:
    🔗 https://github.com/dpath2o/AFIM
    """
    def __init__(self, sim_name,
                 dt0_str                     = None,
                 dtN_str                     = None,
	             ice_concentration_threshold = None,
	             ice_speed_threshold         = None,
	             mean_period                 = None,
                 boolean_window              = None,
                 boolean_min_days            = None,
	             extra_cice_vars             = None,
	             hemisphere                  = None,
	             P_log                       = None,
	             P_json                      = None,
	             overwrite_AF2020_zarr       = False):
        """
		"""
        self.sim_name = sim_name
        if P_json is None:
            P_json = '/home/581/da1339/AFIM/src/AFIM/src/JSONs/afim_cice_analysis.json'
        with open(P_json, 'r') as f:
            self.config = json.load(f)
        self.sim_config     = self.config['sim_dict'][sim_name]
        self.CICE_dict      = self.config.get("CICE_dict", {})
        D_log               = self.config['D_dict']['logs']
        P_log               = P_log                       if P_log                       is not None else Path(D_log, f'SeaIceProcessor_FI_{self.sim_name}.log')
        self.dt0_str        = dt0_str                     if dt0_str                     is not None else self.config.get('dt0_str', '1993-01-01')
        self.dtN_str        = dtN_str                     if dtN_str                     is not None else self.config.get('dtN_str', '1999-12-31')
        self.mean_period    = mean_period                 if mean_period                 is not None else self.config.get('mean_period', 15)
        self.bool_window    = boolean_window              if boolean_window              is not None else self.config.get('bool_window', 7)
        self.bool_min_days  = boolean_min_days            if boolean_min_days            is not None else self.config.get('bool_min_days', 6)
        self.ispd_thresh    = ice_speed_threshold         if ice_speed_threshold         is not None else self.config.get('ice_speed_thresh_hi', 1e-3)
        self.icon_thresh    = ice_concentration_threshold if ice_concentration_threshold is not None else self.config.get('ice_conc_thresh', 0.15)
        self.cice_vars_ext  = extra_cice_vars             if extra_cice_vars             is not None else self.CICE_dict["FI_cice_vars_ext"]
        self.cice_vars_reqd = self.CICE_dict["FI_cice_vars_reqd"]
        self.cice_var_list  = self.cice_vars_reqd + self.cice_vars_ext
        self.D_sim          = Path(self.config['D_dict']['AFIM_out'], sim_name)
        self.D_iceh         = Path(self.D_sim , 'history', 'daily')
        self.D_zarr         = Path(self.D_sim , 'zarr')
        self.doy_vals       = self.config.get("DOY_vals", list(range(1, 366, 15)))
        self.sea_ice_dict   = self.config.get("sea_ice_dict", {})
        self.FIC_scale      = self.config.get('FIC_scale', 1e9)
        self.SIC_scale      = self.config.get('SIC_scale', 1e12)
        self.cm2m_fact      = self.config.get('cm2m_fact', 0.01)
        self.setup_logging(logfile=P_log)
        self.GI_proc= GroundedIcebergProcessor(P_json=P_json, sim_name=sim_name)
        self.GI_proc.load_bgrid()
        self.use_gi = self.GI_proc.use_gi
        if self.use_gi:
            self.GI_proc.compute_grounded_iceberg_area()
        self.GI_total_area = self.GI_proc.G_t['GI_total_area']  if self.use_gi else 0
        #self.GI_P_counts   = self.GI_proc.GI_P_counts if self.use_gi else 'GI not used in this simulation'
        self.P_KMT         = self.GI_proc.P_KMT_mod   if self.use_gi else self.CICE_dict['P_KMT']
        hemisphere         = hemisphere if hemisphere is not None else self.config.get('hemisphere', 'south')
        self.define_hemisphere(hemisphere)
        self.reG_weights_defined        = False
        self.reG_AF2020_weights_defined = False
        self.valid_ispd_types    = ["ispd_B", "ispd_Ta", "ispd_Tx", "ispd_BT"]
        self.valid_zarr_datasets = ['FI_B', 'FI_Ta', 'FI_Tx', 'FI_BT', 'PI_B', 'PI_Ta', 'PI_Tx', 'PI_BT', 'SO']

    def setup_logging(self, logfile=None):
        self.logger = logging.getLogger(self.sim_name)
        self.logger.setLevel(logging.INFO)
        if not self.logger.handlers:
            formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
            ch = logging.StreamHandler()
            ch.setFormatter(formatter)
            self.logger.addHandler(ch)
            if logfile:
                if os.path.exists(logfile):
                    os.remove(logfile)
                fh = logging.FileHandler(logfile)
                fh.setFormatter(formatter)
                self.logger.addHandler(fh)
                self.logger.info(f"log file intialised: {logfile}")

    def define_hemisphere(self, hemisphere):
        if hemisphere.lower() in ['north', 'northern', 'nh', 'n', 'no']:
            self.hemisphere_geographic_extent = [0, 360, 0, 90]
            self.hemisphere_map_extent        = [-180,180,55,90]
            self.hemisphere_projection        = "S0.0/90.0/50/15C"
            self.hemisphere_map_text_location = [-120,56]
            self.hemisphere_abbreviation      = 'NH'
            self.hemisphere_nj_slice          = slice(540,1080)
            self.hemisphere                   = 'north'
        elif hemisphere.lower() in ['south', 'southern', 'sh', 's', 'so']:
            self.hemisphere_geographic_extent = [0, 360, -90, 0]
            self.hemisphere_map_extent        = [-180,180,-90,-55]
            self.hemisphere_projection        = "S0.0/-90.0/50/15C"
            self.hemisphere_map_text_location = [0,-90]
            self.hemisphere_abbreviation      = 'SH'
            self.hemisphere_nj_slice          = slice(0,540)
            self.hemisphere                   = 'south'
        else:
            raise ValueError(f"Invalid hemisphere '{hemisphere}'. Valid options are: "
                             "['north', 'south', 'northern', 'southern', 'sh', 'nh', 'SH', 'NH']")
        self.logger.info(f"hemisphere initialised: {self.hemisphere}")

    def slice_hemisphere(self, var_dict):
        sliced = {k: v.isel(nj=self.hemisphere_nj_slice) if k not in {'FI_OBS_CLI', 'preserved_vars'} else v for k, v in var_dict.items()}
        self.logger.info("hemisphere sliced on 'main' dataset")
        if "preserved_vars" in var_dict and isinstance(var_dict["preserved_vars"], dict):
            preserved = var_dict["preserved_vars"]
            sliced_preserved = { k: v.isel(nj=self.hemisphere_nj_slice) for k, v in preserved.items() }
            sliced["preserved_vars"] = sliced_preserved
            self.logger.info("hemisphere sliced on 'preserved' dataset")
        return sliced

    def normalise_longitudes(self,lon):
        lon_norm = (lon+360)%360
        return lon_norm

    def build_grid_dict(self, lat, lon):
        ny, nx               = lat.shape
        lat_b                = np.zeros((ny + 1, nx + 1))
        lon_b                = np.zeros((ny + 1, nx + 1))
        lat_b[1:-1, 1:-1]    = 0.25 * (lat[:-1, :-1] + lat[:-1, 1:] + lat[1:, :-1] + lat[1:, 1:])
        lon_b[1:-1, 1:-1]    = 0.25 * (lon[:-1, :-1] + lon[:-1, 1:] + lon[1:, :-1] + lon[1:, 1:])
        lat_b[[0, -1], 1:-1] = lat[[0, -1], :-1]
        lon_b[[0, -1], 1:-1] = lon[[0, -1], :-1]
        lat_b[1:-1, [0, -1]] = lat[:-1, [0, -1]]
        lon_b[1:-1, [0, -1]] = lon[:-1, [0, -1]]
        lat_b[0, 0]          = lat[0, 0]
        lat_b[0, -1]         = lat[0, -1]
        lat_b[-1, 0]         = lat[-1, 0]
        lat_b[-1, -1]        = lat[-1, -1]
        lon_b[0, 0]          = lon[0, 0]
        lon_b[0, -1]         = lon[0, -1]
        lon_b[-1, 0]         = lon[-1, 0]
        lon_b[-1, -1]        = lon[-1, -1]
        return {"lat"  : lat,
                "lon"  : self.normalise_longitudes(lon),
                "lat_b": lat_b,
                "lon_b": self.normalise_longitudes(lon_b)}

    def define_reG_weights(self):
        """
        Define xESMF regridding weights from B-grid (u) to T-grid.
        Generates and saves weights if not found on disk.
        """
        #G_u           = self.build_grid_dict(self.GI_proc.G_u['lat'], self.GI_proc.G_u['lon'])
        #G_t           = self.build_grid_dict(self.GI_proc.G_t['lat'], self.GI_proc.G_t['lon'])
        G_u           = xr.Dataset()
        G_u['lat']    = self.GI_proc.G_u['lat']
        G_u['lon']    = self.GI_proc.G_u['lon']
        G_u['lat_b']  = self.GI_proc.G_u['lat_b']
        G_u['lon_b']  = self.GI_proc.G_u['lon_b']
        G_t           = xr.Dataset()
        G_t['lat']    = self.GI_proc.G_t['lat']
        G_t['lon']    = self.GI_proc.G_t['lon']
        G_t['lat_b']  = self.GI_proc.G_t['lat_b']
        G_t['lon_b']  = self.GI_proc.G_t['lon_b']
        G_t["mask"]   = self.GI_proc.G_t['kmt_mod']
        F_weights     = self.CICE_dict["P_reG_u2t_weights"]
        weights_exist = os.path.exists(F_weights)
        self.logger.info(f"{'🔁 Reusing' if weights_exist else '⚙️ Creating'} regrid weights: {F_weights}")
        self.reG = xe.Regridder(G_u, G_t,
                                method            = "bilinear",
                                periodic          = True,
                                ignore_degenerate = True,
                                extrap_method     = "nearest_s2d",
                                reuse_weights     = weights_exist,
                                filename          = F_weights)
        self.reG_weights_defined = True

    def reG_bgrid_to_tgrid_xesmf(self, ds):
        """
        Regrid all B-grid variables (with 'uarea' as cell_measures) to T-grid using batch processing.
        """
        vars_to_reG = []
        for var in ds.data_vars:
            if var == 'tarea':
                continue
            if "cell_measures" in ds[var].attrs and "uarea" in ds[var].attrs["cell_measures"]:
                self.logger.info(f"regridding variable: {var}")
                vars_to_reG.append(var)
        reG_input  = [ds[v] for v in vars_to_reG]
        ds_to_reG  = ds[vars_to_reG]
        reG_out    = self.reG(ds_to_reG, skipna=True)
        for var in reG_out.data_vars:
            reG_out[var].attrs = ds[var].attrs.copy()
        untouched_vars      = {v: ds[v] for v in ds.data_vars if v not in vars_to_reG}
        ds_reG              = xr.Dataset({**reG_out.data_vars, **untouched_vars})
        ds_reG['lat']       = ds['ULAT']
        ds_reG['lon']       = ds['ULON']
        ds_reG['lat'].attrs = ds['ULAT'].attrs
        ds_reG['lon'].attrs = ds['ULON'].attrs
        ds_reG              = ds_reG.drop_vars(['tarea', 'TLAT', 'TLON', 'ULAT', 'ULON'])
        return ds_reG

    def compute_ice_speed_types(self, DS, ispd_type, temporally_average=False, mean_period=None):
        mean_period = mean_period if mean_period is not None else self.mean_period
        if "ispd_B" in ispd_type or "ispd_Ta" in ispd_type or "ispd_BT" in ispd_type:
            self.logger.info("⚡ Computing B-grid ice speed ('ispd_B') and internal ice stress from strintx and strinty ('ists_B')")
            ispd_B = xr.apply_ufunc(np.hypot, DS['uvel']   , DS['vvel']   , dask="allowed", output_dtypes=[DS['uvel'].dtype])
            ists_B = xr.apply_ufunc(np.hypot, DS['strintx'], DS['strinty'], dask="allowed", output_dtypes=[DS['strintx'].dtype])
            if temporally_average:
                self.logger.info("⚡ Temporally-Averaging ispd_B and ists_B")
                ispd_B = ispd_B.rolling(time=mean_period, center=True, min_periods=1).mean()
                ists_B = ists_B.rolling(time=mean_period, center=True, min_periods=1).mean()
            DS['ispd_B'] = ispd_B
            if "ispd_Ta" in ispd_type or "ispd_BT" in ispd_type:
                self.logger.info("🌀 Spatially-Averaging (re-griddding) ispd_B and ists_B to T-grid: CREATING ispd_Ta and ists_Ta")
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
                                                                "TLAT" : (("nj", "ni"), DS.TLAT.values),
                                                                "TLON" : (("nj", "ni"), DS.TLON.values)},
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
                                                                "TLAT" : (("nj", "ni"), DS.TLAT.values),
                                                                "TLON" : (("nj", "ni"), DS.TLON.values)},
                                                    attrs  = {"long_name" : "T-grid interpolated B-grid internal ice stress",
                                                                "units"     : ists_B.attrs.get("units", "N/m^2")})
        if "ispd_Tx" in ispd_type or "ispd_BT" in ispd_type:
            self.logger.info("🌀 xESMF regrid to T-grid")
            if not self.reG_weights_defined:
                self.define_reG_weights()
            DS_reG = self.reG_bgrid_to_tgrid_xesmf(DS)
            self.logger.info("⚡ Computing xESMF-based ice speed ('ispd_Tx') and internal ice stress from strintx and strinty ('ists_Tx')")
            ispd_Tx = xr.apply_ufunc(np.hypot, DS_reG['uvel']   , DS_reG['vvel']   , dask="allowed", output_dtypes=[DS_reG['uvel'].dtype])
            ists_Tx = xr.apply_ufunc(np.hypot, DS_reG['strintx'], DS_reG['strinty'], dask="allowed", output_dtypes=[DS_reG['strintx'].dtype])
            if temporally_average:
                self.logger.info("⚡ Temporally-Averaging ispd_Tx and ists_Tx")
                ispd_Tx = ispd_Tx.rolling(time=mean_period, center=True, min_periods=1).mean()
                ists_Tx = ispd_Tx.rolling(time=mean_period, center=True, min_periods=1).mean()
            DS['ispd_Tx'] = ispd_Tx.copy()
            DS['ists_Tx'] = ists_Tx.copy()
        return DS

    def compute_composite_ice_speed(self, DS):
        """
        Compute element-wise average of all three regridded ice speed fields: ispd_B, ispd_Ta, ispd_Tx.
        Requires these variables to already exist in the dataset.
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

    def reapply_landmask(self, DS):
        kmt_mask = self.GI_proc.G_t['kmt_mod'] == 1  # True for ocean, False for land
        for var in DS.data_vars:
            da = DS[var]
            if {"nj", "ni"}.issubset(da.dims):  # Only apply to spatial fields
                self.logger.debug(f"🏝️ Masking land for variable: {var}")
                DS[var] = da.where(kmt_mask)
        self.logger.info("✅ Applied landmask to rolled dataset")
        return DS

    def create_fast_ice_mask(self, DS, ispd_type, ispd_thresh):
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
        DS_grouped = defaultdict(lambda: {k: [] for k in self.valid_zarr_datasets})        
        for group in masks:
            fi_mask          = masks[group]
            pi_mask          = ~fi_mask
            ds_fi            = DS.where(fi_mask)
            ds_pi            = DS.where(pi_mask)
            ds_fi['FI_mask'] = fi_mask
            ds_pi['PI_mask'] = pi_mask
            pi_group         = group.replace("FI_", "PI_")
            DS_grouped[m_str][group].append(ds_fi)
            DS_grouped[m_str][pi_group].append(ds_pi)
        DS_grouped[m_str]['SO'].append(DS)
        return DS_grouped

    def write_to_zarr(self, DS_grouped, P_zarr_root, ispd_thresh, ispd_type, m_str):
        mask_type_map = {
            "FI_B"  : "fast ice mask based on thresholding ispd_B (native U-grid)",
            "FI_Ta" : "fast ice mask based on thresholding ispd_Ta (spatially averaged ispd_B)",
            "FI_Tx" : "fast ice mask based on thresholding ispd_Tx (xESMF regridded uvel/vvel)",
            "FI_BT" : "fast ice mask based on thresholding ispd_BT (composite of B, Ta, Tx)",
            "PI_B"  : "pack ice mask as complement to FI_B",
            "PI_Ta" : "pack ice mask as complement to FI_Ta",
            "PI_Tx" : "pack ice mask as complement to FI_Tx",
            "PI_BT" : "pack ice mask as complement to FI_BT"
        }
        month_groups = defaultdict(lambda: {k: [] for k in self.valid_zarr_datasets})
        for group, datasets in DS_grouped[m_str].items():
            if not datasets:
                continue  # Skip empty groups
            if group != "SO":
                ispd_var = group.replace("FI_", "ispd_").replace("PI_", "ispd_")
                if ispd_var not in ispd_type:
                    continue
            # Skip if group already exists and overwriting is not allowed
            if not self.overwrite_zarr_group:
                try:
                    zarr_root = xr.open_zarr(P_zarr_root, consolidated=True)
                    if group in zarr_root:
                        self.logger.info(f"⏩ Skipping group {group} (already exists and overwrite disabled)")
                        continue
                except (FileNotFoundError, KeyError, ValueError):
                    pass
            self.logger.info(f"💾 Writing group '{group}' to Zarr store: {P_zarr_root}")
            try:
                ds_monthly = xr.concat(datasets, dim="time").sortby("time")
                if group != "SO":
                    ds_monthly.attrs["ispd_thresh"] = ispd_thresh
                    ds_monthly.attrs["mask_type"]   = mask_type_map.get(group, "unknown")
                ds_monthly = ds_monthly.chunk({"time": -1, "nj": 540, "ni": 1440})
                ds_monthly.to_zarr(P_zarr_root, group=group, mode="a", consolidated=True)
            except Exception as e:
                self.logger.error(f"❌ Failed to write group {group}: {e}")
        DS_grouped[m_str].clear()


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
            raise FileNotFoundError(f"No iceh files found between {dt0_str} and {dtN_str} in {D_iceh}")
        return sorted(files)

    def drop_unwanted_ispd_vars(self, ds, ispd_type_requested):
        keep_vars = set(ispd_type_requested) | set(ds.data_vars) - set(self.valid_ispd_types)
        drop_vars = [v for v in self.valid_ispd_types if v in ds and v not in keep_vars]
        if drop_vars:
            self.logger.info(f"🧹 Dropping unused ice speed variables: {drop_vars}")
            return ds.drop_vars(drop_vars)
        return ds


    def process_daily_cice(self,
                           sim_name             = None,
                           dt0_str              = None,
                           dtN_str              = None,
                           D_iceh               = None,
                           ispd_thresh          = None,
                           ispd_type            = None,
                           overwrite_zarr_group = False):
        sim_name      = sim_name    or self.sim_name
        dt0_str       = dt0_str     or self.dt0_str
        dtN_str       = dtN_str     or self.dtN_str
        D_iceh        = D_iceh      or self.D_iceh
        ispd_type_req = ispd_type   or self.valid_ispd_types
        ispd_thresh   = float(ispd_thresh) if ispd_thresh is not None else self.ispd_thresh
        ispd_thresh_str = f"{ispd_thresh:.1e}".replace("e-0", "e-")
        self.overwrite_zarr_group = overwrite_zarr_group
        if isinstance(ispd_type_req, str):
            ispd_type_req = [ispd_type_req]
        assert all(t in set(self.valid_ispd_types) for t in ispd_type_req), f"❌ Invalid requested sea ice speed 'type': {ispd_type_req}.\n\tMust be one or more of these valid types: {self.valid_ispd_types}"
        dts  = pd.date_range(dt0_str, dtN_str, freq="D")
        m_DS = defaultdict(lambda: {k: [] for k in self.valid_zarr_datasets})
        for dt in dts:
            date_str = dt.strftime("%Y-%m-%d")
            m_str    = dt.strftime("%Y-%m")
            P_iceh   = Path(D_iceh, f"iceh.{date_str}.nc")
            P_zarr   = Path(self.D_zarr, f"ispd_thresh_{ispd_thresh_str}", f"cice_daily_{m_str}.zarr")
            if not P_iceh.exists():
                self.logger.info(f"❌ Missing CICE file: {P_iceh}")
                continue
            self.logger.info(f"📂 Loading: {P_iceh}")
            CICE                = xr.open_dataset(P_iceh, engine="netcdf4")[list(self.cice_var_list)]
            CICE_ispd           = self.compute_ice_speed_types( CICE , ispd_type_req )
            CICE_reM            = self.reapply_landmask( CICE_ispd )
            CICE_reM['ispd_BT'] = self.compute_composite_ice_speed( CICE_reM )
            CICE_reM            = self.drop_unwanted_ispd_vars(CICE_reM, ispd_type_req)
            self.logger.info("🌊 Subsetting Ocean into either southern or northern hemisphere (default: southern)")
            CICE_SO = CICE_reM.isel(nj=self.hemisphere_nj_slice)
            self.logger.info("❄️ create fast ice masks")
            masks = self.create_fast_ice_mask( CICE_SO , ispd_type_req , ispd_thresh )
            self.logger.info("❄️ apply fast ice masks to dataset")
            CICE_grouped = self.groupby_fast_ice_masks( CICE_SO, masks, m_str)
            for key in CICE_grouped[m_str]:
                m_DS[m_str][key].extend(CICE_grouped[m_str][key])
            if dt == dts[-1] or dt.month != (dt + pd.Timedelta(days=1)).month:
                self.write_to_zarr(m_DS, P_zarr, ispd_thresh, ispd_type_req, m_str)
                m_DS[m_str].clear()

    def process_rolling_cice(self,
                             sim_name             = None,
                             dt0_str              = None,
                             dtN_str              = None,
                             D_iceh               = None,
                             mean_period          = None,
                             ispd_thresh          = None,
                             ispd_type            = None,
                             overwrite_zarr_group = False):
        sim_name     = sim_name    or self.sim_name
        dt0_str      = dt0_str     or self.dt0_str
        dtN_str      = dtN_str     or self.dtN_str
        D_iceh       = D_iceh      or self.D_iceh
        mean_period  = mean_period or self.mean_period
        ispd_type    = ispd_type   or self.valid_ispd_types
        ispd_thresh  = float(ispd_thresh) if ispd_thresh is not None else self.ispd_thresh
        ispd_thresh_str = f"{ispd_thresh:.1e}".replace("e-0", "e-")
        self.overwrite_zarr_group = overwrite_zarr_group
        if isinstance(ispd_type, str):
            ispd_type = [ispd_type]
        valid_types = set(self.valid_ispd_types)
        assert all(t in valid_types for t in ispd_type), f"❌ Invalid ispd_type: {ispd_type}. Must be one or more of {self.valid_ispd_types}"
        dts = pd.date_range(dt0_str, dtN_str, freq="D")
        self.logger.info(f"📆 Rolling mean first, then fast ice masking between {dt0_str} and {dtN_str}")
        self.logger.info("📚 Loading entire daily dataset with variables: %s", self.cice_var_list)
        try:
            files = self.get_cice_files_between_dates(D_iceh, dt0_str, dtN_str)
            self.logger.info(f"📦 Loading {len(files)} CICE history files from {files[0].name} to {files[-1].name}")
            def _preprocess_cice(ds):
                return ds[list(self.cice_var_list)]
            CICE_all = xr.open_mfdataset(files,
                                         engine="scipy",#"netcdf4",
                                         parallel=True, 
                                         preprocess=_preprocess_cice, 
                                         combine="by_coords", 
                                         cache=True, 
                                         chunks={})
            CICE_all = CICE_all.chunk({'time': mean_period * 2})
        except Exception as e:
            self.logger.error(f"Failed to load datasets: {e}")
            return
        CICE_ispd = self.compute_ice_speed_types( CICE_all , ispd_type , temporally_average=True, mean_period=mean_period)
        self.logger.info("🌀 Computing selective rolling means")
        CICE_roll_vars = {}
        for var in CICE_ispd.data_vars:
            da = CICE_ispd[var]
            # Skip ice speed fields
            if var.endswith(("_B", "_Ta", "_Tx")):
                self.logger.debug(f"⏭️ Skipping rolling mean for {var} (already derived)")
                CICE_roll_vars[var] = da
                continue
            # Skip vars with uarea-based cell measures
            cell_meas = da.attrs.get("cell_measures", "")
            if "area: uarea" in cell_meas:
                self.logger.info(f"⏭️ Skipping temporal mean for {var} due to 'cell_measures = {cell_meas}'")
                CICE_roll_vars[var] = da
                continue
            # Apply rolling mean
            self.logger.info(f"📉 Rolling mean on variable: {var}")
            CICE_roll_vars[var] = da.rolling(time=mean_period, center=True, min_periods=1).mean()
        CICE_roll = xr.Dataset(CICE_roll_vars, coords=CICE_ispd.coords)
        # Adjust for CICE time convention: each time represents the day *after* the average
        CICE_roll['time'] = CICE_roll['time'] - np.timedelta64(1, 'D')
        CICE_roll = CICE_roll.where(~np.isnan(CICE_roll['aice']), drop=False).compute()
        CICE_roll = CICE_roll.dropna(dim="time", how="all", subset=["aice"])  # drop only fully-NaN time steps
        CICE_reM  = self.reapply_landmask( CICE_roll )
        CICE_reM['ispd_BT'] = self.compute_composite_ice_speed( CICE_reM )
        self.logger.info(f"📊 Original time steps: {CICE_all.time.size}")
        self.logger.info(f"📊 Rolled time steps: {CICE_reM.time.size}")        
        time_vals = pd.to_datetime(CICE_reM.time.values)
        mo_str_da = xr.DataArray(time_vals.strftime("%Y-%m"), coords={"time": CICE_reM.time}, dims="time")
        mo_uni    = np.unique(mo_str_da.values)
        for m_str in mo_uni:
            P_zarr = Path(self.D_zarr, f"ispd_thresh_{ispd_thresh_str}", f"cice_rolling_{m_str}.zarr")
            CICE_month = CICE_reM.sel(time=mo_str_da == m_str)
            if CICE_month.time.size == 0:
                self.logger.warning(f"⚠️ No data for month: {m_str}")
                continue
            CICE_SO = CICE_month.isel(nj=self.hemisphere_nj_slice)
            self.logger.info(f"📆 Rolling monthly group: {m_str} with {CICE_SO.time.size} time steps")
            masks = self.create_fast_ice_mask(CICE_SO, ispd_type, ispd_thresh)
            CICE_grouped = self.groupby_fast_ice_masks(CICE_SO, masks, m_str)
            self.write_to_zarr(CICE_grouped, P_zarr, ispd_thresh, ispd_type, m_str)

    def load_processed_cice(self, 
                            sim_name    = None,
                            rolling     = False,
                            ispd_thresh = None,
                            ice_type    = "FI_B",
                            dt0_str     = None,
                            dtN_str     = None,
                            D_zarr      = None,
                            chunks      = None):
        sim_name        = sim_name    or self.sim_name
        ispd_thresh     = ispd_thresh or self.ispd_thresh
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
        if not P_zarrs:
            self.logger.warning(f"⚠️ No Zarr datasets found in {D_zarr}")
            return None
        datasets = []
        for P_zarr in P_zarrs:
            self.logger.debug(f"attempting to load: {P_zarr}")
            try:
                ds = xr.open_zarr(P_zarr, group=ice_type, consolidated=True)
                datasets.append(ds)
            except (OSError, KeyError) as e:
                self.logger.warning(f"⚠️ Skipping {P_zarr} ({ice_type}): {e}")
        if not datasets:
            self.logger.warning(f"⚠️ No {ice_type} datasets found in any Zarr group")
            return None
        ds_all = xr.concat(datasets, dim="time")
        ds_all = ds_all.chunk(chunks)
        self.logger.info(f"✅ Loaded {ice_type}: {len(ds_all.time)} time steps from {len(datasets)} files")
        return ds_all

    def compute_rolling_mean_on_dataset(self, ds, mean_period=None):
        mean_period = mean_period if mean_period is not None else self.mean_period
        return ds.rolling(time=mean_period, center=True, min_periods=1).mean()

    def boolean_fast_ice(self, FI_mask, dim="time", window=7, min_count=6):
        window    = window    if window    is not None else self.bool_window
        min_count = min_count if min_count is not None else self.bool_min_days
        self.logger.info(f"🔁 Rolling boolean presence: window = {window}, min_count = {min_count}")
        FI_roll_mask = FI_mask.rolling({dim: window}, center=True).construct(f"{dim}_window").sum(dim=f"{dim}_window")
        FI_bool      = (FI_roll_mask >= min_count).persist()
        return FI_bool

    def compute_ice_area(self, SIC, GC_area, ice_area_scale=None):
        ice_area_scale = ice_area_scale if ice_area_scale is not None else self.FIC_scale
        self.logger.info(f"🧮 Spatially-integrating the product of sea ice concentrations and grid cell areas")
        IA = ((SIC * GC_area).sum(dim=("nj", "ni"))).persist()
        IA = IA+self.GI_total_area
        IA = IA/ice_area_scale
        return IA

    def compute_variable_aggregate(self, da, time_coord_name='time'):
        return da.sum(dim=time_coord_name) / da[time_coord_name].sizes.get(time_coord_name, 1)

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
        sectors     = ['circumpolar', 'IOsector', 'WPOsector', 'RSsector', 'BASsector', 'WSsector']
        self.logger.info(f"loaded: {csv_path}")
        return row[sectors].reset_index(drop=True)

    def define_AF2020_reG_weights(self, FI_obs_native):
        from pyproj import CRS, Transformer
        self.logger.info("define model grid")
        G_t           = xr.Dataset()
        G_t['lat']    = self.GI_proc.G_t['lat']
        G_t['lon']    = self.GI_proc.G_t['lon']
        G_t['lat_b']  = self.GI_proc.G_t['lat_b']
        G_t['lon_b']  = self.GI_proc.G_t['lon_b']
        G_t["mask"]   = self.GI_proc.G_t['kmt_mod']
        self.logger.info("defining AF2020 regridder weights")
        F_weights     = self.sea_ice_dict["AF_reG_weights"]
        weights_exist = os.path.exists(F_weights)
        self.logger.info("convert 'AF_FI_OBS_2020db' Cartesian coordinates to spherical coordindates")
        crs_obs          = CRS.from_epsg(self.sea_ice_dict["projection_FI_obs"]) #unique to observations
        crs_spherical    = CRS.from_epsg(self.sea_ice_dict["projection_wgs84"])  #spherical coordinates
        transformer      = Transformer.from_crs(crs_obs, crs_spherical, always_xy=True)
        X, Y             = np.meshgrid(FI_obs_native['x'].isel(time=0).values, FI_obs_native['y'].isel(time=0).values)
        lon_obs, lat_obs = transformer.transform(X,Y)
        self.G_obs       = self.GI_proc.build_grid_dict(lat_obs, lon_obs)
        self.logger.info(f"{'🔁 Reusing' if weights_exist else '⚙️ Creating'} regrid weights: {F_weights}")
        self.reG_AF2020 = xe.Regridder(self.G_obs, G_t,
                                       method            = "bilinear",
                                       periodic          = True,
                                       ignore_degenerate = True,
                                       extrap_method     = "nearest_s2d",
                                       reuse_weights     = weights_exist,
                                       filename          = F_weights)

    def load_AF2020_FI_org_netcdf(self, P_orgs):
        FI_obs = xr.open_mfdataset(P_orgs, engine='netcdf4')
        self.define_AF2020_reG_weights(FI_obs)
        self.logger.info("*** Regridding 'AF_FI_OBS_2020db' to CICE T-grid *** ")
        self.FI_obs_reG_lon = self.reG_AF2020(self.G_obs["lon"]).compute()
        self.FI_obs_reG_lat = self.reG_AF2020(self.G_obs["lat"]).compute()
        self.FI_obs_reG_dat = self.reG_AF2020( FI_obs["Fast_Ice_Time_series"] ).compute()
        mask     = xr.where(self.FI_obs_reG_dat >= 4, 1, np.nan)
        FI       = (('t_FI_obs', 'nj', 'ni'),
                    self.FI_obs_reG_dat.where(mask).values,
                    {'long_name': FI_obs['Fast_Ice_Time_series'].attrs['long_name']})
        t_alt    = (('t_FI_obs'),
                    FI_obs.date_alt.values,
                    {'long_name'   : FI_obs.date_alt.attrs['long_name'],
                     'description' : FI_obs.date_alt.attrs['description']})
        t_coords = (('t_FI_obs'),
                    FI_obs.time.values,
                    {'description' : "Start date of 15- or 20-day image mosaic window.",
                     'units'       : "days since 2000-1-1 0:0:0"})
        x_coords = (('nj','ni'),
                    self.FI_obs_reG_lon.values,
                    {'long_name': 'longitude',
                     'units'    : 'degrees_north'})
        y_coords = (('nj','ni'),
                    self.FI_obs_reG_lat.values,
                    {'long_name': 'latitude',
                     'units'    : 'degrees_east'})
        self.logger.info("converted AF2020 database for use with SeaIceProcessor")
        return xr.Dataset({'FI'       : FI,
                           'FI_t_alt' : t_alt },
                          coords=dict(t_FI_obs=t_coords, obs_lon=x_coords, obs_lat=y_coords)).compute()

    def filter_AF2020_FI_by_date(self, start_date, end_date):
        D_obs     = Path(self.sea_ice_dict['D_AF2020_db_org'])
        yrs_reqd  = set([start_date.year, end_date.year])
        P_orgs    = [D_obs / f"FastIce_70_{yr}.nc" for yr in yrs_reqd]
        ds        = self.load_AF2020_FI_org_netcdf(P_orgs)
        alt_dates = pd.to_datetime(ds['FI_t_alt'].values.astype(str), format='%Y%m%d')
        valid_idx = (alt_dates >= pd.to_datetime(start_date)) & (alt_dates <= pd.to_datetime(end_date))
        matched   = ds.sel(t_FI_obs=valid_idx).persist()
        if matched.dims['t_FI_obs'] == 0:
            self.logger.warning(f"No matching observational dates found between {start_date} and {end_date}")
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
        alt_dates = pd.to_datetime(ds_all['FI_t_alt'].values.astype(str), format='%Y%m%d')
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
