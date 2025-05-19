import json, os, shutil, sys, time, logging
import xarray as xr
import pandas as pd
import numpy  as np
#import xesmf  as xe
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
        self.gi_processor = GroundedIcebergProcessor(P_json=P_json, sim_name=sim_name)
        self.gi_processor.load_grid_and_landmask()
        self.use_gi = self.gi_processor.use_gi
        if self.use_gi:
            self.gi_processor.load_AFIM_GI()
        self.GI_total_area = self.gi_processor.total_area  if self.use_gi else 0
        self.GI_P_counts   = self.gi_processor.GI_P_counts if self.use_gi else 'GI not used in this simulation'
        self.P_KMT         = self.gi_processor.P_KMT_mod   if self.use_gi else self.CICE_dict['P_KMT']
        hemisphere         = hemisphere if hemisphere is not None else self.config.get('hemisphere', 'south')
        self.define_hemisphere(hemisphere)
        self.reG_weights_defined        = False
        self.reG_AF2020_weights_defined = False
        self.datasets_by_group = {'FI_B'  : defaultdict(list),
                                  'FI_Ta' : defaultdict(list),
                                  'FI_Tx' : defaultdict(list),
                                  'PI_B'  : defaultdict(list),
                                  'PI_Ta' : defaultdict(list),
                                  'PI_Tx' : defaultdict(list),
                                  'SO'    : defaultdict(list)}
        self.mask_type_map = {"FI_B"  : "fast ice mask based on thresholding sea ice speed *without* regridding (i.e. on the native B-grid)",
                              "FI_Ta" : "fast ice mask based on thresholding sea ice speed that has been first spatial-averaged to the T-grid",
                              "FI_Tx" : "fast ice mask based on thresholding sea ice speed that has been first re-gridded to the T-grid based on xESMF re-gridding weights",
                              "PI_B"  : "pack ice mask as complement to FI_B",
                              "PI_Ta" : "pack ice mask as complement to FI_Ta",
                              "PI_Tx" : "pack ice mask as complement to FI_Tx"}        

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
        G_u           = self.build_grid_dict(self.gi_processor.G_u['lat'], self.gi_processor.G_u['lon'])
        G_t           = self.build_grid_dict(self.gi_processor.G_t['lat'], self.gi_processor.G_t['lon'])
        G_u["mask"]   = self.gi_processor.KMT_mod
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

    def compute_ice_speed_types(self, DS, ispd_type):
        if "ispd_B" in ispd_type:
            self.logger.info("⚡ Computing B-grid ispd")
            ispd_B       = xr.apply_ufunc(np.hypot, DS['uvel'], DS['vvel'],
                                          dask          = "allowed", 
                                          output_dtypes = [DS['uvel'].dtype])
            DS['ispd_B'] = ispd_B
        if "ispd_Ta" in ispd_type:
            self.logger.info("🌀 Spatial average regrid to T-grid")
            ispd_Ta                = xr.full_like(ispd_B, fill_value=np.nan)
            ispd_stack             = xr.concat([ispd_B.isel(nj=slice(0, -1), ni=slice(0, -1)),
                                                ispd_B.isel(nj=slice(0, -1), ni=slice(1, None)),
                                                ispd_B.isel(nj=slice(1, None), ni=slice(0, -1)),
                                                ispd_B.isel(nj=slice(1, None), ni=slice(1, None)) ], dim="offset")
            ispd_avg               = ispd_stack.mean(dim="offset", skipna=True)    
            ispd_Ta[..., :-1, :-1] = ispd_avg.data  # or .values
            DS['ispd_Ta']          = xr.DataArray(ispd_Ta,
                                                dims   = ("time", "nj", "ni"),
                                                coords = {"time" : (("time"),DS.time.values),
                                                            "nj"   : np.arange(ispd_Ta.sizes["nj"]),
                                                            "ni"   : np.arange(ispd_Ta.sizes["ni"]),
                                                            "TLAT" : (("nj", "ni"), DS.TLAT.values),
                                                            "TLON" : (("nj", "ni"), DS.TLON.values)},
                                                attrs  = {"long_name" : "T-grid interpolated B-grid ice speed",
                                                            "units"     : ispd_B.attrs.get("units", "m/s")})
        if "ispd_Tx" in ispd_type:
            self.logger.info("🌀 xESMF regrid to T-grid")
            if not self.reG_weights_defined:
                self.define_reG_weights()
            DS_reG = self.reG_bgrid_to_tgrid_xesmf(DS)
            self.logger.info("⚡ Computing xESMF-based ispd")
            ispd_Tx         = xr.apply_ufunc(np.hypot, DS_reG['uvel'], CICE_reG['vvel'],
                                            dask          = "allowed",
                                            output_dtypes = [DS_reG['uvel'].dtype])
            DS['ispd_Tx'] = ispd_Tx.copy()
        return DS

    def process_daily_cice(self, sim_name=None, dt0_str=None, dtN_str=None, D_iceh=None, ispd_thresh=None, ispd_type=None):
        sim_name    = sim_name    or self.sim_name
        dt0_str     = dt0_str     or self.dt0_str
        dtN_str     = dtN_str     or self.dtN_str
        D_iceh      = D_iceh      or self.D_iceh
        ispd_thresh = float(ispd_thresh) if ispd_thresh is not None else self.ispd_thresh
        ispd_type   = ispd_type   or ["ispd_B", "ispd_Ta", "ispd_Tx"]
        if isinstance(ispd_type, str):
            ispd_type = [ispd_type]
        valid_types = {"ispd_B", "ispd_Ta", "ispd_Tx"}
        assert all(t in valid_types for t in ispd_type), f"Invalid ispd_type: {ispd_type}"
        dts = pd.date_range(dt0_str, dtN_str, freq="D")                          
        for dt in dts:
            date_str = dt.strftime("%Y-%m-%d")
            m_str    = dt.strftime("%Y-%m")
            P_iceh   = Path(D_iceh, f"iceh.{date_str}.nc")
            if not P_iceh.exists():
                self.logger.info(f"❌ Missing CICE file: {P_iceh}")
                continue
            self.logger.info(f"📂 Loading: {P_iceh}")
            CICE      = xr.open_dataset(P_iceh, engine="netcdf4")[list(self.cice_var_list)]
            CICE_ispd = self.compute_ice_speed_types( CICE , ispd_type )
            # if "ispd_B" in ispd_type:
            #     self.logger.info("⚡ Computing B-grid ispd")
            #     ispd_B         = xr.apply_ufunc(np.hypot, CICE['uvel'], CICE['vvel'],
            #                                     dask          = "allowed", 
            #                                     output_dtypes = [CICE['uvel'].dtype])
            #     CICE['ispd_B'] = ispd_B
            # if "ispd_Ta" in ispd_type:
            #     self.logger.info("🌀 Spatial average regrid to T-grid")
            #     ispd_Ta                = xr.full_like(ispd_B, fill_value=np.nan)
            #     ispd_stack             = xr.concat([ispd_B.isel(nj=slice(0, -1), ni=slice(0, -1)),
            #                                         ispd_B.isel(nj=slice(0, -1), ni=slice(1, None)),
            #                                         ispd_B.isel(nj=slice(1, None), ni=slice(0, -1)),
            #                                         ispd_B.isel(nj=slice(1, None), ni=slice(1, None)) ], dim="offset")
            #     ispd_avg               = ispd_stack.mean(dim="offset", skipna=True)    
            #     ispd_Ta[..., :-1, :-1] = ispd_avg.data  # or .values
            #     CICE['ispd_Ta']        = xr.DataArray(ispd_Ta,
            #                                           dims   = ("time", "nj", "ni"),
            #                                           coords = {"time" : (("time"),CICE.time.values),
            #                                                     "nj"   : np.arange(ispd_Ta.sizes["nj"]),
            #                                                     "ni"   : np.arange(ispd_Ta.sizes["ni"]),
            #                                                     "TLAT" : (("nj", "ni"), CICE.TLAT.values),
            #                                                     "TLON" : (("nj", "ni"), CICE.TLON.values)},
            #                                          attrs  = {"long_name" : "T-grid interpolated B-grid ice speed",
            #                                                     "units"     : ispd_B.attrs.get("units", "m/s")})
            # if "ispd_Tx" in ispd_type:
            #     self.logger.info("🌀 xESMF regrid to T-grid")
            #     if not self.reG_weights_defined:
            #         self.define_reG_weights()
            #     CICE_reG = self.reG_bgrid_to_tgrid_xesmf(CICE)
            #     self.logger.info("⚡ Computing xESMF-based ispd")
            #     ispd_Tx         = xr.apply_ufunc(np.hypot, CICE_reG['uvel'], CICE_reG['vvel'],
            #                                      dask          = "allowed",
            #                                      output_dtypes = [CICE_reG['uvel'].dtype])
            #     CICE['ispd_Tx'] = ispd_Tx.copy()
            self.logger.info("🌊 Subsetting Ocean into either southern or northern hemisphere (default: southern)")
            CICE_SO = CICE_ispd.isel(nj=self.hemisphere_nj_slice)
            self.logger.info("❄️ create fast ice masks")
            sic_mask = CICE_SO['aice'] > self.icon_thresh
            masks    = {}
            if "ispd_B" in ispd_type:
                masks['FI_B']  = sic_mask & (CICE_SO['ispd_B']  <= ispd_thresh)
            if "ispd_Ta" in ispd_type:
                masks['FI_Ta'] = sic_mask & (CICE_SO['ispd_Ta'] <= ispd_thresh)
            if "ispd_Tx" in ispd_type:
                masks['FI_Tx'] = sic_mask & (CICE_SO['ispd_Tx'] <= ispd_thresh)
            self.logger.info("❄️ apply fast ice masks to dataset")
            for key in masks:
                mask = masks[key]
                ds_fi = CICE_SO.where(mask)
                ds_fi['FI_mask'] = mask
                self.datasets_by_group[key][m_str].append(ds_fi)
                ds_pi = CICE_SO.where(~mask)
                ds_pi['PI_mask'] = ~mask
                self.datasets_by_group[key.replace("FI", "PI")][m_str].append(ds_pi)
            self.datasets_by_group['SO'][m_str].append(CICE_SO)
            if dt == dts[-1] or dt.month != (dt + pd.Timedelta(days=1)).month:
                ispd_thresh_str = f"{ispd_thresh:.1e}".replace("e-0", "e-") 
                P_zarr_root     = Path(self.D_zarr, f"ispd_thresh_{ispd_thresh_str}", f"cice_daily_{m_str}.zarr")
                for group, data_dict in self.datasets_by_group.items():
                    datasets = data_dict[m_str]
                    if not datasets:
                        continue
                    self.logger.info(f"💾 Writing group {group} to {P_zarr_root}")
                    ds_monthly = xr.concat(datasets, dim="time").sortby("time")
                    if group != "SO":
                        ds_monthly.attrs["ispd_thresh"] = ispd_thresh
                        ds_monthly.attrs["mask_type"]   = self.mask_type_map.get(group, "unknown")
                    ds_monthly.to_zarr(P_zarr_root, group=group, mode="w", consolidated=True)
                    data_dict[m_str].clear()

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

    def process_rolling_cice(self,
                             sim_name       = None,
                             dt0_str        = None,
                             dtN_str        = None,
                             D_iceh         = None,
                             mean_period    = None,
                             ispd_thresh    = None,
                             ispd_type      = None):
        sim_name     = sim_name or self.sim_name
        dt0_str      = dt0_str or self.dt0_str
        dtN_str      = dtN_str or self.dtN_str
        D_iceh       = D_iceh or self.D_iceh
        mean_period  = mean_period or self.mean_period
        ispd_thresh  = float(ispd_thresh) if ispd_thresh is not None else self.ispd_thresh
        ispd_type    = ispd_type or ["ispd_B", "ispd_Ta", "ispd_Tx"]
        if isinstance(ispd_type, str):
            ispd_type = [ispd_type]
        valid_types = {"ispd_B", "ispd_Ta", "ispd_Tx"}
        assert all(t in valid_types for t in ispd_type), f"Invalid ispd_type: {ispd_type}"
        dts = pd.date_range(dt0_str, dtN_str, freq="D")
        month_groups = defaultdict(lambda: {k: [] for k in ['FI_B', 'FI_Ta', 'FI_Tx', 'PI_B', 'PI_Ta', 'PI_Tx', 'SO']})
        self.logger.info(f"📆 Rolling mean first, then fast ice masking between {dt0_str} and {dtN_str}")
        self.logger.info("📚 Loading entire daily dataset with variables: %s", self.cice_var_list)
        try:
            files = self.get_cice_files_between_dates(D_iceh, dt0_str, dtN_str)
            self.logger.info(f"📦 Loading {len(files)} CICE history files from {files[0].name} to {files[-1].name}")
            def _preprocess_cice(ds):
                return ds[list(self.cice_var_list)]
            CICE_all = xr.open_mfdataset(files, engine="netcdf4", parallel=True, preprocess=_preprocess_cice, combine="by_coords")
            CICE_all = CICE_all.chunk({'time': mean_period * 2})
        except Exception as e:
            self.logger.error(f"Failed to load datasets: {e}")
            return
        self.logger.info("🌀 Computing rolling means")
        CICE_roll = CICE_all.rolling(time=mean_period, center=True, min_periods=1).mean().compute()
        for dt in dts:
            date_str = dt.strftime("%Y-%m-%d")
            m_str = dt.strftime("%Y-%m")
            try:
                CICE_day = CICE_roll.sel(time=dt, method="nearest")
            except Exception as e:
                self.logger.warning(f"⛔ Could not find rolling average data for {date_str}: {e}")
                continue
            self.logger.info(f"📆 {date_str}: computing speeds and masks")
            masks = {}
            # if "ispd_B" in ispd_type or "ispd_Ta" in ispd_type:
            #     ispd_B = xr.apply_ufunc(np.hypot, CICE_day['uvel'], CICE_day['vvel'])
            #     CICE_day['ispd_B'] = ispd_B
            # if "ispd_Ta" in ispd_type:
            #     self.logger.info("🌀 Spatial average regrid to T-grid")
            #     ispd_Ta                = xr.full_like(ispd_B, fill_value=np.nan)
            #     ispd_stack             = xr.concat([ispd_B.isel(nj=slice(0, -1), ni=slice(0, -1)),
            #                                         ispd_B.isel(nj=slice(0, -1), ni=slice(1, None)),
            #                                         ispd_B.isel(nj=slice(1, None), ni=slice(0, -1)),
            #                                         ispd_B.isel(nj=slice(1, None), ni=slice(1, None)) ], dim="offset")
            #     ispd_avg               = ispd_stack.mean(dim="offset", skipna=True)    
            #     ispd_Ta[..., :-1, :-1] = ispd_avg.data  # or .values
            #     CICE_day['ispd_Ta']    = xr.DataArray(ispd_Ta,
            #                                           dims   = ("time", "nj", "ni"),
            #                                           coords = {"time" : (("time"),CICE_day.time.values),
            #                                                     "nj"   : np.arange(ispd_Ta.sizes["nj"]),
            #                                                     "ni"   : np.arange(ispd_Ta.sizes["ni"]),
            #                                                     "TLAT" : (("nj", "ni"), CICE_day.TLAT.values),
            #                                                     "TLON" : (("nj", "ni"), CICE_day.TLON.values)},
            #                                          attrs  = {"long_name" : "T-grid interpolated B-grid ice speed",
            #                                                     "units"     : ispd_B.attrs.get("units", "m/s")})                
            # if "ispd_Tx" in ispd_type:
            #     if not self.reG_weights_defined:
            #         self.define_reG_weights()
            #     try:
            #         CICE_reG = self.reG_bgrid_to_tgrid_xesmf(CICE_day)
            #         ispd_Tx = xr.apply_ufunc(np.hypot, CICE_reG['uvel'], CICE_reG['vvel'])
            #         CICE_day['ispd_Tx'] = ispd_Tx
            #     except Exception as e:
            #         self.logger.error(f"❌ Regridding failed: {e}")
            #         continue
            CICE_ispd = self.compute_ice_speed_types( CICE_day , ispd_type )
            CICE_SO   = CICE_ispd.isel(nj=self.hemisphere_nj_slice)
            sic_mask  = CICE_SO['aice'] > self.icon_thresh
            if "ispd_B" in ispd_type:
                masks['FI_B'] = sic_mask & (CICE_SO['ispd_B'] <= ispd_thresh)
            if "ispd_Ta" in ispd_type:
                masks['FI_Ta'] = sic_mask & (CICE_SO['ispd_Ta'] <= ispd_thresh)
            if "ispd_Tx" in ispd_type:
                masks['FI_Tx'] = sic_mask & (CICE_SO['ispd_Tx'] <= ispd_thresh)
            for group in masks:
                fi_mask = masks[group]
                pi_mask = ~fi_mask
                ds_fi = CICE_SO.where(fi_mask)
                ds_pi = CICE_SO.where(pi_mask)
                ds_fi['FI_mask'] = fi_mask
                ds_pi['PI_mask'] = pi_mask
                pi_group = group.replace("FI_", "PI_")
                month_groups[m_str][group].append(ds_fi)
                month_groups[m_str][pi_group].append(ds_pi)
            month_groups[m_str]['SO'].append(CICE_SO)
            if dt == dts[-1] or dt.month != (dt + pd.Timedelta(days=1)).month:
                ispd_thresh_str = f"{ispd_thresh:.1e}".replace("e-0", "e-")
                P_zarr_root = Path(self.D_zarr, f"ispd_thresh_{ispd_thresh_str}", f"cice_rolling_{m_str}.zarr")
                for group, datasets in month_groups[m_str].items():
                    if group != "SO" and group.replace("FI_", "ispd_").replace("PI_", "ispd_") not in ispd_type:
                        continue
                    if not datasets:
                        continue
                    self.logger.info(f"💾 Writing group {group} to {P_zarr_root}")
                    try:
                        ds_monthly = xr.concat(datasets, dim="time").sortby("time")
                        if group != "SO":
                            ds_monthly.attrs["ispd_thresh"] = ispd_thresh
                            ds_monthly.attrs["mask_type"] = self.mask_type_map.get(group, "unknown")
                        ds_monthly.to_zarr(P_zarr_root, group=group, mode="w", consolidated=True)
                    except Exception as e:
                        self.logger.error(f"Failed to write group {group}: {e}")
                month_groups[m_str].clear()

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
        assert ice_type in ['FI_B', 'FI_Ta', 'FI_Tx', 'PI_B', 'PI_Ta', 'PI_Tx', 'SO'], f"Invalid ice_type: {ice_type}"
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
        IA = (((SIC * GC_area).sum(dim=("nj", "ni"))/ice_area_scale)+self.GI_total_area/1e4).persist()
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

    def define_AF2020_reG_weights(self):
        from pyproj import CRS, Transformer
        self.logger.info("defining AF2020 regridder weights")
        F_weights        = self.sea_ice_dict["AF_reG_weights"]
        self.G_u         = self.build_grid_dict(self.gi_processor.G_u['lat'], self.gi_processor.G_u['lon'])
        self.G_u["mask"] = self.gi_processor.KMT_mod
        self.logger.info("convert 'AF_FI_OBS_2020db' Cartesian coordinates to spherical coordindates")
        crs_obs          = CRS.from_epsg(self.sea_ice_dict["projection_FI_obs"]) #unique to observations
        crs_spherical    = CRS.from_epsg(self.sea_ice_dict["projection_wgs84"])  #spherical coordinates
        transformer      = Transformer.from_crs(crs_obs, crs_spherical, always_xy=True)
        X, Y             = np.meshgrid(FI_obs_native['x'].isel(time=0).values, FI_obs_native['y'].isel(time=0).values)
        lon_obs, lat_obs = transformer.transform(X,Y)
        G_obs            = self.build_grid_dict(lat_obs, lon_obs)
        self.logger.info("load/create regridder into memory")
        if self.sea_ice_dict['overwrite_weights'] or os.path.exists(F_weights) is None:
            self.reG_AF2020 = xe.Regridder(self.G_u, self.G_obs,
                                          method            = "bilinear",
                                          periodic          = True,
                                          ignore_degenerate = True,
                                          extrap_method     = "nearest_s2d",
                                          reuse_weights     = True,
                                          weights           = F_weights if os.path.exists(F_weights) else None,
                                          filename          = None if os.path.exists(F_weights) else F_weights)
        self.reG_AF2020_weights_defined = True

    def regrid_AF2020_to_ugrid(self, FI_obs_native):
        if not self.reG_AF2020_weights_defined:
            self.defin_AF2020_reG_weights()
        self.logger.info("*** Regridding 'AF_FI_OBS_2020db' to CICE U-grid *** ")
        self.FI_obs_reG_lon = self.reG_AF2020(self.G_obs["lon"]).compute()
        self.FI_obs_reG_lat = self.reG_AF2020(self.G_obs["lat"]).compute()
        self.FI_obs_reG_dat = self.reG_AF2020( FI_obs_native["Fast_Ice_Time_series"] ).compute()

    def load_AF2020_FI_org_netcdf(self, P_orgs):
        FI_obs = xr.open_mfdataset(P_orgs, engine='netcdf4')
        self.regrid_AF2020_to_ugrid(FI_obs)
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
