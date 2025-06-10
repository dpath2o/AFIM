import os, sys, time, json, imageio, shutil, pygmt, imageio, shutil, re, zarr, logging, time
import xarray             as xr
import pandas             as pd
import numpy              as np
import geopandas          as gpd
import matplotlib.pyplot  as plt
import matplotlib.dates   as mdates
from pathlib              import Path
from datetime             import datetime, timedelta
from dask.distributed     import Client, LocalCluster
from collections          import defaultdict
from pathlib              import Path
sys.path.insert(0, '/home/581/da1339/AFIM/src/AFIM/src')
_dask_client = None

class SeaIceScience:
    """

    A class for processing and analysing Antarctic sea ice output from CICE simulations.

    This processor is tailored to support the Antarctic Fast Ice Model (AFIM) workflow, enabling
    high-resolution diagnostics of landfast and pack ice using sea ice concentration, speed,
    and stress fields. It provides methods to generate fast ice masks, persistence fields,
    climatological metrics, and observational comparisons via Fraser et al. (2020) and NSIDC data.

    Core capabilities include:
    + Rolling and boolean-based fast ice classification
    + Fast and pack ice area computation (FIA, PIA)
    + Fast ice persistence metrics (FIP) and climatological summaries
    + Zarr-based storage of processed outputs and metric summaries
    + Support for automated regional map generation and time series plots

    This class is tightly coupled with:
    + `GroundedIcebergProcessor` for sub-grid grounded icebergs and modified landmasking
    + Observational fast ice climatology from Fraser et al. 2020: https://doi.org/10.5194/essd-12-2987-2020
    + NSIDC sea ice concentration data (optional, for observed pack ice comparisons)
    + A JSON configuration file defining simulation paths and parameter settings:
        https://github.com/dpath2o/AFIM/blob/main/src/AFIM/src/JSONs/afim_cice_analysis.json

    Parameters are inferred from an AFIM configuration file and simulation-specific keys,
    including ice speed thresholds, rolling window size, hemispheric focus, output directories,
    and optional data thinning or landmask overrides.

    This class supports both batch-style PBS workflows and interactive Jupyter sessions.

    ATTRIBUTES:
       sim_name            : str; name of the CICE simulation to analyze.
       dt0_str, dtN_str    : str; start and end dates of the analysis period (YYYY-MM-DD).
       icon_thresh         : float; sea ice concentration threshold for landfast ice presence.
       ispd_thresh         : float; ice speed threshold (in m/s) below which ice is considered fast.
       mean_period         : int; window (in days) used for temporal rolling means.
       bool_window         : int; window size for boolean fast ice presence filter.
       bool_min_days       : int; minimum number of fast ice days in window to qualify as persistent.
       GI_proc             : GroundedIcebergProcessor; object handling grounded iceberg identification and masking.
       D_iceh              : Path; directory containing raw daily CICE ice history NetCDF files.
       D_zarr              : Path; output directory for all derived Zarr datasets.
       valid_ispd_types    : list of str; ccepted ice speed types ("ispd_B", "ispd_Ta", "ispd_Tx", "ispd_BT").
       valid_zarr_datasets : list of str; output dataset types written to Zarr ("FI_B", etc.).
       hemisphere_nj_slice : slice; latitudinal range to use for hemisphere-specific masking (e.g., Southern Ocean).

    NOTES:
    + Fast ice masks are computed using ice concentration and speed thresholds and may optionally
      incorporate grounded iceberg presence and landmask modifications.
    + Outputs can be used for model sensitivity studies, climatology diagnostics, or comparisons to observations.
    + Rolling and boolean methods align with definitions in Fraser et al. (2020, ESSD).
    + NSIDC sea ice metrics (SIA, SIE, SIC) can be computed in coordination with this class via PackIceProcessor.

    EXAMPLES:
    >>> SI_proc = SeaIceProcessor(sim_name="gi-mid")
    >>> FI = SI_proc.process_daily_cice(ispd_thresh=1e-3, ispd_type="ispd_BT")
    >>> FI_roll = SI_proc.process_rolling_cice(ispd_thresh=1e-3, ispd_type="ispd_BT", mean_period=15)
    >>> SI_proc.sea_ice_metrics_wrapper(ice_type="FI_BT", ispd_thresh=1e-3)

    For more knowledge go to AFIM repository: https://github.com/dpath2o/AFIM

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
	             overwrite_AF2020_zarr       = False,
                 overwrite_saved_figs        = False,
                 save_new_figs               = False,
                 show_figs                   = False):
        self.sim_name = sim_name
        if P_json is None:
            P_json = '/home/581/da1339/AFIM/src/AFIM/src/JSONs/afim_cice_analysis.json'
        with open(P_json, 'r') as f:
            self.config = json.load(f)
        self.CICE_dict      = self.config.get("CICE_dict" , {})
        self.GI_dict        = self.config.get('GI_dict'   , {})
        self.D_dict         = self.config.get('D_dict'    , {})
        self.NSIDC_dict     = self.config.get('NSIDC_dict', {})
        D_log               = self.D_dict['logs']
        P_log               = P_log                       if P_log                       is not None else Path(D_log, f'SeaIceProcessor_FI_{self.sim_name}.log')
        self.dt0_str        = dt0_str                     if dt0_str                     is not None else self.config.get('dt0_str', '1993-01-01')
        self.dtN_str        = dtN_str                     if dtN_str                     is not None else self.config.get('dtN_str', '1999-12-31')
        self.mean_period    = mean_period                 if mean_period                 is not None else self.config.get('mean_period', 15)
        self.bool_window    = boolean_window              if boolean_window              is not None else self.config.get('bool_window', 7)
        self.bool_min_days  = boolean_min_days            if boolean_min_days            is not None else self.config.get('bool_min_days', 6)
        self.ispd_thresh    = ice_speed_threshold         if ice_speed_threshold         is not None else self.config.get('ice_speed_thresh_hi', 1e-3)
        self.icon_thresh    = ice_concentration_threshold if ice_concentration_threshold is not None else self.config.get('ice_conc_thresh', 0.15)
        self.cice_vars_ext  = extra_cice_vars             if extra_cice_vars             is not None else self.CICE_dict["FI_cice_vars_ext"]
        self.ow_fig         = overwrite_saved_figs        if overwrite_saved_figs        is not None else False
        self.save_fig       = save_new_figs               if save_new_figs               is not None else False
        self.show_fig       = show_figs                   if show_figs                    is not None else True
        hemisphere          = hemisphere                  if hemisphere                  is not None else self.config.get('hemisphere', 'south')
        if not os.path.exists(P_log):
            os.system(f"touch {P_log}")
        self.setup_logging(logfile=P_log)
        self.define_hemisphere(hemisphere)
        self.valid_ispd_types    = ["ispd_B", "ispd_Ta", "ispd_Tx", "ispd_BT"]
        self.valid_zarr_datasets = ['FI_B', 'FI_Ta', 'FI_Tx', 'FI_BT']
        self.cice_vars_reqd      = self.CICE_dict["FI_cice_vars_reqd"]
        self.cice_var_list       = self.cice_vars_reqd + self.cice_vars_ext
        self.D_sim               = Path(self.D_dict['AFIM_out'], sim_name)
        self.D_iceh              = Path(self.D_sim , 'history', 'daily')
        self.D_zarr              = Path(self.D_sim , 'zarr')
        self.D_graph             = Path(self.config['D_dict']['graph'], 'AFIM')
        self.doy_vals            = self.config.get("DOY_vals", list(range(1, 366, 15)))
        self.sea_ice_dict        = self.config.get("sea_ice_dict", {})
        self.pygmt_dict          = self.config.get("pygmt_dict", {})
        self.plot_var_dict       = self.config.get("plot_var_dict", {})
        self.reg_dict            = self.config.get('AF_regions', {})
        self.FIC_scale           = self.config.get('FIC_scale', 1e9)
        self.SIC_scale           = self.config.get('SIC_scale', 1e12)
        self.cm2m_fact           = self.config.get('cm2m_fact', 0.01)
        self.P_KMT_org           = Path(self.GI_dict["D_GI_thin"],self.GI_dict['KMT_org_fmt'])
        self.sim_config          = self.parse_simulation_metadata()
        self.GI_thin_fact        = self.sim_config.get('GI_thin_fact')
        self.GI_thin_str         = f"{self.GI_thin_fact:0.2f}".replace('.', 'p')
        self.GI_version          = self.sim_config.get('GI_version')
        self.GI_version_str      = f"{self.GI_version:.2f}".replace('.', 'p')
        if self.GI_thin_fact>0:
            self.use_gi    = True
            self.P_KMT_mod = os.path.join(self.GI_dict['D_GI_thin'],
                                          self.GI_dict['KMT_mod_fmt'].format(GI_thin_fact = self.GI_thin_str,
                                                                             version      = self.GI_version_str))
        else:
            self.use_gi = False
        self.ispd_thresh_metrics        = self.interpret_ice_speed_threshold(lat_thresh=-60)
        self.reG_weights_defined        = False
        self.reG_AF2020_weights_defined = False

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

    def parse_simulation_metadata(self):
        """
        Parse CICE diagnostic text file to extract numerical modeling parameters.

        Parameters
        ----------
        diag_path : str or Path
            Path to the diagnostics text file (e.g., ice_diag.d).

        Returns
        -------
        dict
            Dictionary of extracted parameters.
            Includes inferred GI thinning fraction if modified `kmt_file` is used.
        """
        P_diag     = Path(self.D_sim,"ice_diag.d")
        PARAM_KEYS = [ "dt", "ndtd", "ndte", "kdyn", "revised_evp", "e_yieldcurve", "e_plasticpot",
                       "Ktens", "kstrength", "Pstar", "Cstar", "Cf", "visc_method", "kmt_file" ]
        PATTERNS   = {key: re.compile(rf"{key}\s*=\s*(.+?)\s*:") if key != "kmt_file"
                      else re.compile(rf"{key}\s*=\s*(.+)$")
                      for key in PARAM_KEYS }
        result     = {key: "" for key in PARAM_KEYS}
        with open(P_diag, "r", encoding="utf-8", errors="replace") as f:
            for i, line in enumerate(f):
                if i > 500:
                    break
                for key, pattern in PATTERNS.items():
                    if result[key] == "":
                        match = pattern.search(line)
                        if match:
                            val = match.group(1).strip()
                            if key == "kmt_file":
                                val = Path(val).name
                            result[key] = val
        kmt_name = result["kmt_file"]
        if "kmt_mod_thinned-" in kmt_name:
            try:
                thin_str  = kmt_name.split("thinned-")[1].split("_")[0]  # e.g., "0p85"
                thin_frac = float(thin_str.replace("p", "."))
                result["GI_thin_fact"] = thin_frac #round(1.0 - thin_frac, 2)
            except Exception:
                result["GI_thin_fact"] = "?"
            try:
                version_str = kmt_name.split("_v")[1].split(".")[0]  # e.g., "1p50"
                version_float = float(version_str.replace("p", "."))
                result["GI_version"] = version_float
            except Exception:
                result["GI_version"] = "?"
        else:
            result["GI_thin_fact"] = 0.0
            result["GI_version"] = 0.0
        return result

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

    def interpret_ice_speed_threshold(self, ispd_thresh=None, lat_thresh=-60):
        ispd_thresh  = ispd_thresh if ispd_thresh is not None else self.ispd_thresh
        m_per_day    = ispd_thresh * 86400  # meters/day
        area_da      = xr.open_dataset( self.CICE_dict["P_G"] )['tarea']                         # [m^2]
        lat_da       = self.radians_to_degrees( xr.open_dataset(self.CICE_dict["P_G"])['tlat'] ) # [degrees]
        mask         = lat_da < lat_thresh
        area_vals    = area_da.where(mask).values
        grid_lengths = np.sqrt(area_vals)
        grid_lengths = grid_lengths[np.isfinite(grid_lengths)]
        if len(grid_lengths) == 0:
            raise ValueError("No valid grid cells found south of the specified latitude.")
        GC_len_median = np.median(grid_lengths)
        pct_GC_disp   = m_per_day / GC_len_median
        days_per_GC   = GC_len_median / m_per_day
        self.logger.info(f"Ice speed threshold                        : {ispd_thresh:.1e} m/s → {m_per_day:.1f} m/day")
        self.logger.info(f"Median grid cell length below {lat_thresh}°: {GC_len_median:.1f} m")
        self.logger.info(f"→ Displacement                             = {pct_GC_disp*100:.2f}% of grid cell per day")
        self.logger.info(f"→ Days to fully traverse one grid cell     : {days_per_GC:.2f} days")
        return {"ice_speed_thresh_m_per_s"    : ispd_thresh,
                "displacement_m_per_day"      : m_per_day,
                "median_grid_cell_length_m"   : GC_len_median,
                "percent_displacement_per_day": pct_GC_disp,
                "days_per_grid_cell"          : days_per_GC}

    def get_dir_size(self, path):
        size_gb = sum(f.stat().st_size for f in path.rglob("*") if f.is_file()) / (1024**3)
        self.logger.info(f"Disk-usage (size) of directory {path}: {size_gb:.2f} GB")

    def count_zarr_files(self, path):
        total_files = sum(len(files) for _, _, files in os.walk(path))
        self.logger.info(f"{path} contains {total_files} files")

    def radians_to_degrees(self, da):
        return (da * 180) / np.pi

    def normalise_longitudes(self,lon):
        return ((lon + 360) % 360) - 180

    def build_grid_corners(self, lat_rads, lon_rads, grid_res=0.25):
        lon                  = self.radians_to_degrees(lon_rads)
        lat                  = self.radians_to_degrees(lat_rads)
        ny, nx               = lat.shape
        lat_b                = np.zeros((ny + 1, nx + 1))
        lon_b                = np.zeros((ny + 1, nx + 1))
        lat_b[1:-1, 1:-1]    = grid_res * (lat[:-1, :-1] + lat[:-1, 1:] + lat[1:, :-1] + lat[1:, 1:])
        lon_b[1:-1, 1:-1]    = grid_res * (lon[:-1, :-1] + lon[:-1, 1:] + lon[1:, :-1] + lon[1:, 1:])
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
        lon_b                = self.normalise_longitudes(lon_b)
        return lon_b, lat_b

    def define_cice_grid(self, grid_type='t', mask=False, build_grid_corners=False):
        std_dim_names = ("nj", "ni")
        G             = xr.open_dataset(self.CICE_dict["P_G"])
        lon_rads      = G[f"{grid_type}lon"].values
        lat_rads      = G[f"{grid_type}lat"].values
        ny, nx        = lon_rads.shape
        lon           = self.normalise_longitudes(self.radians_to_degrees(lon_rads))
        lat           = self.radians_to_degrees(lat_rads)
        area          = G[f'{grid_type}area'].values
        data_vars     = {"area": (std_dim_names, area),
                         "lon":  (std_dim_names, lon),
                         "lat":  (std_dim_names, lat)}
        coords        = {"nj": np.arange(ny),
                         "ni": np.arange(nx)}
        if build_grid_corners:
            crn_dim_names      = ("nj_b", "ni_b")
            lon_b, lat_b       = self.build_grid_corners(lon_rads, lat_rads)
            data_vars["lon_b"] = (crn_dim_names, lon_b)
            data_vars["lat_b"] = (crn_dim_names, lat_b)
            coords["nj_b"]     = np.arange(ny + 1)
            coords["ni_b"]     = np.arange(nx + 1)
        if mask:
            mask_data = xr.open_dataset(self.P_KMT_mod if self.use_gi else self.P_KMT_org).kmt.values
            data_vars["mask"] = (std_dim_names, mask_data)
        return xr.Dataset(data_vars=data_vars, coords=coords)

    def define_reG_weights(self):
        G_u           = self.define_cice_grid( grid_type='u'             , build_grid_corners=True )
        G_t           = self.define_cice_grid( grid_type='t' , mask=True , build_grid_corners=True )
        F_weights     = self.CICE_dict["P_reG_u2t_weights"]
        weights_exist = os.path.exists(F_weights)
        self.logger.info(f"{'Reusing' if weights_exist else 'Creating'} regrid weights: {F_weights}")
        self.reG = xe.Regridder(G_u, G_t,
                                method            = "bilinear",
                                periodic          = True,
                                ignore_degenerate = True,
                                extrap_method     = "nearest_s2d",
                                reuse_weights     = weights_exist,
                                filename          = F_weights)
        self.reG_weights_defined = True

    def reG_bgrid_to_tgrid_xesmf(self, ds):
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
