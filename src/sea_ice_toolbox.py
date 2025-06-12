import os, sys, time, json, imageio, shutil, pygmt, imageio, shutil, re, zarr, logging, time
import xarray             as xr
import xesmf              as xe
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
from sea_ice_models       import SeaIceModels
from sea_ice_plotter      import SeaIcePlotter
from sea_ice_icebergs     import SeaIceIcebergs
from sea_ice_observations import SeaIceObservations
_dask_client = None

#class SeaIceToolbox:
class SeaIceToolbox(SeaIceModels, SeaIcePlotter, SeaIceIcebergs, SeaIceObservations):
    """
    Unified toolbox for processing and analysing Antarctic sea ice from CICE simulations
    as part of the Antarctic Fast Ice Model (AFIM) workflow.

    This superclass consolidates functionality from:
    + SeaIceModels: for simulation I/O and fast/pack ice masking
    + SeaIcePlotter: for PyGMT and timeseries visualization
    + SeaIceIcebergs: for grounded iceberg thinning/masking
    + SeaIceObservations: for Fraser et al. (2020) and NSIDC comparisons

    + Observational fast ice climatology from Fraser et al. 2020: https://doi.org/10.5194/essd-12-2987-2020
    + NSIDC sea ice concentration data (optional, for observed pack ice comparisons)
    + A JSON configuration file defining simulation paths and parameter settings:
        https://github.com/dpath2o/AFIM/blob/main/src/AFIM/src/JSONs/afim_cice_analysis.json    

    Core capabilities:
    + Fast and pack ice area computation (FIA, PIA)
    + Rolling and boolean-based fast ice classification
    + Fast ice persistence and climatological metrics (FIP)
    + Regional map generation and timeseries analysis
    + Zarr-based output, with optional observational overlays

    External dependencies:
    - GroundedIcebergProcessor for KMT and GI masks
    - Fraser et al. (2020) observational fast ice climatology
    - NSIDC daily sea ice concentration fields (v4)
    - Configured via a user-defined JSON file

    Parameters
    ----------
    P_json : str or Path, optional
        Path to AFIM JSON configuration file. Defaults to project root if None.
    sim_name : str, required
        Simulation name (must match entry under `AFIM_out` in config).
    dt0_str, dtN_str : str, optional
        Analysis start and end dates (YYYY-MM-DD).
    ice_concentration_threshold : float, optional
        Concentration cutoff for fast ice. Default: 0.15
    ice_speed_threshold : float, optional
        Speed threshold (m/s) below which ice is considered fast. Default: 5e-4
    ice_speed_type : str or list of str, optional
        Valid speed types: 'ispd_B', 'ispd_Ta', 'ispd_Tx', 'ispd_BT'.
    ice_type : str or list of str, optional
        Ice classification scheme (e.g., 'FI_BT', 'FI_BT_bool').
    mean_period : int, optional
        N-day rolling window for averaging (default: 15).
    boolean_window : int, optional
        Window size for boolean filtering (default: 7).
    boolean_min_days : int, optional
        Minimum required fast ice days in window (default: 6).
    extra_cice_vars : list of str, optional
        Additional CICE variables required for output (default: config-defined).
    hemisphere : str, optional
        'south' or 'north'. Controls region masks, slicing.
    P_log : Path, optional
        File to write runtime logs to.
    overwrite_zarr : bool, optional
        Whether to overwrite any Zarr file group(s).
    overwrite_AF2020_zarr : bool, optional
        Whether to overwrite regridded AF2020 dataset.
    overwrite_saved_figs : bool, optional
        Overwrite previously saved figures.
    save_new_figs : bool, optional
        Whether to save new figures to disk.
    show_figs : bool, optional
        Whether to plot figures interactively.

    Usage
    -----
    >>> SI_tools = SeaIceToolbox(sim_name='elps-min', dt0_str='1993-01-01', dtN_str='1999-12-31')
    >>> FI = SI_tools.process_daily_cice(ispd_thresh=5e-4)
    >>> SI_tools.plot_ice_area(FI)

    Repository:
    https://github.com/dpath2o/AFIM

    See also
    --------
    SeaIceModels, SeaIcePlotter, SeaIceIcebergs, SeaIceObservations
    """
    def summary(self):
        """Print a summary of key simulation setup and configuration."""
        print("--- SeaIceToolbox Summary ---")
        print(f"Simulation Name     : {self.sim_name}")
        print(f"Analysis Start Date : {self.dt0_str}")
        print(f"Analysis End Date   : {self.dtN_str}")
        print(f"Speed Threshold     : {self.ispd_thresh:.1e} m/s")
        print(f"Speed Type(s)       : {self.ispd_type}")
        print(f"Ice Type(s)         : {self.ice_type}")
        print(f"Mean Period         : {self.mean_period} days")
        print(f"Bool Window         : {self.bool_window} days")
        print(f"Bool Min Days       : {self.bool_min_days}")
        print(f"Using GI?           : {self.use_gi}")
        print(f"Overwrite Zarr      : {self.overwrite_zarr_group}")
        print(f"Save Figures        : {self.save_fig}")
        print(f"Show Figures        : {self.show_fig}")
        print(f"Hemisphere          : {self.hemisphere}")
        print("------------------------------")

    def help(self):
        """Print an overview of available methods grouped by module."""
        print("Available methods by submodule:")
        print("- SeaIceModels:")
        print("    • process_daily_cice")
        print("    • process_rolling_cice")
        print("    • compute_fast_ice_outputs")
        print("    • compute_pack_ice_outputs")
        print("    • sea_ice_metrics_wrapper")
        print("- SeaIcePlotter:")
        print("    • plot_ice_area")
        print("    • plot_timeseries_compare")
        print("    • plot_map")
        print("    • pygmt_map_plot_one_var")
        print("- SeaIceIcebergs:")
        print("    • load_GI_lon_lats")
        print("    • thin_grounded_icebergs")
        print("    • modify_landmask")
        print("- SeaIceObservations:")
        print("    • process_nsidc_obs")
        print("    • filter_AF2020_FI_by_date")
        print("    • regrid_AF2020_to_model")
        print("------------------------------")

    def __init__(self,
                 P_json                      = None,# the configuration file for which there are many dependencies
                                                    # that this toolbox relies upon
                 sim_name                    = None,# valid name of a model simulation; essentially 'valid' means
                                                    # any name given underneath the directory in the config file
                                                    # named 'AFIM_out'; the sea_ice_model class underneath the hood
                                                    # of this super-class relies on this name for processing and
                                                    # loading of simulation data
                 dt0_str                     = None,# the start period over which many methods underneath use
                                                    # format is YYYY-MM-DD; default is 1993-01-01
                 dtN_str                     = None,# the end period over which many methods underneath use
                                                    # format is YYYY-MM-DD; default is 1999-12-31
	             ice_concentration_threshold = None,# almost should never be changed from default value of
                                                    # 0.15 (grid cell concentration)
	             ice_speed_threshold         = None,# a significantly important value in the determination
                                                    # and classification (masking) of fast ice; defualt value
                                                    # 5e-4 m/s
                 ice_speed_type              = None,# a valid ice_speed_type or list thereof
                 ice_type                    = None,# a valid ice_type or list thereof
	             mean_period                 = None,# rolling average, N-days
                 boolean_window              = None,# the window with which to apply the minimum number of days
                 boolean_min_days            = None,# minimum number of days binary-days
	             extra_cice_vars             = None,# these will be included in the fast ice mask
                                                    # that is, in addtion to those listed in the
                                                    # in config file 'cice_vars_reqd'; default is
                                                    # also a list in config file 'FI_cice_vars_ext'
	             hemisphere                  = None,# used in many ares of the toolbox to define the hemisphere
                                                    # that the user is interested; unfortunately, the toolbox
                                                    # does not allow for a user to be interested in both at the 
                                                    # same time; either 'south' or 'north'; defualt is 'south'
	             P_log                       = None,# the log file to send print statements to	             
                 overwrite_zarr              = None,# whether or not to overwrite a zarr; default is false
	             overwrite_AF2020_zarr       = None,# whether or not to overwrite AF2020db zarr; default is false
                 overwrite_saved_figs        = None,# whether or not to overwite saved figures; default is false
                 save_new_figs               = None,# whether or not to write new figures to disk; default is true
                 show_figs                   = None,# whether or not to show/print figures to screen; default is false
                 delete_original_cice_iceh_nc= None,# whether or not to delete the original CICE ice history 
                 **kwargs):
        SeaIceModels.__init__(self, sim_name, **kwargs)
        SeaIcePlotter.__init__(self, **kwargs)
        SeaIceIcebergs.__init__(self, **kwargs)
        SeaIceObservations.__init__(self, **kwargs)
        self.sim_name = sim_name if sim_name is not None else 'test'
        if P_json is None:
            P_json = '/home/581/da1339/AFIM/src/AFIM/src/JSONs/sea_ice_config.json'
        with open(P_json, 'r') as f:
            self.config = json.load(f)
        self.CICE_dict          = self.config.get("CICE_dict" , {})
        self.GI_dict            = self.config.get('GI_dict'   , {})
        self.D_dict             = self.config.get('D_dict'    , {})
        self.NSIDC_dict         = self.config.get('NSIDC_dict', {}) 
        self.AF_FI_dict         = self.config.get("AF_FI_dict", {})
        self.Sea_Ice_Obs_dict   = self.config.get("Sea_Ice_Obs_dict", {})
        self.AOM2_dict          = self.config.get("AOM2_dict", {})
        self.MOM_dict           = self.config.get("MOM_dict", {})
        self.ERA5_dict          = self.config.get("ERA5_dict", {})
        self.ORAS_dict          = self.config.get("ORAS_dict", {})
        self.pygmt_dict         = self.config.get("pygmt_dict", {})
        self.plot_var_dict      = self.config.get("plot_var_dict", {})
        self.hemispheres_dict   = self.config.get("hemispheres_dict", {})
        self.Ant_8sectors       = self.config.get('Ant_8sectors', {})
        self.Ant_2sectors       = self.config.get('Ant_2sectors', {})
        self.plot_ice_area_dict = self.config.get('plot_ice_area_dict', {})
        D_log                   = self.D_dict['logs']
        P_log                   = P_log if P_log is not None else Path(D_log, f'SeaIceProcessor_FI_{self.sim_name}.log')
        if not os.path.exists(P_log):
            os.system(f"touch {P_log}")
        self.setup_logging(logfile=P_log)
        self.dt0_str              = dt0_str                     if dt0_str                     is not None else self.config.get('dt0_str', '1993-01-01')
        self.dtN_str              = dtN_str                     if dtN_str                     is not None else self.config.get('dtN_str', '1999-12-31')
        self.mean_period          = mean_period                 if mean_period                 is not None else self.config.get('mean_period', 15)
        self.bool_window          = boolean_window              if boolean_window              is not None else self.config.get('bool_window', 7)
        self.bool_min_days        = boolean_min_days            if boolean_min_days            is not None else self.config.get('bool_min_days', 6)
        self.ispd_thresh          = ice_speed_threshold         if ice_speed_threshold         is not None else self.config.get('ice_speed_thresh_hi', 1e-3)
        self.ispd_type            = ice_speed_type              if ice_speed_type              is not None else self.config.get('ice_speed_type', 'ispd_BT')
        self.ice_type             = ice_type                    if ice_type                    is not None else self.config.get('ice_type', 'FI_BT')
        self.icon_thresh          = ice_concentration_threshold if ice_concentration_threshold is not None else self.config.get('ice_conc_thresh', 0.15)
        self.cice_vars_ext        = extra_cice_vars             if extra_cice_vars             is not None else self.CICE_dict["FI_cice_vars_ext"]
        self.overwrite_zarr_group = overwrite_zarr              if overwrite_zarr              is not None else False
        self.ow_fig               = overwrite_saved_figs        if overwrite_saved_figs        is not None else False
        self.save_fig             = save_new_figs               if save_new_figs               is not None else True
        self.show_fig             = show_figs                   if show_figs                   is not None else False
        self.delete_original_cice_iceh_nc = delete_original_cice_iceh_nc if delete_original_cice_iceh_nc is not None else False
        hemisphere                = hemisphere                  if hemisphere                  is not None else self.config.get('hemisphere', 'south')
        self.define_hemisphere(hemisphere)
        self.ispd_thresh_str     = f"{self.ispd_thresh:.1e}".replace("e-0", "e-")
        self.D_sim               = Path(self.D_dict['AFIM_out'], sim_name)
        self.D_iceh              = Path(self.D_sim , 'history', 'daily')
        self.D_zarr              = Path(self.D_sim , 'zarr')
        self.D_graph             = Path(self.config['D_dict']['graph'], 'AFIM')
        self.D_metrics           = Path(self.D_zarr, f"ispd_thresh_{self.ispd_thresh_str}", "metrics")
        self.sim_config          = self.parse_simulation_metadata()   
        self.valid_ispd_types    = self.config.get("valid_ispd_types", [])
        self.valid_ice_types     = self.config.get("valid_ice_types", [])
        self.cice_vars_reqd      = self.CICE_dict["FI_cice_vars_reqd"]
        self.cice_var_list       = self.cice_vars_reqd + self.cice_vars_ext
        self.FIC_scale           = self.config.get('FIC_scale', 1e9)
        self.SIC_scale           = self.config.get('SIC_scale', 1e12)
        self.P_KMT_org           = Path(self.GI_dict["D_GI_thin"],self.GI_dict['KMT_org_fmt'])
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
        self.reG_weights_defined       = False
        self.modified_landmask_aligned = False
        self.bgrid_loaded              = False

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
            self.hemisphere_dict = self.hemispheres_dict['north']
        elif hemisphere.lower() in ['south', 'southern', 'sh', 's', 'so']:
            self.hemisphere_dict = self.hemispheres_dict['south']
        else:
            raise ValueError(f"Invalid hemisphere '{hemisphere}'. Valid options are: "
                              "['north', 'south', 'northern', 'southern', 'sh', 'nh', 'SH', 'NH']")
        self.hemisphere_dict['nj_slice'] = slice(self.hemisphere_dict['nj_slice'][0],self.hemisphere_dict['nj_slice'][1])
        self.hemisphere                  = hemisphere.lower()
        self.logger.info(f"hemisphere initialised: {self.hemisphere_dict['abbreviation']}")

    def slice_hemisphere(self, var_dict):
        sliced = {k: v.isel(nj=self.hemisphere_dict['nj_slice']) if k not in {'FI_OBS_CLI', 'preserved_vars'} else v for k, v in var_dict.items()}
        self.logger.info("hemisphere sliced on 'main' dataset")
        if "preserved_vars" in var_dict and isinstance(var_dict["preserved_vars"], dict):
            preserved = var_dict["preserved_vars"]
            sliced_preserved = { k: v.isel(nj=self.hemisphere_dict['nj_slice']) for k, v in preserved.items() }
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

    def get_static_grid(self, var):
        if "time" in var.dims:
            return var.isel(time=0)
        return var

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

    def load_bgrid(self):
        """
        Load the B-grid (t-grid and u-grid) and create cooordinates which are converted to degrees and standardized to [-180, 180].
        Sets attributes `self.G_t`, `self.G_u`.
        """
        G        = xr.open_dataset(self.CICE_dict['P_G'])
        KMT_org  = xr.open_dataset(self.P_KMT_org).kmt.data
        if self.use_gi:
            KMT_mod  = xr.open_dataset(self.P_KMT_mod).kmt.data
        else:
            KMT_mod = KMT_org
        TLAT     = self.radians_to_degrees(G['tlat'].data)
        TLON     = self.radians_to_degrees(G['tlon'].data)
        ULAT     = self.radians_to_degrees(G['ulat'].data)
        ULON     = self.radians_to_degrees(G['ulon'].data)
        TCOORDS  = self.build_grid_dict( TLAT , TLON )
        UCOORDS  = self.build_grid_dict( ULAT , ULON )
        T_ANGLE  = self.radians_to_degrees(G['angleT'].data)
        U_ANGLE  = self.radians_to_degrees(G['angle'].data)
        TAREA    = G['tarea'].data
        UAREA    = G['uarea'].data
        # Dimensions
        nj, ni           = TLAT.shape
        nj_b, ni_b       = nj + 1, ni + 1
        native_dim_names = ('nj','ni')
        native_dims      = (nj, ni)
        extend_dim_names = ('nj_b','ni_b')
        extend_dims      = (nj_b, ni_b)
        self.G_t = xr.Dataset(data_vars = { 'lat'     : (native_dim_names, TCOORDS['lat'].data   , {'units': 'degrees'}),
                                            'lat_b'   : (extend_dim_names, TCOORDS['lat_b'].data , {'units': 'degrees'}),
                                            'lon'     : (native_dim_names, TCOORDS['lon'].data   , {'units': 'degrees'}),
                                            'lon_b'   : (extend_dim_names, TCOORDS['lon_b'].data , {'units': 'degrees'}),
                                            'angle'   : (native_dim_names, T_ANGLE               , {'units': 'degrees'}),
                                            'area'    : (native_dim_names, TAREA                 , {'units': 'm^2'}),
                                            'kmt_org' : (native_dim_names, KMT_org               , {'units'      : 'binary',
                                                                                                    'description': '1=land, 0=ocean',
                                                                                                    'long_name'  : 'original landmask on t-grid'}),
                                            'kmt_mod' : (native_dim_names, KMT_mod               , {'units'      : 'binary',
                                                                                                    'description': '1=land, 0=ocean',
                                                                                                    'long_name'  : 'modified t-grid-landmask to simulate grounded icebergs'})},
                               coords   = { 'nj'   : np.arange(nj),
                                            'ni'   : np.arange(ni),
                                            'nj_b' : np.arange(nj_b),
                                            'ni_b' : np.arange(ni_b)})
        self.G_u = xr.Dataset(data_vars = { 'lat'     : (native_dim_names, UCOORDS['lat'].data   , {'units': 'degrees'}),
                                            'lat_b'   : (extend_dim_names, UCOORDS['lat_b'].data , {'units': 'degrees'}),
                                            'lon'     : (native_dim_names, UCOORDS['lon'].data   , {'units': 'degrees'}),
                                            'lon_b'   : (extend_dim_names, UCOORDS['lon_b'].data , {'units': 'degrees'}),
                                            'angle'   : (native_dim_names, U_ANGLE               , {'units': 'degrees'}),
                                            'area'    : (native_dim_names, UAREA                 , {'units': 'm^2'}),
                                            'kmt_org' : (native_dim_names, KMT_org               , {'units'      : 'binary',
                                                                                                    'description': '1=land, 0=ocean',
                                                                                                    'long_name'  : 'original landmask on t-grid'}),
                                            'kmt_mod' : (native_dim_names, KMT_mod               , {'units'      : 'binary',
                                                                                                    'description': '1=land, 0=ocean',
                                                                                                    'long_name'  : 'modified t-grid-landmask to simulate grounded icebergs'})},
                              coords    = { 'nj'   : np.arange(nj),
                                            'ni'   : np.arange(ni),
                                            'nj_b' : np.arange(nj_b),
                                            'ni_b' : np.arange(ni_b)})
        self.bgrid_loaded = True

    def align_modified_landmask(self):
        """
        Compute the grounded iceberg mask and coordinates from the difference
        between original and modified KMT arrays. Store in self.G_t.
        """
        if not self.bgrid_loaded:
            self.load_bgrid()
        kmt_mod = self.G_t['kmt_mod'].data
        kmt_org = self.G_t['kmt_org'].data
        lat     = self.G_t['lat'].data
        lon     = self.G_t['lon'].data
        area    = self.G_t['area'].data
        # Difference: grounded icebergs are cells that changed from ocean (1) to land (0)
        diff_mask           = (kmt_org == 1) & (kmt_mod == 0)
        self.G_t['GI_mask'] = (('nj', 'ni'), diff_mask)
        # Get coordinates of affected cells (shifted west by one ni index to match B-grid layout)
        nj_idx, ni_idx = np.where(diff_mask)
        ni_idx_shifted = ni_idx - 1
        valid          = ni_idx_shifted >= 0
        nj_idx         = nj_idx[valid]
        ni_idx         = ni_idx_shifted[valid]
        GI_lat         = lat[nj_idx, ni_idx]
        GI_lon         = lon[nj_idx, ni_idx]
        # Save 1D arrays with iceberg IDs
        iceberg_id         = np.arange(len(GI_lat))
        self.G_t['GI_lat'] = (('iceberg_id',), GI_lat)
        self.G_t['GI_lon'] = (('iceberg_id',), GI_lon)
        self.G_t['GI_lat'].attrs.update({'long_name': 'Grounded iceberg latitude'})
        self.G_t['GI_lon'].attrs.update({'long_name': 'Grounded iceberg longitude'})
        #print(f"{len(iceberg_id)} circumpolar grounded icebergs associated with {self.sim_name}")
        self.modified_landmask_aligned = True
        return self.G_t

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
        af2020_df    = pd.read_csv(self.AF_FI_dict['P_AF2020_cli_csv'])
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
        kmt             = kmt.isel(nj=self.hemisphere_dict['nj_slice']).values
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
