import os, sys, time, json, logging, re
import xarray             as xr
import xesmf              as xe
import pandas             as pd
import numpy              as np
import matplotlib.pyplot  as plt
from pathlib              import Path
from datetime             import datetime, timedelta
from dask.distributed     import Client, LocalCluster, get_client
from distributed.client   import _get_global_client
sys.path.insert(0, '/home/581/da1339/AFIM/src/AFIM/src')
from sea_ice_classification import SeaIceClassification
from sea_ice_metrics        import SeaIceMetrics
from sea_ice_plotter        import SeaIcePlotter
from sea_ice_icebergs       import SeaIceIcebergs
from sea_ice_observations   import SeaIceObservations
from sea_ice_ACCESS         import SeaIceACCESS

class SeaIceToolbox(SeaIceClassification, SeaIceMetrics, SeaIcePlotter, SeaIceIcebergs, SeaIceObservations, SeaIceACCESS):
    """
    Unified toolbox for processing and analysing Antarctic sea ice from CICE simulations
    as part of the Antarctic Fast Ice Model (AFIM) workflow.

    This superclass consolidates functionality from:
    + SeaIceClassification : fast/pack ice masking and simulation I/O
    + SeaIceMetrics        : computing sea ice measurements and statistics 
    + SeaIcePlotter        : PyGMT and timeseries visualization
    + SeaIceIcebergs       : grounded iceberg thinning/masking
    + SeaIceObservations   : Fraser et al. (2020) and NSIDC comparisons
    + SeaIceACCESS         : ACCESS-OM sea ice generated result; data extraction

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
    ice_vector_type : str or list of str, optional
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
        print(f"Speed Type(s)       : {self.ivec_type}")
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
                 ice_vector_type              = None,# a valid ice_vector_type or list thereof
                 ice_type                    = None,# a valid ice_type or list thereof
	             mean_period                 = None,# rolling average, N-days
                 boolean_window              = None,# the window with which to apply the minimum number of days
                 boolean_min_days            = None,# minimum number of days binary-days
	             extra_cice_vars             = None,# these will be included in the fast ice mask
                                                    # that is, in addtion to those listed in the
                                                    # in config file 'cice_vars_reqd'; default is
                                                    # also a list in config file 'FI_cice_vars_ext'
                                                    # can be set to True and will use those listed in the JSON config file
	             hemisphere                  = None,# used in many ares of the toolbox to define the hemisphere
                                                    # that the user is interested; unfortunately, the toolbox
                                                    # does not allow for a user to be interested in both at the 
                                                    # same time; either 'south' or 'north'; defualt is 'south'
	             P_log                       = None,# the log file to send print statements to
                 log_level                   = None,# the logging level (see python logging doc for more info)
                 overwrite_zarr              = None,# whether or not to overwrite a zarr; default is false
	             overwrite_AF2020_zarr       = None,# whether or not to overwrite AF2020db zarr; default is false
                 overwrite_saved_figs        = None,# whether or not to overwite saved figures; default is false
                 save_new_figs               = None,# whether or not to write new figures to disk; default is true
                 show_figs                   = None,# whether or not to show/print figures to screen; default is false
                 delete_original_cice_iceh_nc= None,# whether or not to delete the original CICE ice history
                 client                      = None,# dask distributed client, can be externally passed here
                 **kwargs):
        self.sim_name = sim_name if sim_name is not None else 'test'
        if P_json is None:
            P_json = '/home/581/da1339/AFIM/src/AFIM/src/JSONs/sea_ice_config.json'
        with open(P_json, 'r') as f:
            self.config = json.load(f)
        self.D_dict             = self.config.get('D_dict'            , {})
        D_log                   = self.D_dict['logs']
        P_log                   = P_log     if P_log     is not None else Path(D_log, f'SeaIceProcessor_FI_{self.sim_name}.log')
        log_level               = log_level if log_level is not None else logging.INFO
        if not os.path.exists(P_log):
            os.system(f"touch {P_log}")
        self.setup_logging(logfile=P_log, log_level=log_level)
        self.initialise_dask_client(client=client)
        self.CICE_dict          = self.config.get("CICE_dict"         , {})
        self.GI_dict            = self.config.get('GI_dict'           , {})
        self.NSIDC_dict         = self.config.get('NSIDC_dict'        , {}) 
        self.AF_FI_dict         = self.config.get("AF_FI_dict"        , {})
        self.Sea_Ice_Obs_dict   = self.config.get("Sea_Ice_Obs_dict"  , {})
        self.AOM2_dict          = self.config.get("AOM2_dict"         , {})
        self.MOM_dict           = self.config.get("MOM_dict"          , {})
        self.ERA5_dict          = self.config.get("ERA5_dict"         , {})
        self.ORAS_dict          = self.config.get("ORAS_dict"         , {})
        self.pygmt_dict         = self.config.get("pygmt_dict"        , {})
        self.plot_var_dict      = self.config.get("plot_var_dict"     , {})
        self.hemispheres_dict   = self.config.get("hemispheres_dict"  , {})
        self.Ant_8sectors       = self.config.get('Ant_8sectors'      , {})
        self.Ant_2sectors       = self.config.get('Ant_2sectors'      , {})
        self.plot_ice_area_dict = self.config.get('plot_ice_area_dict', {})
        self.dt0_str              = dt0_str                     if dt0_str                     is not None else self.config.get('dt0_str', '1993-01-01')
        self.dtN_str              = dtN_str                     if dtN_str                     is not None else self.config.get('dtN_str', '1999-12-31')
        self.mean_period          = mean_period                 if mean_period                 is not None else self.config.get('mean_period', 15)
        self.bool_window          = boolean_window              if boolean_window              is not None else self.config.get('bool_window', 7)
        self.bool_min_days        = boolean_min_days            if boolean_min_days            is not None else self.config.get('bool_min_days', 6)
        self.ispd_thresh          = ice_speed_threshold         if ice_speed_threshold         is not None else self.config.get('ice_speed_thresh_hi', 1e-3)
        self.ivec_type            = ice_vector_type              if ice_vector_type              is not None else self.config.get('ice_vector_type', 'ispd_BT')
        self.ice_type             = ice_type                    if ice_type                    is not None else self.config.get('ice_type', 'FI_BT')
        self.icon_thresh          = ice_concentration_threshold if ice_concentration_threshold is not None else self.config.get('ice_conc_thresh', 0.15)
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
        self.D_tmp               = Path(self.config['D_dict']['tmp'])
        self.D_metrics           = Path(self.D_zarr, f"ispd_thresh_{self.ispd_thresh_str}", "metrics")
        self.sim_config          = self.parse_simulation_metadata()   
        self.valid_ivec_types    = self.config.get("valid_ivec_types", [])
        self.valid_ice_types     = self.config.get("valid_ice_types", [])
        self.cice_vars_reqd      = self.CICE_dict["cice_vars_reqd"]
        self.spatial_dims        = self.CICE_dict["spatial_dims"]
        if extra_cice_vars is not None:
            if extra_cice_vars:
                self.cice_var_list = self.cice_vars_reqd + self.CICE_dict["cice_vars_ext"]
            else:
                self.cice_var_list = self.cice_vars_reqd + extra_cice_vars
        else:
            self.cice_var_list = self.cice_vars_reqd
        self.FIC_scale           = self.config.get('FIC_scale', 1e9)
        self.SIC_scale           = self.config.get('SIC_scale', 1e12)
        self.P_KMT_org           = Path(self.GI_dict["D_GI_thin"],self.GI_dict['KMT_org_fmt'])
        self.GI_thin             = self.sim_config.get('GI_thin_fact')
        self.GI_version          = self.sim_config.get('GI_version')
        self.GI_iteration        = self.sim_config.get("GI_iter")
        if self.GI_thin>0:
            self.use_gi    = True
            GI_thin_str    = f"{self.GI_thin:0.2f}".replace('.', 'p')
            GI_vers_str    = f"{self.GI_version:0.2f}".replace('.', 'p')
            self.P_KMT_mod = os.path.join(self.GI_dict['D_GI_thin'],
                                          self.GI_dict['KMT_mod_fmt'].format(GI_thin   = GI_thin_str,
                                                                             version   = GI_vers_str,
                                                                             iteration = self.GI_iteration))
        else:
            self.use_gi = False
        self.reG_weights_defined       = False
        self.modified_landmask_aligned = False
        self.bgrid_loaded              = False
        SeaIceClassification.__init__(self, sim_name, **kwargs)
        SeaIceMetrics.__init__(self, **kwargs)
        SeaIcePlotter.__init__(self, **kwargs)
        SeaIceIcebergs.__init__(self, **kwargs)
        SeaIceObservations.__init__(self, **kwargs)
        SeaIceACCESS.__init__(self, **kwargs)

    def setup_logging(self, logfile=None, log_level=logging.INFO):
        self.logger = logging.getLogger(self.sim_name)
        self.logger.setLevel(log_level)
        if not self.logger.handlers:
            formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
            ch = logging.StreamHandler()
            ch.setFormatter(formatter)
            ch.setLevel(log_level)
            self.logger.addHandler(ch)
            if logfile:
                if os.path.exists(logfile):
                    os.remove(logfile)
                fh = logging.FileHandler(logfile)
                fh.setFormatter(formatter)
                fh.setLevel(log_level)
                self.logger.addHandler(fh)
            self.logger.info(f"log file initialised: {logfile}")
    
    def initialise_dask_client(self, client=None, threads_per_worker=1):
        # Case 1: User has passed a client explicitly
        if client is not None:
            self.client = client
            self.logger.info("Using externally provided Dask client.")
        # Case 2: No client passed — try to detect existing global client
        elif _get_global_client() is not None:
            existing_client = _get_global_client()
            self.logger.warning("Dask client already exists but was not passed to SeaIceToolbox.")
            self.logger.warning("Please explicitly pass the existing Dask client using `client=...` to avoid confusion.")
            self.client = existing_client
        # Case 3: No client at all — create a new one
        else:
            self.client = Client(threads_per_worker=threads_per_worker, memory_limit="16GB")
            self.logger.info("Initialized new Dask client.")
        self.logger.info(f"Dask distributed client can be accessed at url {self.client.dashboard_link}")

    def parse_simulation_metadata(self, force_recompile=False):
        """
        Parse or load simulation metadata for CICE AFIM runs.

        If a cached JSON config file exists in the simulation directory and
        `force_recompile=False`, this method will load and return the metadata from that file.
        Otherwise, it will parse the `ice_diag.d` file and store the extracted parameters
        in a structured JSON file for future use.

        Parameters
        ----------
        force_recompile : bool
            If True, re-parse the raw ice_diag.d file even if a cached JSON exists.

        Returns
        -------
        dict
            Dictionary of extracted parameters including inferred grounded iceberg info.
        """
        sim_name   = Path(self.D_sim).name
        P_diag     = Path(self.D_sim, "ice_diag.d")
        P_json     = Path(self.D_sim, f"ice_in_AFIM_subset_{sim_name}.json")
        # Load from cache if it exists
        if P_json.exists() and not force_recompile:
            with open(P_json, "r") as f:
                return json.load(f)
        # Keys and patterns
        PARAM_KEYS = ["dt", "ndtd", "ndte", "kdyn", "revised_evp", "e_yieldcurve", "e_plasticpot",
                      "Ktens", "kstrength", "Pstar", "Cstar", "Cf", "visc_method", "kmt_file"]
        PATTERNS = {key: re.compile(rf"{key}\s*=\s*(.+?)\s*:") if key != "kmt_file"
                    else re.compile(rf"{key}\s*=\s*(.+)$")
                    for key in PARAM_KEYS}
        result = {key: "" for key in PARAM_KEYS}
        # Parse ice_diag.d
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
        # Additional derived metadata
        kmt_name = result["kmt_file"]
        if "kmt_mod_thinned-" in kmt_name:
            try:
                thin_str  = kmt_name.split("thinned-")[1].split("_")[0]  # e.g., "0p85"
                result["GI_thin_fact"] = float(thin_str.replace("p", "."))
            except Exception:
                result["GI_thin_fact"] = "?"
            try:
                version_match = re.search(r"_v(\d+p\d+)", kmt_name)
                if version_match:
                    version_str = version_match.group(1)  # '1p50'
                    result["GI_version"] = float(version_str.replace("p", "."))
                else:
                    self.logger.info(f"[DEBUG] No version match in: {kmt_name}")
                    result["GI_version"] = "?"
            except Exception as e:
                self.logger.info(f"[ERROR] GI_version failed for {kmt_name}: {e}")
                result["GI_version"] = "?"
            except Exception:
                result["GI_version"] = "?"
            try:
                if "iter" in kmt_name:
                    iter_str = kmt_name.split("iter")[1].split(".")[0]   # e.g., "0"
                    result["GI_iter"] = int(iter_str)
                else:
                    result["GI_iter"] = 0
            except Exception:
                result["GI_iter"] = None
        else:
            result["GI_thin_fact"] = 0.0
            result["GI_version"]   = 0.0
            result["GI_iter"]      = 0
        # Save to JSON
        self.logger.info(f"[CHECK] Parsed GI_version = {result.get('GI_version')} from {kmt_name}")
        with open(P_json, "w") as f:
            json.dump(result, f, indent=2)
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

    def reG_bgrid_to_tgrid_xesmf(self, da):
        """
        Regrid a single B-grid DataArray to the T-grid using pre-defined xESMF regridder.

        The input `da` must have coordinates 'ULAT' and 'ULON', which are renamed to 'lat' and 'lon'
        temporarily for compatibility with the xESMF regridder (which was defined using G_u/G_t).

        INPUTS:
        da : xr.DataArray; b-grid variable with coordinates 'ULAT' and 'ULON'.

        OUTPUTS:
        xr.DataArray; re-gridded DataArray on the T-grid with dimensions (time, nj, ni).
        """
        if "ULAT" not in da.coords or "ULON" not in da.coords:
            self.logger.error("❌ Cannot regrid: 'ULAT'/'ULON' not found in coordinates.")
            return None
        da_tmp = da.rename({"ULAT": "lat", "ULON": "lon"})
        try:
            da_reG = self.reG(da_tmp)
        except Exception as e:
            self.logger.error(f"❌ Regridding failed: {e}")
            return None
        return da_reG

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
        TLON_b,TLAT_b = self.build_grid_corners( TLAT , TLON )
        ULON_b,ULAT_b = self.build_grid_corners( ULAT , ULON )
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
        self.G_t = xr.Dataset(data_vars = { 'lat'     : (native_dim_names, TLAT   , {'units': 'degrees'}),
                                            'lat_b'   : (extend_dim_names, TLAT_b , {'units': 'degrees'}),
                                            'lon'     : (native_dim_names, TLON   , {'units': 'degrees'}),
                                            'lon_b'   : (extend_dim_names, TLON_b , {'units': 'degrees'}),
                                            'angle'   : (native_dim_names, T_ANGLE, {'units': 'degrees'}),
                                            'area'    : (native_dim_names, TAREA  , {'units': 'm^2'}),
                                            'kmt_org' : (native_dim_names, KMT_org, {'units'      : 'binary',
                                                                                     'description': '1=land, 0=ocean',
                                                                                     'long_name'  : 'original landmask on t-grid'}),
                                            'kmt_mod' : (native_dim_names, KMT_mod, {'units'      : 'binary',
                                                                                     'description': '1=land, 0=ocean',
                                                                                     'long_name'  : 'modified t-grid-landmask to simulate grounded icebergs'})},
                               coords   = { 'nj'   : np.arange(nj),
                                            'ni'   : np.arange(ni),
                                            'nj_b' : np.arange(nj_b),
                                            'ni_b' : np.arange(ni_b)})
        self.G_u = xr.Dataset(data_vars = { 'lat'     : (native_dim_names, ULAT   , {'units': 'degrees'}),
                                            'lat_b'   : (extend_dim_names, ULAT_b , {'units': 'degrees'}),
                                            'lon'     : (native_dim_names, ULON   , {'units': 'degrees'}),
                                            'lon_b'   : (extend_dim_names, ULON_b , {'units': 'degrees'}),
                                            'angle'   : (native_dim_names, U_ANGLE, {'units': 'degrees'}),
                                            'area'    : (native_dim_names, UAREA  , {'units': 'm^2'}),
                                            'kmt_org' : (native_dim_names, KMT_org, {'units'      : 'binary',
                                                                                     'description': '1=land, 0=ocean',
                                                                                     'long_name'  : 'original landmask on t-grid'}),
                                            'kmt_mod' : (native_dim_names, KMT_mod, {'units'      : 'binary',
                                                                                     'description': '1=land, 0=ocean',
                                                                                     'long_name'  : 'modified t-grid-landmask to simulate grounded icebergs'})},
                              coords    = { 'nj'   : np.arange(nj),
                                            'ni'   : np.arange(ni),
                                            'nj_b' : np.arange(nj_b),
                                            'ni_b' : np.arange(ni_b)})
        self.bgrid_loaded = True