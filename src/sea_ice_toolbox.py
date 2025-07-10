import os, sys, time, json, logging, re, psutil
import xarray             as xr
import xesmf              as xe
import pandas             as pd
import numpy              as np
import matplotlib.pyplot  as plt
from collections          import defaultdict
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
        P_log                   = P_log     if P_log     is not None else Path(D_log, f'SeaIceToolbox_{self.sim_name}.log')
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
        self.ivec_type            = ice_vector_type             if ice_vector_type             is not None else self.config.get('ice_vector_type', 'ispd_BT')
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
        if self.sim_config is not None:
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
        else:
            self.GI_thin             = None
            self.GI_version          = None
            self.GI_iteration        = None
            self.use_gi              = None
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
        if P_json.exists() and not force_recompile:
            with open(P_json, "r") as f:
                return json.load(f)
        PARAM_KEYS = ["dt", "ndtd", "ndte", "kdyn", "revised_evp", "e_yieldcurve", "e_plasticpot",
                      "Ktens", "kstrength", "Pstar", "Cstar", "Cf", "visc_method", "kmt_file"]
        PATTERNS = {key: re.compile(rf"{key}\s*=\s*(.+?)\s*:") if key != "kmt_file"
                    else re.compile(rf"{key}\s*=\s*(.+)$")
                    for key in PARAM_KEYS}
        result = {key: "" for key in PARAM_KEYS}
        try:
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
        except Exception:
            return None
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

    def create_monthly_strings(self, dt0_str=None, dtN_str=None):
        dt0_str = dt0_str if dt0_str is not None else self.dt0_str
        dtN_str = dtN_str if dtN_str is not None else self.dtN_str
        dts     = pd.date_range(dt0_str, dtN_str, freq="D")
        return sorted(set(dt.strftime("%Y-%m") for dt in dts))

    def create_empty_valid_DS_dictionary(self, valid_zarr_DS_list=None):
        valid_DS_list = valid_zarr_DS_list if valid_zarr_DS_list is not None else self.valid_ice_types
        return defaultdict(lambda: {k: [] for k in valid_DS_list})

    def count_zarr_files(self, path):
        total_files = sum(len(files) for _, _, files in os.walk(path))
        self.logger.info(f"{path} contains {total_files} files")

    def radians_to_degrees(self, da):
        return (da * 180) / np.pi

    def get_static_grid(self, var):
        if "time" in var.dims:
            return var.isel(time=0)
        return var

    def normalise_longitudes(self,lon):
        return ((lon + 360) % 360) - 180

    def build_grid_corners(self, lat_rads, lon_rads, grid_res=0.25, source_in_radians=True):
        if source_in_radians:
            lon = self.radians_to_degrees(lon_rads)
            lat = self.radians_to_degrees(lat_rads)
        else:
            lon = lon_rads
            lat = lat_rads
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

    def dict_to_ds(self, data_dict):
        return xr.Dataset({k: v for k, v in data_dict.items()})

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
        chunks          = chunks      or self.CICE_dict['FI_chunks']
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
            self.logger.info(f"searching for files between {dt0} and {dtN} in {Path(D_zarr,f'ispd_thresh_{ispd_thresh_str}')}")
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
                ds = xr.open_zarr(P_zarr, group=ice_type, consolidated=True)
                DS_list.append(ds)
            except (OSError, KeyError) as e:
                self.logger.warning(f"Skipping {P_zarr} ({ice_type}): {e}")
        if not DS_list:
            self.logger.warning(f"No {ice_type} datasets found in any Zarr group")
            return None
        DS_FI = xr.concat(DS_list, dim="time", coords='minimal').chunk(chunks)
        self.logger.info(f"Loaded {ice_type}: {len(DS_FI.time)} time steps from {len(DS_list)} files")
        self.logger.info(f"Memory after Zarr load: {psutil.virtual_memory().percent:.1f}% used")
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
                                     parallel   = True)
            CICE = CICE.chunk(chunks)
            if slice_hem:
                CICE = CICE.isel(nj=self.hemisphere_dict['nj_slice'])
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

    def coarsen_and_align_simulated_FI_to_observed_FI(self, sim_ds, obs_ds, doy_vals=None, method="mean", obs_time_coord="time"):
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
        obs_times = pd.to_datetime(obs_ds[obs_time_coord].values)
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
                ds_agg = ds_agg.expand_dims({obs_time_coord: [dt_start]})
                grouped.append(ds_agg)
        if not grouped:
            raise ValueError("❌ No observation periods matched simulation data")
        ds_aligned = xr.concat(grouped, dim=obs_time_coord)
        self.logger.info(f"✅ Aligned model output to {len(ds_aligned[obs_time_coord])} obs windows")
        return ds_aligned

    def align_time_coordinate_of_three_arrays(self, ds1, ds2, ds3, time_coord="time"):
        for da in [ds1, ds2, ds3]:
            da[time_coord] = pd.to_datetime(da[time_coord].values).normalize()
        t_common = np.intersect1d(np.intersect1d(ds1[time_coord].values, ds2[time_coord].values), ds3[time_coord].values)
        return ds1.sel(time=t_common), ds2.sel(time=t_common), ds3.sel(time=t_common)