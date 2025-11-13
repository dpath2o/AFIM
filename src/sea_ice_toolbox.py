import os, sys, time, json, logging, re, psutil, dask, tempfile
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
import warnings
warnings.filterwarnings( "ignore", message="Sending large graph of size", category=UserWarning, module="distributed.client")
sys.path.insert(0, '/home/581/da1339/AFIM/src/AFIM/src')
from sea_ice_classification import SeaIceClassification
from sea_ice_metrics        import SeaIceMetrics
from sea_ice_plotter        import SeaIcePlotter
from sea_ice_icebergs       import SeaIceIcebergs
from sea_ice_observations   import SeaIceObservations
from sea_ice_ACCESS         import SeaIceACCESS
from sea_ice_gridwork       import SeaIceGridWork
from sea_ice_regridder      import SeaIceRegridder

__all__ = ["SeaIceToolbox", "SeaIceToolboxManager"]

##################################################################################################

class SeaIceToolboxManager:
    """
    Lightweight factory for a shared Dask client and preconfigured `SeaIceToolbox`.

    Parameters
    ----------
    P_log : str | Path
        Path to a log file to attach to the toolbox logger.

    Notes
    -----
    - `get_shared_dask_client()` creates a class-level client (cached) so multiple
      toolboxes reuse the same cluster.
    - `get_toolbox(sim_name, **kwargs)` returns a `SeaIceToolbox` bound to the
      shared client and your log file.

    See Also
    --------
    SeaIceToolbox
        The main analysis class composed from SeaIce* modules.
    """
    _shared_client = None

    def __init__(self, P_log,
                 n_workers : int  = 4,
                 n_threads : int  = 1,
                 mem_lim   : str  = "16GB",
                 process   : bool = True,
                 D_dask    : str  = None):
        D_dask     = D_dask if D_dask is not None else os.environ.get("DASK_TEMPORARY_DIRECTORY", tempfile.gettempdir())
        self.P_log = P_log
        if SeaIceToolboxManager._shared_client is None:
            LocCls = LocalCluster(n_workers         = n_workers, 
                                  threads_per_worker = n_threads, # use all CPUs
                                  processes          = process,   # threads, not processes
                                  memory_limit       = mem_lim,   # no nanny hard limit
                                  local_directory    = D_dask,
                                  dashboard_address  = None)
            SeaIceToolboxManager._shared_client = Client(LocCls)

    def shutdown(self):
        if SeaIceToolboxManager._shared_client is not None:
            SeaIceToolboxManager._shared_client.close()
            print("Dask client shut down.")
            SeaIceToolboxManager._shared_client = None
        logger = logging.getLogger("SeaIceToolbox")
        for h in logger.handlers[:]:
            if isinstance(h, logging.FileHandler):
                h.flush(); h.close(); logger.removeHandler(h)
                print(f"Closed log file handler: {getattr(h, 'baseFilename', h)}")

    def get_toolbox(self, sim_name, **kwargs):
        '''
        This requires a simulation name and can pass any of the accepted agruments that SeaIceToolbox takes
        '''
        return SeaIceToolbox(sim_name = sim_name,
                             client   = SeaIceToolboxManager._shared_client,  # explicit
                             P_log    = self.P_log, **kwargs)

####################################################################################################################

class SeaIceToolbox(SeaIceClassification, SeaIceMetrics, SeaIcePlotter, 
                    SeaIceIcebergs, SeaIceObservations, SeaIceACCESS, SeaIceGridWork, SeaIceRegridder):
    """
    Unified toolbox for processing and analysing Antarctic sea ice from CICE simulations
    as part of the Antarctic Fast Ice Model (AFIM) workflow.

    Composition
    -----------
    - SeaIceClassification : fast/pack ice masking and simulation I/O
    - SeaIceMetrics        : metrics and statistics
    - SeaIcePlotter        : PyGMT maps and time series
    - SeaIceIcebergs       : grounded iceberg thinning/masking
    - SeaIceObservations   : Fraser et al. (2020) & NSIDC integration
    - SeaIceACCESS         : ACCESS-OM sea-ice helpers
    - SeaIceGridWork       : 
    - SeaIceRegridder      : B→T regridding helpers

    Key external data
    -----------------
    - Fraser et al. (2020) fast-ice climatology (ESSD)
    - NSIDC daily sea-ice concentration
    - AFIM JSON config (paths, parameters)

    Parameters
    ----------
    P_json : str | Path, optional
        Path to AFIM JSON configuration file (defaults to project root if None).
    sim_name : str
        Simulation name (must match an entry under `AFIM_out` in the config).
    dt0_str, dtN_str : str, optional
        Analysis date bounds in `YYYY-MM-DD`.
    ice_concentration_threshold : float, optional
        Concentration cutoff for fast ice (default 0.15).
    ice_speed_threshold : float, optional
        Speed threshold (m/s) below which ice is considered fast (default 5e-4).
    ice_vector_type : str | list[str], optional
        Which speed field(s) to use: `'B'`, `'Ta'`, `'Tx'`, or `'BT'`.
    ice_type : str | list[str], optional
        Classification (“FI”, “PI”, etc.), combined later into `self.ice_class`.
    mean_period : int, optional
        N-day rolling window for speed smoothing (default 15).
    bin_win_days, bin_min_days : int, optional
        Binary-days window size and minimum days (defaults per config).
    hemisphere : {'south','north'}, optional
        Controls hemisphere slicing throughout the pipeline.
    P_log : Path, optional
        Log file to attach to the toolbox logger.
    overwrite_zarr / overwrite_AF2020_zarr / overwrite_saved_figs : bool, optional
        Overwrite behaviours for outputs and figures.
    save_new_figs, show_figs : bool, optional
        Figure output/interactive toggles.
    client : dask.distributed.Client, optional
        Pre-existing Dask client. If omitted, use `SeaIceToolboxManager`.

    Notes
    -----
    - The constructor normalises configuration, paths, thresholds, and date
      bounds, then prints a short summary to the logger.
    """
    def summary(self):
        """Print a summary of key simulation setup and configuration."""
        self.logger.info("--- SeaIceToolbox Summary ---")
        self.logger.info(f"Simulation Name     : {self.sim_name}")
        self.logger.info(f"Analysis Start Date : {self.dt0_str}")
        self.logger.info(f"Analysis End Date   : {self.dtN_str}")
        self.logger.info(f"Speed Threshold     : {self.ispd_thresh:.1e} m/s")
        self.logger.info(f"Speed Type(s)       : {self.BT_composite_grids}")
        self.logger.info(f"Ice Speed Name      : {self.ivec_type}")
        self.logger.info(f"Ice Type(s)         : {self.ice_type}")
        self.logger.info(f"Mean Period         : {self.mean_period} days")
        self.logger.info(f"Binary-days Window  : {self.bin_win_days} days")
        self.logger.info(f"Binary-days Min-Days: {self.bin_min_days}")
        self.logger.info(f"Using GI?           : {self.use_gi}")
        self.logger.info(f"Overwrite Zarr      : {self.overwrite_zarr_group}")
        self.logger.info(f"Save Figures        : {self.save_fig}")
        self.logger.info(f"Show Figures        : {self.show_fig}")
        self.logger.info(f"Hemisphere          : {self.hemisphere}")
        self.logger.info("------------------------------")

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
                 list_of_composite_grids     = None,# select list of ["B", "Ta", "Tb", "Tx"];  must contain two;
                                                    # default ["Ta","Tx"]
                 iceh_frequency              = None,# 'hourly', 'daily', 'monthly', 'yearly'
                                                    # defines the history files that will be used by this toolbox
	             ice_concentration_threshold = None,# almost should never be changed from default value of
                                                    # 0.15 (grid cell concentration)
	             ice_speed_threshold         = None,# a significantly important value in the determination
                                                    # and classification (masking) of fast ice; defualt value
                                                    # 5e-4 m/s
                 ice_vector_type             = None,# a valid ice_vector_type or list thereof
                 ice_type                    = None,# a valid ice_type or list thereof
	             mean_period                 = None,# rolling average, N-days
                 bin_win_days                = None,# the window of days with which to apply binary-days method
                 bin_min_days                = None,# minimum number of days binary-days
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
                 dask_memory_limit           = None,# provide the memory limit to dask, default is 16GB
                 overwrite_zarr              = None,# whether or not to overwrite a zarr; default is false
	             overwrite_AF2020_zarr       = None,# whether or not to overwrite AF2020db zarr; default is false
                 overwrite_saved_figs        = None,# whether or not to overwite saved figures; default is false
                 save_new_figs               = None,# whether or not to write new figures to disk; default is true
                 show_figs                   = None,# whether or not to show/print figures to screen; default is false
                 delete_original_cice_iceh_nc= None,# whether or not to delete the original CICE ice history
                 client                      = None,# dask distributed client, can be externally passed here
                 force_recompile_ice_in      = False,# reinitialise ice_in JSON file; see help self.parse_simulation_metadata()
                 **kwargs):
        """
        Construct the unified AFIM sea-ice toolbox for a given simulation/config.

        The constructor normalises configuration from a JSON file (`P_json`) and
        simulation metadata, sets directory paths, thresholds, and plotting/Dask
        options, and wires up mixins (classification, metrics, plotting, icebergs,
        observations). It also initialises hemisphere slicing and output naming.

        Parameters
        ----------
        P_json : str or pathlib.Path, optional
            Path to the AFIM configuration JSON. If None, uses project defaults.
        sim_name : str, optional
            Simulation name (must match a key under `AFIM_out` in the config).
        dt0_str, dtN_str : str, optional
            Inclusive analysis window in ``YYYY-MM-DD``. Defaults are pulled from
            the config if not provided.
        list_of_composite_grids : list[str], optional
            Which speed components to average into a composite (`"B"`, `"Ta"`, `"Tb"`, `"Tx"`).
            Defaults to ``['Ta','Tx']`` via config.
        iceh_frequency : {"hourly","daily","monthly","yearly"}, optional
            Source history cadence for CICE iceh inputs. Default from config is ``"daily"``.
        ice_concentration_threshold : float, optional
            Sea-ice concentration threshold used for fast-ice masking (default 0.15).
        ice_speed_threshold : float, optional
            Speed threshold (m/s) below which ice is treated as fast-ice. Defaults to
            a high/“strict” threshold from config (e.g., ``5e-4``).
        ice_vector_type : str or list[str], optional
            Speed field(s) to use (e.g., ``"B"``, ``"Ta"``, ``"Tx"``, or composite like
            ``"BT"``). Defaults from config.
        ice_type : str or list[str], optional
            Ice classification type(s) to compute (e.g., ``"FI"`` fast-ice). Defaults from config.
        mean_period : int, optional
            N-day rolling window length for speed smoothing (default 15 via config).
        bin_win_days, bin_min_days : int, optional
            Binary-days window size and required number of days (defaults from config).
        extra_cice_vars : list[str] or bool, optional
            Additional CICE variables to include in outputs. If True, uses the
            config’s `cice_vars_ext`; if a list, uses that list; if None, only the
            required minimal set is used.
        hemisphere : {'south','north','sh','nh', ...}, optional
            Hemisphere of interest. Various synonyms are accepted; default from config.
        P_log : str or pathlib.Path, optional
            Log file to attach a file handler to. If omitted, logging goes to the
            preconfigured handlers only.
        log_level : int or str, optional
            Python logging level for this toolbox’s logger (e.g., ``logging.INFO``).
        dask_memory_limit : str or int, optional
            Memory limit to pass to a created Dask client (e.g., ``"16GB"``).
        overwrite_zarr, overwrite_AF2020_zarr : bool, optional
            Overwrite control flags for Zarr outputs.
        overwrite_saved_figs : bool, optional
            Whether to overwrite previously saved figures.
        save_new_figs, show_figs : bool, optional
            Figure saving and interactive display toggles.
        delete_original_cice_iceh_nc : bool, optional
            If True, delete original daily iceh NetCDFs after writing monthly Zarr.
        client : dask.distributed.Client, optional
            Reuse an existing Dask client; otherwise one may be obtained externally.
        force_recompile_ice_in : bool, default False
            Force regeneration/reparse of the simulation’s `ice_in` JSON metadata.
        **kwargs
            Extra fields are attached verbatim to `self`, allowing downstream mixins
            to see additional tuning parameters.

        Notes
        -----
        - Directory attributes like `self.D_sim`, `self.D_iceh`, `self.D_zarr`,
        and metrics directories are derived from the config and `sim_name`.
        - `self.define_hemisphere()` is called during initialisation.
        - Mixins are initialised at the end to ensure shared attributes are in place.
        """
        # essentially high-level administrative work:
        self.sim_name = sim_name if sim_name is not None else 'test'
        if P_json is None:
            P_json = '/home/581/da1339/AFIM/src/AFIM/src/JSONs/sea_ice_config.json'
        with open(P_json, 'r') as f:
            self.config = json.load(f)
        self.D_dict = self.config.get('D_dict'            , {})
        D_log       = self.D_dict['logs']
        P_log       = P_log     if P_log     is not None else Path(D_log, f'SeaIceToolbox_{self.sim_name}.log')
        log_level   = log_level if log_level is not None else logging.INFO
        if not os.path.exists(P_log):
            os.system(f"touch {P_log}")
        if not hasattr(self, 'logger'):
            self.setup_logging(logfile=P_log, log_level=log_level)
        dask_memory_limit = dask_memory_limit if dask_memory_limit is not None else "7GB"
        if client is not None:
            self.client = client
        if not hasattr(self, "client"):
            raise ValueError("Dask client must be provided explicitly or via manager.")
        self.logger.info(f"Dask Client Connected\n"
                         f"  Dashboard      : {self.client.dashboard_link}\n"
                         f"  Threads        : {sum(w['nthreads'] for w in client.scheduler_info()['workers'].values())}\n"
                         f"  Threads/Worker : {[w['nthreads'] for w in client.scheduler_info()['workers'].values()]}\n"
                         f"  Total Memory   : {sum(w['memory_limit'] for w in client.scheduler_info()['workers'].values()) / 1e9:.2f} GB\n")
        if sim_name=="__SI-toolbox-mgr__":
            return
        # now for the technical and sim-specific configurations
        self.leap_year          = self.config.get("leap_year"         , 1996)
        self.CICE_dict          = self.config.get("CICE_dict"         , {})
        self.GI_dict            = self.config.get('GI_dict'           , {})
        self.NSIDC_dict         = self.config.get('NSIDC_dict'        , {})
        self.BAS_dict           = self.config.get('BAS_dict'          , {}) 
        self.AF_FI_dict         = self.config.get("AF_FI_dict"        , {})
        self.Sea_Ice_Obs_dict   = self.config.get("Sea_Ice_Obs_dict"  , {})
        self.AOM2_dict          = self.config.get("AOM2_dict"         , {})
        self.MOM_dict           = self.config.get("MOM_dict"          , {})
        self.ERA5_dict          = self.config.get("ERA5_dict"         , {})
        self.ORAS_dict          = self.config.get("ORAS_dict"         , {})
        self.plot_var_dict      = self.config.get("plot_var_dict"     , {})
        self.hemispheres_dict   = self.config.get("hemispheres_dict"  , {})
        self.Ant_8sectors       = self.config.get('Ant_8sectors'      , {})
        self.Ant_2sectors       = self.config.get('Ant_2sectors'      , {})
        self.pygmt_dict         = self.config.get("pygmt_dict"        , {})
        self.pygmt_FIA_dict     = self.config.get('pygmt_FIA_dict'    , {})
        self.pygmt_FI_panel     = self.config.get('pygmt_FI_panel'    , {})
        self.dt0_str              = dt0_str                     if dt0_str                     is not None else self.config.get('dt0_str', '1993-01-01')
        self.dtN_str              = dtN_str                     if dtN_str                     is not None else self.config.get('dtN_str', '1999-12-31')
        self.BT_composite_grids   = list_of_composite_grids     if list_of_composite_grids     is not None else self.config.get('BT_composite_grids', ['Ta','Tx'])
        self.iceh_freq            = iceh_frequency              if iceh_frequency              is not None else self.config.get('iceh_freq', 'daily')
        self.mean_period          = mean_period                 if mean_period                 is not None else self.config.get('mean_period', 15)
        self.bin_win_days         = bin_win_days                if bin_win_days                is not None else self.config.get('bin_win_days', 11)
        self.bin_min_days         = bin_min_days                if bin_min_days                is not None else self.config.get('bin_min_days', 9)
        self.ispd_thresh          = ice_speed_threshold         if ice_speed_threshold         is not None else self.config.get('ice_speed_thresh_hi', 5.0e-4)
        self.ivec_type            = ice_vector_type             if ice_vector_type             is not None else self.config.get('ice_vector_type', 'BT')
        self.ice_type             = ice_type                    if ice_type                    is not None else self.config.get('ice_type', 'FI')
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
        self.sim_config          = self.parse_simulation_metadata(force_recompile=force_recompile_ice_in)   
        self.valid_ivec_types    = self.config.get("valid_ivec_types", [])
        self.valid_ice_types     = self.config.get("valid_ice_types", [])
        self.cice_vars_reqd      = self.CICE_dict["cice_vars_reqd"]
        self.spatial_dims        = self.CICE_dict["spatial_dims"]
        self.define_ispd_thresh_dir()
        self.define_ice_class_name()
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
            if self.GI_thin is not None:
                self.use_gi    = True
                GI_thin_str    = f"{self.GI_thin:0.2f}".replace('.', 'p')
                GI_vers_str    = f"{self.GI_version:0.2f}".replace('.', 'p')
                self.P_KMT_mod = os.path.join(self.GI_dict['D_GI_thin'],
                                              self.GI_dict['KMT_mod_fmt'].format(GI_thin   = GI_thin_str,
                                                                                 version   = GI_vers_str,
                                                                                 iteration = self.GI_iteration))
                self.P_GI_thin = os.path.join(self.GI_dict['D_GI_thin'],
                                              self.GI_dict['GI_thin_fmt'].format(GI_thin   = GI_thin_str,
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
        self.summary()

    def setup_logging(self, logfile=None, log_level=logging.INFO):
        logger_name = "sea_ice_classification"
        self.logger = logging.getLogger(logger_name)
        self.logger.setLevel(log_level)
        self.logger.propagate = False
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        # === Remove any old file handlers pointing to other files ===
        if logfile:
            for h in list(self.logger.handlers):
                if isinstance(h, logging.FileHandler):
                    self.logger.removeHandler(h)
                    h.close()
        # === Add stream handler if none exists ===
        if not any(isinstance(h, logging.StreamHandler) for h in self.logger.handlers):
            ch = logging.StreamHandler()
            ch.setFormatter(formatter)
            ch.setLevel(log_level)
            self.logger.addHandler(ch)
        # === Add (new) file handler ===
        if logfile:
            fh = logging.FileHandler(logfile, mode='a')  # always attach new one
            fh.setFormatter(formatter)
            fh.setLevel(log_level)
            self.logger.addHandler(fh)
            self.logger.info(f"log file connected: {logfile}")

    def _method_name(self):
        import inspect
        return inspect.currentframe().f_back.f_code.co_name

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
        self.logger.info(f"reading {P_diag} to construct {P_json}")
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
                    result["GI_iter"] = None
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
        """
        Initialise the hemisphere configuration for analysis.

        Maps common inputs (e.g. 'north', 'NH', 'sh') to the internal
        `self.hemisphere_dict` used for slicing, and converts the configured
        `nj_slice` tuple to a Python `slice`.

        Sets
        ----
        self.hemisphere_dict : dict
            Slicing and naming metadata (includes `nj_slice` and abbreviation).
        self.hemisphere : str
            Canonical lowercase hemisphere name.

        Raises
        ------
        ValueError
            If the input does not match any accepted names.
        """
        key = hemisphere.lower()
        if key in ['north', 'northern', 'nh', 'n', 'no']:
            self.hemisphere_dict = self.hemispheres_dict['north']
        elif key in ['south', 'southern', 'sh', 's', 'so']:
            self.hemisphere_dict = self.hemispheres_dict['south']
        else:
            raise ValueError(f"Invalid hemisphere '{hemisphere}'. Valid options are: "
                             "['north', 'south', 'northern', 'southern', 'sh', 'nh', 'SH', 'NH']")
        start, stop = self.hemisphere_dict['nj_slice']
        self.hemisphere_dict['nj_slice'] = slice(start, stop)
        self.hemisphere = key
        self.logger.info(f"hemisphere initialised: {self.hemisphere_dict['abbreviation']}")

    @staticmethod
    def _expand_bounds_slice(slc: slice, nb_len: int) -> slice:
        """
        Convert a centre-row slice (length nb_len-1) to a corner-row slice (length nb_len)
        by expanding the stop bound by +1, clipped to nb_len.
        """
        if not isinstance(slc, slice):
            raise TypeError(f"hemisphere nj_slice must be a Python slice, got {type(slc)}")
        start = 0 if slc.start is None else slc.start
        stop_nj = nb_len - 1              # last valid nj index is nb_len-2
        stop = stop_nj if slc.stop is None else slc.stop
        step = slc.step
        stop_b = min(stop + 1, nb_len)    # include boundary row at stop
        return slice(start, stop_b, step)

    @staticmethod
    def _indexers_for(obj, y_dim: str, nj_slice: slice):
        """
        Build an `.isel()` indexer dict for an object with optional centre/corner dims.

        - Applies `nj_slice` to `y_dim` if present.
        - If a matching corner dim `'nj_b'` exists, expands the slice using
          `_expand_bounds_slice()` to include the boundary row.
        """
        idx = {}
        dims = getattr(obj, "dims", {})
        if y_dim in dims:
            idx[y_dim] = nj_slice
        if "nj_b" in dims:
            nb_len = obj.sizes["nj_b"]
            idx["nj_b"] = SeaIceToolbox._expand_bounds_slice(nj_slice, nb_len)
        return idx

    def slice_hemisphere(self, var):
        """
        Apply the configured hemisphere slice to a Dataset/DataArray (and corners).

        Uses `self.hemisphere_dict['nj_slice']` on the cell-centre row dimension
        (e.g., `nj`). If a matching corner dimension is present (e.g., `nj_b`),
        expands the stop bound by +1 to include the boundary row.

        Parameters
        ----------
        var : xr.Dataset | xr.DataArray | dict[str, xr.Dataset | xr.DataArray]
            Object(s) to slice. Dict values are sliced if their dims match.

        Returns
        -------
        xr.Dataset | xr.DataArray | dict
            Sliced object of the same type as the input.

        Notes
        -----
        - Cell-centre dimension name is taken from `self.CICE_dict["y_dim"]`.
        - Corner dim is assumed to be `'nj_b'`. If you parameterise corners in
          your config, adapt `_indexers_for()` accordingly.
        """
        y_dim    = self.CICE_dict["y_dim"]          # typically 'nj'
        nj_slice = self.hemisphere_dict['nj_slice'] # a Python slice
        def _apply(obj):
            idx = SeaIceToolbox._indexers_for(obj, y_dim, nj_slice)
            return obj.isel(idx) if idx else obj
        if isinstance(var, dict):
            out = {}
            for k, v in var.items():
                out[k] = _apply(v) if isinstance(v, (xr.Dataset, xr.DataArray)) else v
            which = f"{y_dim} and/or nj_b"
            self.logger.debug(f"Hemisphere slice applied on dict members where dims matched ('{which}').")
            return out
        elif isinstance(var, (xr.Dataset, xr.DataArray)):
            sliced = _apply(var)
            which = " & ".join([d for d in (y_dim, "nj_b") if d in getattr(var, "dims", {})]) or "none"
            self.logger.debug(f"Hemisphere slice applied on dims: {which}.")
            return sliced
        else:
            raise ValueError(f"Unsupported input type: {type(var)}. Must be dict, Dataset, or DataArray.")

    def define_iceh_dirs(self, D_sim=None, iceh_freq=None):
        D_sim              = D_sim     or self.D_sim
        iceh_freq          = iceh_freq or self.iceh_freq
        self.D_iceh_netcdf = Path(D_sim, "history", iceh_freq)
        self.D_iceh_zarr   = Path(D_sim, "zarr", f"iceh_{iceh_freq}.zarr")

    def define_ispd_thresh_dir(self, D_zarr=None, ispd_thresh=None):
        """
        Define the output directory path for a given ice speed threshold.

        Parameters
        ----------
        D_zarr : str or Path, optional
            Base path to the parent Zarr directory. Defaults to `self.D_zarr`.
        ispd_thresh : float, optional
            Ice speed threshold to be included in the subdirectory name. Defaults to `self.ispd_thresh`.

        Sets
        ----
        self.D_ispd_thresh : Path
            Full path to the threshold-specific subdirectory, e.g., "ispd_thresh_5.0e-4".
        """
        D_zarr             = D_zarr      or self.D_zarr
        ispd_thresh        = ispd_thresh or self.ispd_thresh
        ispd_thresh_str    = f"{ispd_thresh:.1e}".replace("e-0", "e-").replace("e+0", "e+")
        self.D_ispd_thresh = Path(D_zarr, f"ispd_thresh_{ispd_thresh_str}")

    def _check_ivec_type(self,ivec_type):
        if isinstance(ivec_type, str):
            ivec_type = [ivec_type]
        assert all(v in self.valid_ivec_types for v in ivec_type), f"Invalid ivec_type: {ivec_type}"

    def _check_ice_type(self,ice_type):
        if isinstance(ice_type, str):
            ice_type = [ice_type]
        assert all(v in self.valid_ice_types for v in ice_type), f"Invalid ice_type: {ice_type}"

    def _as_da_mask(self, x):
        """
        Normalize a fast-ice mask input to an xarray.DataArray.

        Accepts either:
        • a DataArray that already is the binary fast-ice mask, or
        • a Dataset that contains a variable named 'FI_mask'.

        Returns
        -------
        xarray.DataArray
            The 'FI_mask' mask with its original coordinates/dtype preserved.

        Raises
        ------
        ValueError
            If a Dataset is provided but it does not contain 'FI_mask'.
        TypeError
            If the input is neither a DataArray nor a Dataset with 'FI_mask'.

        Notes
        -----
        This function does not cast dtype; callers may wish to `.astype('i1')`
        (or similar) if they want a compact integer mask.
        """
        if isinstance(x, xr.Dataset):
            if "FI_mask" in x:
                x = x["FI_mask"]
            else:
                raise ValueError("I_mask is a Dataset but lacks variable 'FI_mask'.")
        if not isinstance(x, xr.DataArray):
            raise TypeError("I_mask must be an xarray.DataArray or Dataset containing 'FI_mask'.")
        return x

    def _as_da_area(self, x):
        """
        Normalize a grid-cell area input to an xarray.DataArray.

        Accepts either:
        • a DataArray that already is the cell-area field, or
        • a Dataset containing one of {'tarea','area','TAREA'} (first match wins).

        Returns
        -------
        xarray.DataArray
            The area field with its original coordinates/dtype preserved.

        Raises
        ------
        ValueError
            If a Dataset is provided but no recognized area variable is present.
        TypeError
            If the input is neither a DataArray nor a Dataset with an area var.

        Notes
        -----
        This function does not alter dimensions. If the area has a 'time' dim,
        strip it upstream or downstream (e.g., `A.isel(time=0)`).
        """
        if isinstance(x, xr.Dataset):
            for k in ("tarea", "area", "TAREA"):
                if k in x:
                    x = x[k]
                    break
            else:
                raise ValueError("A is a Dataset but no area variable found (looked for 'tarea','area','TAREA').")
        if not isinstance(x, xr.DataArray):
            raise TypeError("A must be an xarray.DataArray or Dataset containing an area variable.")
        return x

    def _to_float_scalar(self, x):
        """
        Convert a scalar-like value (NumPy/xarray) to a Python float.

        Intended for outputs of reductions (e.g., `.sum()`, `.max()`) or
        `dask.compute(...)` that yield 0-D arrays.

        Parameters
        ----------
        x : Any
            NumPy scalar, 0-D ndarray, or 0-D xarray object.

        Returns
        -------
        float
            The scalar as a Python float (including NaN if present).

        Notes
        -----
        If a higher-dimensional array is passed accidentally, `float(np.asarray(x))`
        will raise; callers should only pass scalar-like values.
        """
        # Accept numpy scalars, xarray 0-d arrays, etc.
        if hasattr(x, "values"):
            x = x.values
        try:
            return float(x.item())  # numpy scalar or 0-d array
        except Exception:
            return float(np.asarray(x))
        
    def define_ice_class_name(self, ice_type=None , ivec_type=None ):
        """
        Define the classification name string for ice type and vector component type.

        Parameters
        ----------
        ice_type : str or list of str, optional
            Type(s) of ice classification (e.g., 'fast', 'drift'). Defaults to `self.ice_type`.
        ivec_type : str or list of str, optional
            Type(s) of ice velocity vector used (e.g., 'B', 'Ta', 'Tx'). Defaults to `self.ivec_type`.

        Raises
        ------
        AssertionError
            If any provided `ice_type` or `ivec_type` is not in the list of valid types.

        Sets
        ----
        self.ice_class : str
            Combined classification string in the format "{ice_type}_{ivec_type}".
        """
        ice_type  = ice_type  or self.ice_type
        ivec_type = ivec_type or self.ivec_type
        self._check_ivec_type(ivec_type) 
        self._check_ice_type(ice_type)  
        self.ice_class = f"{ice_type}_{ivec_type}"
        self.logger.info(f" self.ice_class defined as {self.ice_class}")

    def define_ice_speed_name(self, ivec_type=None):
        """
        Set the canonical name for the selected ice-speed vector type.

        Parameters
        ----------
        ivec_type : str, optional
            One of the valid vector types (e.g., ``"B"``, ``"Ta"``, ``"Tx"``, or
            composites depending on your config). Defaults to ``self.ivec_type``.

        Sets
        ----
        self.ispd_name : str
            Name used throughout outputs/paths, formatted as ``f"ispd_{ivec_type}"``.

        Raises
        ------
        ValueError
            If `ivec_type` is invalid (validated by `_check_ivec_type`).
        """
        ivec_type = ivec_type or self.ivec_type
        self._check_ivec_type(ivec_type) 
        self.ispd_name = f"ispd_{ivec_type}"

    def define_datetime_vars(self, dt0_str=None, dtN_str=None):
        """
        Define date range attributes from start and end date strings.

        Parameters
        ----------
        dt0_str : str, optional
            Start date string (e.g., '1994-01-01'). Defaults to `self.dt0_str`.
        dtN_str : str, optional
            End date string (e.g., '1999-12-31'). Defaults to `self.dtN_str`.

        Sets
        ----
        self.dt0 : pd.Timestamp
            Parsed start date.
        self.dtN : pd.Timestamp
            Parsed end date.
        self.yrs_mos : np.ndarray
            Array of 'YYYY-MM' strings for each month in the range.
        self.ymd_strs : np.ndarray
            Array of 'YYYY-MM-DD' strings for each day in the range.
        """
        from pandas.tseries.offsets import MonthEnd
        dt0_str       = dt0_str or self.dt0_str
        dtN_str       = dtN_str or self.dtN_str
        self.dt0      = pd.to_datetime(dt0_str)
        self.dtN      = pd.to_datetime(dtN_str)
        self.dt_range = pd.date_range(self.dt0, self.dtN, freq="D")
        self.ymd_strs = self.dt_range.strftime("%Y-%m-%d")
        self.mos0     = pd.date_range(self.dt0, self.dtN, freq="MS")
        self.mosN     = pd.date_range(self.dt0, self.dtN, freq="ME")
        self.yrs_mos0 = self.mos0.strftime("%Y-%m")
        self.mo0_strs = self.mos0.strftime("%Y-%m-%d").tolist()
        self.moN_strs = self.mosN.strftime("%Y-%m-%d").tolist()#(self.mos + MonthEnd(1)).strftime("%Y-%m-%d").tolist()
        self.yrs0     = pd.date_range(self.dt0, self.dtN, freq="YS")
        self.yrsN     = pd.date_range(self.dt0, self.dtN, freq="YE")
        self.yr0_strs = self.yrs0.strftime("%Y-%m-%d")
        self.yrN_strs = self.yrsN.strftime("%Y-%m-%d")

    def define_month_first_last_dates(self, year_month_str):
        """
        Given a 'YYYY-MM' string, return the first and last dates of that month.

        Parameters
        ----------
        year_month_str : str
            Year and month string in 'YYYY-MM' format.

        Returns
        -------
        tuple of str
            First and last day of the month in 'YYYY-MM-DD' format.
        """
        m0_str = f"{year_month_str}-01"
        mN_str = (pd.to_datetime(m0_str) + pd.offsets.MonthEnd()).strftime("%Y-%m-%d")
        return m0_str, mN_str

    def interpret_ice_speed_threshold(self, ispd_thresh=None, lat_thresh=-60):
        """
        Translate the ice speed threshold into intuitive grid-scale metrics.

        Computes:
        - meters per day at the given speed,
        - the **median** grid-cell edge length south of `lat_thresh` on the model
            grid (`self.CICE_dict['P_G']`),
        - displacement as a fraction of a grid cell per day,
        - days required to traverse a grid cell at the threshold speed.

        Parameters
        ----------
        ispd_thresh : float, optional
            Threshold in m/s. Defaults to ``self.ispd_thresh``.
        lat_thresh : float, default -60
            Latitude (degrees) used to select the polar region for the median cell size.

        Returns
        -------
        dict
            Summary metrics with keys:
            ``{'ice_speed_thresh_m_per_s', 'displacement_m_per_day',
            'median_grid_cell_length_m', 'percent_displacement_per_day',
            'days_per_grid_cell'}``.

        Notes
        -----
        - Uses CICE grid variables from ``self.CICE_dict['P_G']`` (radians) and
        converts to degrees for masking.
        - Logs a human-readable summary via `self.logger`.
        """
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
        if Path(P_iceh_zarr, ".zgroup").exists():
            self.logger.info(f"Deleting original NetCDF files for {m_str}")
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
            Zarr write. Defaults to ``self.delete_original_cice_iceh_nc``.

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
                                             engine   = "scipy",
                                             parallel = True,
                                             combine  = "by_coords",
                                             cache    = True,
                                             chunks   = {})
                CICE_all = CICE_all.chunk(self.CICE_dict['FI_chunks'])
                self.logger.info(f"Subtracting one day from original dataset as CICE reports one day ahead for daily-averages")
                CICE_all["time"] = CICE_all["time"] - np.timedelta64(1, "D")
                self.logger.info(f"Writing {P_iceh_zarr} and group ('YYYY-MM'): {m_str}")
                CICE_all.to_zarr(P_iceh_zarr, group=m_str, mode="w", consolidated=True)
                self.get_dir_size(P_iceh_zarr_group)
                self.count_zarr_files(P_iceh_zarr_group)
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

    def load_cice_zarr(self,
                    sim_name  = None,
                    dt0_str   = None,
                    dtN_str   = None,
                    D_alt     = None,
                    variables = None,
                    slice_hem = False):
        """
        Open monthly Zarr groups for the requested period and concatenate in time.

        Parameters
        ----------
        sim_name : str, optional
            Simulation name. Defaults to `self.sim_name`.
        dt0_str, dtN_str : str, optional
            Clamp the load to the intersection of requested and available data.
        D_alt : str | Path, optional
            Alternative simulation directory (overrides `self.D_sim`).
        variables : list[str], optional
            Subset to these variables; groups missing any are skipped with a warning.
        slice_hem : bool, default False
            If True, apply `slice_hemisphere()` after concatenation.

        Returns
        -------
        xr.Dataset
            Concatenated dataset with time cropped to `[user_dt0, user_dtN]`.

        Notes
        -----
        - Uses conservative Dask settings for slicing and fusion to keep graphs
          small (`split_large_chunks=True`, `chunk-size=256MiB`).
        - Requires Zarr groups of the form `YYYY-MM`.
        """
        import dask
        sim_name = sim_name or self.sim_name
        dt0_str  = dt0_str  or self.dt0_str
        dtN_str  = dtN_str  or self.dtN_str
        D_alt    = D_alt    or self.D_sim
        self.define_iceh_dirs(D_alt)
        # Step 1: Identify available YYYY-MM groups
        zarr_root = self.D_iceh_zarr
        available_groups = sorted([ p.name for p in zarr_root.glob("????-??") if (zarr_root / p.name / ".zgroup").exists() ])
        if not available_groups:
            raise FileNotFoundError(f"No Zarr groups found in {zarr_root}.")
        # Step 2: Load first and last datasets to get true data bounds
        ds0 = xr.open_zarr(zarr_root, group=available_groups[0], consolidated=False)
        dsN = xr.open_zarr(zarr_root, group=available_groups[-1], consolidated=False)
        available_dt0 = pd.to_datetime(ds0.time.values[0])
        available_dtN = pd.to_datetime(dsN.time.values[-1])
        # Step 3: Clamp user request to data availability
        user_dt0 = max(pd.to_datetime(dt0_str), available_dt0)
        user_dtN = min(pd.to_datetime(dtN_str), available_dtN)
        # Step 4: Select required YYYY-MM groups
        required_groups = [g for g in available_groups
                          if pd.to_datetime(f"{g}-01") <= user_dtN and pd.to_datetime(f"{g}-01") + pd.offsets.MonthEnd(1) >= user_dt0 ]
        self.logger.info(f"Loading Zarr groups between {user_dt0.date()} and {user_dtN.date()}")
        with dask.config.set({'array.slicing.split_large_chunks': True,
                            'array.chunk-size': '256MiB',
                            'optimization.fuse.active': False}):
            ds_list = []
            for g in required_groups:
                self.logger.debug(f"  - opening group {g}")
                ds = xr.open_zarr(zarr_root, group=g, consolidated=False)
                if variables:
                    missing = [v for v in variables if v not in ds]
                    if missing:
                        self.logger.warning(f"  > Skipping {g}, missing: {missing}")
                        continue
                    ds = ds[variables]
                ds_list.append(ds)
            if not ds_list:
                raise ValueError("No datasets to concatenate after filtering.")
            ds_all = xr.concat(ds_list, dim="time", coords="minimal", compat="override")
            ds_all = ds_all.sel(time=slice(user_dt0, user_dtN))
        if slice_hem:
            self.logger.info("  slicing hemisphere")
            ds_all = self.slice_hemisphere(ds_all)
        return ds_all

    @staticmethod
    def _drop_duplicate_coords(ds: xr.Dataset, dim: str = "ni") -> xr.Dataset:
        """
        Drop duplicate coordinate values along `dim` (keeping the first occurrence).
        Useful when concatenating yearly groups that may contain duplicated x-indices.
        """
        if dim in ds.coords:
            _, index = np.unique(ds[dim], return_index=True)
            ds = ds.isel({dim: sorted(index)})
        return ds

    def load_classified_ice(self,
                            sim_name    : str   = None,
                            bin_days    : bool  = True,
                            roll_mean   : bool  = False,
                            ispd_thresh : float = None,
                            ice_type    : str   = None,
                            ivec_type   : str   = None,
                            variables   : list  = None,
                            dt0_str     : str   = None,
                            dtN_str     : str   = None,
                            D_zarr      : str   = None,
                            chunks      : dict  = None,
                            persist     : bool  = False):
        """
        Load classified ice products (daily, binary-day, or rolling) from yearly Zarr groups.

        Parameters
        ----------
        sim_name : str, optional
        bin_days : bool, default True
            If True, load the `_bin.zarr` product; if False and `roll_mean=True`,
            load `_roll.zarr`; else load the daily product.
        roll_mean : bool, default False
        ispd_thresh : float, optional
            Threshold embedded in the Zarr path (defaults to `self.ispd_thresh`).
        ice_type, ivec_type : str, optional
        variables : list[str], optional
            Optional variable subset to select on open.
        dt0_str, dtN_str : str, optional
            Used only to determine which yearly groups to open.
        D_zarr : str | Path, optional
            Override base zarr directory.
        chunks : dict, optional
            Chunking to apply after concat (default `self.CICE_dict["FI_chunks"]`).
        persist : bool, default False
            Persist the concatenated dataset in memory.

        Returns
        -------
        xr.Dataset

        Notes
        -----
        - Drops duplicate `ni` coords per year before concatenation to avoid
          index collisions along x (see `_drop_duplicate_coords()`).
        """
        sim_name     = sim_name     or self.sim_name
        ispd_thresh  = ispd_thresh  or self.ispd_thresh
        ice_type     = ice_type     or self.ice_type
        ivec_type    = ivec_type    or self.ivec_type
        dt0_str      = dt0_str      or self.dt0_str
        dtN_str      = dtN_str      or self.dtN_str
        chunks       = chunks       or self.CICE_dict["FI_chunks"]
        ispd_thresh_str = f"{ispd_thresh:.1e}".replace("e-0", "e-")
        D_class = Path(D_zarr or Path(self.config['D_dict']['AFIM_out'], sim_name, "zarr", f"ispd_thresh_{ispd_thresh_str}"))
        if bin_days:
            zarr_store = D_class / f"{ice_type}_{ivec_type}_bin.zarr"
        elif roll_mean:
            zarr_store = D_class / f"{ice_type}_{ivec_type}_roll.zarr"
        else:
            zarr_store = D_class / f"{ice_type}_{ivec_type}.zarr"
        # === Loop over years, not months ===
        dt0 = datetime.strptime(dt0_str, "%Y-%m-%d")
        dtN = datetime.strptime(dtN_str, "%Y-%m-%d")
        years = list(range(dt0.year, dtN.year + 1))
        datasets = []
        for yr in years:
            try:
                ds_yr = xr.open_zarr(zarr_store, group=str(yr), consolidated=False)
                ds_yr = SeaIceToolbox._drop_duplicate_coords(ds_yr, dim="ni")
                datasets.append(ds_yr)
            except Exception as e:
                self.logger.warning(f"Skipping year {yr}: {e}")
        if not datasets:
            raise FileNotFoundError(f"No valid Zarr groups found for {ice_type} in {zarr_store}")
        ds = xr.concat(datasets, dim="time")
        if variables:
            ds = ds[variables]
        if persist:
            ds = ds.persist()
        return ds 

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

    def compute_doy_climatology(self, da, leap_year=None, time_coord=None):
        """
        Compute day-of-year (DOY) climatology statistics from a time series.

        This method calculates the climatological mean, minimum, maximum, and standard deviation 
        of a given time series DataArray, grouped by day-of-year. The result is returned as a 
        dictionary of Pandas Series indexed by a datetime index constructed using a reference 
        `leap_year`.

        Parameters
        ----------
        da : xarray.DataArray
            The input time series with a time coordinate. Can be daily or sub-daily resolution,
            but must be regular and span multiple years for meaningful climatology.

        leap_year : int, optional
            The reference leap year to use when reconstructing the datetime index for the output.
            Defaults to `self.leap_year` if not provided.

        time_coord : str, optional
            The name of the time coordinate in `da`. Defaults to `self.CICE_dict['time_dim']` if not specified.

        Returns
        -------
        dict
            A dictionary with keys:
            - 'mean' : pd.Series of climatological mean values
            - 'min'  : pd.Series of climatological minimum values
            - 'max'  : pd.Series of climatological maximum values
            - 'std'  : pd.Series of climatological standard deviation

            Each Series is indexed by datetime values (from the specified `leap_year`) corresponding to
            days 1–366.

        Notes
        -----
        - The output includes 366 days if data contains leap years; otherwise, it includes up to 365.
        - The use of a leap year for index construction ensures that the DOY mapping to dates is valid,
        especially for plotting or seasonal alignment.
        - The DataArray is fully loaded into memory before processing.
        """
        leap_year  = leap_year  if leap_year  is not None else self.leap_year
        time_coord = time_coord if time_coord is not None else self.CICE_dict['time_dim']
        da         = da.load()
        df         = pd.DataFrame({"time" : pd.to_datetime(da[time_coord].values),
                                   "data" : da.values})
        df["doy"]  = df["time"].dt.dayofyear
        data_clim  = df.groupby("doy")["data"]
        data_min   = data_clim.min()
        data_max   = data_clim.max()
        data_mean  = data_clim.mean()
        data_std   = data_clim.std()
        t_idx      = pd.to_datetime(data_mean.index - 1, unit="D", origin=pd.Timestamp(f"{leap_year}-01-01"))
        data_min.index = data_max.index = data_mean.index = data_std.index = t_idx
        return {'min'  : data_min,
                'max'  : data_max,
                'std'  : data_std,
                'mean' : data_mean}

    def align_time_coordinate_of_three_arrays(self, ds1, ds2, ds3, time_coord="time"):
        for da in [ds1, ds2, ds3]:
            da[time_coord] = pd.to_datetime(da[time_coord].values).normalize()
        t_common = np.intersect1d(np.intersect1d(ds1[time_coord].values, ds2[time_coord].values), ds3[time_coord].values)
        return ds1.sel(time=t_common), ds2.sel(time=t_common), ds3.sel(time=t_common)

    def cosine_vector_similarity(self, uo, vo, um, vm, eps=1e-12):
        """
        Compute the cosine similarity between two vector fields (e.g., observed vs. modelled velocity vectors).

        This metric quantifies the directional alignment of the two vector fields without considering magnitude. 
        A value of:
        - +1.0 means the vectors point in exactly the same direction,
        -  0.0 means the vectors are orthogonal (90° apart),
        - -1.0 means the vectors point in opposite directions.

        The cosine similarity is computed as the dot product of the two vectors, divided by the product of their magnitudes.

        Parameters
        ----------
        uo, vo : xarray.DataArray
            Components of the observed vector field (e.g., `u` and `v` velocity components) in units of m/s.
        um, vm : xarray.DataArray
            Components of the modelled vector field (same units as `uo`, `vo`).
        eps : float, optional
            Small constant to prevent division by zero in regions where either vector magnitude is near-zero. Default is 1e-12.

        Returns
        -------
        xarray.DataArray
            Cosine similarity between vectors, dimensionless, in the range [-1, 1].

        Notes
        -----
        - NaNs are returned where either the observed or modelled vector magnitude is near-zero.
        - This metric is **scale-invariant** — it compares **direction only**, not speed.
        - Particularly useful for evaluating sea ice drift direction skill, regardless of speed bias.
        """
        dot_prod = um * uo + vm * vo
        obs_mag  = np.sqrt(uo**2 + vo**2)
        mod_mag  = np.sqrt(um**2 + vm**2)
        return dot_prod / xr.where((obs_mag * mod_mag) < eps, np.nan, obs_mag * mod_mag)

    def vector_angle_diff(self, uo, vo, um, vm):
        """
        Compute the signed angular difference (in radians) between two vector fields.

        This metric measures the angle by which the modelled vector differs from the observed vector.
        The difference is returned as a signed value in the range [-π, π], where:
        -  0    indicates perfect alignment,
        - +pi/2 indicates model is rotated 90° counter-clockwise from observation,
        - -pi/2 indicates 90° clockwise rotation,
        - +/-pi indicates vectors are anti-parallel.

        Parameters
        ----------
        uo, vo : xarray.DataArray
            Components of the observed vector field (e.g., `u` and `v` velocity components), in m/s.
        um, vm : xarray.DataArray
            Components of the modelled vector field (in m/s).

        Returns
        -------
        xarray.DataArray
            Signed angular difference (model - obs) in radians, bounded in [-pi, pi].

        Notes
        -----
        - The angular difference is computed using `arctan2` on the vector components.
        - Use `np.rad2deg()` to convert output to degrees, if desired.
        - This metric is useful for directional error analysis in drift or flow fields.
        """
        ang_o = np.arctan2(vo, uo)
        ang_m = np.arctan2(vm, um)
        d     = ang_m - ang_o
        return (d + np.pi) % (2 * np.pi) - np.pi