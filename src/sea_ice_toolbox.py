from __future__ import annotations
import os, sys, logging
import xarray as xr
import pandas as pd
import numpy  as np
from pathlib  import Path
_dev_path = "/home/581/da1339/AFIM/src/AFIM/src"
if os.path.isdir(_dev_path) and _dev_path not in sys.path:
    sys.path.insert(0, _dev_path)
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
    Factory for creating `SeaIceToolbox` instances backed by a shared Dask client.

    `SeaIceToolboxManager` maintains a **class-level** Dask client so that multiple
    toolboxes (e.g., multiple simulations or multiple analyses) reuse a single
    scheduler and worker pool. This is particularly useful on HPC login/compute nodes
    where repeatedly creating local clusters is slow and can fragment resources.

    Parameters
    ----------
    P_log : str or pathlib.Path
        Path to a log file. The created toolboxes attach a file handler to this path.
    n_workers : int, default 4
        Number of Dask workers for the shared LocalCluster.
    n_threads : int, default 1
        Threads per worker.
    mem_lim : str, default "16GB"
        Worker memory limit passed to Dask (e.g., "8GB", "16GB").
    process : bool, default True
        Whether to use separate worker processes (True) or threads-only workers (False).
    D_dask : str, optional
        Local directory for Dask worker spill and temporary files. If None, uses
        `DASK_TEMPORARY_DIRECTORY` or the system temp directory.

    Notes
    -----
    - The shared client is cached on the class (`SeaIceToolboxManager._shared_client`).
    - `shutdown()` closes the shared client and attempts to close file handlers.
    - `get_toolbox(sim_name, **kwargs)` returns a preconfigured `SeaIceToolbox`
    attached to the shared client and log file.

    See Also
    --------
    SeaIceToolbox
        Unified AFIM sea-ice analysis class composed of multiple SeaIce* mixins.
    """
    _shared_client = None

    def __init__(self, P_log,
                 n_workers : int  = 4,
                 n_threads : int  = 1,
                 mem_lim   : str  = "16GB",
                 process   : bool = True,
                 D_dask    : str  = None):
        """
        Create (or reuse) a shared Dask client for AFIM sea-ice analysis.

        If a shared client has not yet been created, this initialiser creates a
        `dask.distributed.LocalCluster` and a corresponding `Client`, storing it at the
        class level so that subsequent instances reuse the same client.

        Parameters
        ----------
        P_log : str or pathlib.Path
            Log file path to pass through to created `SeaIceToolbox` instances.
        n_workers : int, default 4
            Number of Dask workers.
        n_threads : int, default 1
            Threads per worker.
        mem_lim : str, default "16GB"
            Worker memory limit (Dask "memory_limit").
        process : bool, default True
            If True, use multiple processes. If False, use threads-only workers.
        D_dask : str, optional
            Local directory for worker temporary files and spill. If None, uses
            `DASK_TEMPORARY_DIRECTORY` or the system temp directory.

        Notes
        -----
        - This method does not return a `Client`; access it via
        `SeaIceToolboxManager._shared_client`.
        - The shared client is created once per Python process.
        """
        import tempfile
        from dask.distributed  import Client, LocalCluster
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
        """
        Shut down the shared Dask client and close file log handlers.

        Closes the class-level Dask client if it exists and resets the cached reference.
        Also attempts to close any `logging.FileHandler` instances attached to the
        `SeaIceToolbox` logger.

        Notes
        -----
        - Intended for interactive sessions to avoid orphaned local clusters.
        - If multiple loggers/handlers are in use, only file handlers on the named
        logger are closed.
        """
        import logging
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
        """
        Create a `SeaIceToolbox` for a given simulation using the shared Dask client.

        Parameters
        ----------
        sim_name : str
            Simulation name (must correspond to a configured simulation directory).
        **kwargs
            Additional keyword arguments forwarded to `SeaIceToolbox(...)`. Use this to
            override date ranges, thresholds, plotting options, etc.

        Returns
        -------
        SeaIceToolbox
            Toolbox instance bound to the shared Dask client and configured log file.

        Raises
        ------
        ValueError
            If the shared Dask client has not been created successfully.
        """
        return SeaIceToolbox(sim_name = sim_name,
                             client   = SeaIceToolboxManager._shared_client,  # explicit
                             P_log    = self.P_log, **kwargs)

####################################################################################################################

class SeaIceToolbox(SeaIceClassification, SeaIceMetrics, SeaIcePlotter, 
                    SeaIceIcebergs, SeaIceObservations, SeaIceACCESS, SeaIceGridWork, SeaIceRegridder):
    """
    Unified AFIM toolbox for processing and analysing Antarctic sea ice from CICE.

    `SeaIceToolbox` composes several functional mixins into a single interface for
    reading CICE output, classifying fast/pack ice, computing metrics, regridding,
    working with grounded-iceberg masks, and producing plots.

    Inheritance / Composition
    -------------------------
    SeaIceToolbox inherits from the following modules:

    - SeaIceClassification : classification masks and simulation I/O helpers
    - SeaIceMetrics        : time-series and spatial metrics, skill statistics
    - SeaIcePlotter        : PyGMT maps and time series plotting utilities
    - SeaIceIcebergs       : grounded-iceberg thinning/masking and GI datasets
    - SeaIceObservations   : Fraser et al. (2020) and NSIDC observational utilities
    - SeaIceACCESS         : ACCESS-OM related helpers (where applicable)
    - SeaIceGridWork       : grid geometry, landmask application, hemisphere slicing
    - SeaIceRegridder      : B-grid → T-grid and swath/gridded regridding utilities

    Configuration
    -------------
    The toolbox is configured primarily by an AFIM JSON file. The constructor loads
    that JSON, sets commonly used directories (simulation output, zarr, metrics,
    figures), defines hemisphere behaviour, and stores thresholds used throughout
    the pipeline.

    Parameters
    ----------
    P_json : str or pathlib.Path, optional
        Path to the AFIM configuration JSON. If None, a project default path is used.
    P_CICE_grid : str or pathlib.Path, optional
        Override path to the CICE grid file; otherwise uses config entry `CICE_dict['P_G']`.
    sim_name : str
        Simulation name (must exist under the configured AFIM output root).
    dt0_str, dtN_str : str, optional
        Inclusive analysis window bounds in ``YYYY-MM-DD``.
    list_of_BorC2T : list[str], optional
        Speed/vector products to use (e.g., ["Tb"], ["Ta","Tx"], or ["B"]).
    iceh_frequency : {"hourly","daily","monthly","yearly"}, optional
        Which CICE history cadence to use for iceh inputs.
    ice_concentration_threshold : float, optional
        Concentration threshold used for masking / metrics (default from config; often 0.15).
    ice_speed_threshold : float, optional
        Speed threshold (m/s) below which ice is treated as “fast” (default from config).
    ice_type : str or list[str], optional
        Ice classification type(s) to process (e.g., "FI", "PI", "SI").
    mean_period : int, optional
        Rolling mean window length (days) used in some classification products.
    bin_win_days, bin_min_days : int, optional
        Binary-days window length and minimum count used for persistence classification.
    extra_cice_vars : list[str] or bool, optional
        Additional variables to include when loading/processing beyond `cice_vars_reqd`.
    hemisphere : str, optional
        Hemisphere selector. Common aliases are accepted (e.g., "south", "SH", "nh").
    P_log : str or pathlib.Path, optional
        Log file path to attach a `FileHandler` to.
    log_level : int or str, optional
        Python logging level for the toolbox logger.
    dask_memory_limit : str, optional
        Memory limit (informational here unless you create the client externally).
    overwrite_zarr : bool, optional
        If True, overwrite Zarr groups when writing classification/metrics products.
    overwrite_saved_figs : bool, optional
        If True, overwrite existing figure files.
    save_new_figs, show_figs : bool, optional
        Figure saving and interactive display toggles.
    delete_original_cice_iceh_nc : bool, optional
        If True, delete original NetCDF inputs after conversion to monthly Zarr (where implemented).
    client : dask.distributed.Client, optional
        Dask client to use. In this implementation, a client must already exist and be passed in.
    force_recompile_ice_in : bool, default False
        Force regeneration of derived simulation metadata (e.g., parsing/rebuilding ice_in JSON).
    **kwargs
        Additional keyword arguments are attached to the instance so mixins can
        access specialised tunables without changing the constructor signature.

    Attributes Set (high-level)
    ---------------------------
    - Configuration dicts: `self.CICE_dict`, `self.GI_dict`, `self.NSIDC_dict`, etc.
    - Directory paths: `self.D_sim`, `self.D_zarr`, `self.D_metrics`, etc.
    - Hemisphere metadata: `self.hemisphere_dict`, `self.hemisphere`
    - Thresholds/settings: `self.ispd_thresh`, `self.icon_thresh`, `self.mean_period`, etc.
    - State flags: `self.grid_loaded`, `self.reG_weights_defined`, etc.

    Notes
    -----
    - The constructor performs substantial configuration and path normalisation and
    therefore has side effects (file IO, logger configuration).
    - Grids and regridders are not loaded/built until needed (`load_cice_grid`,
    `define_reG_weights`, etc.).
    """
    def summary(self):
        """
        Log a concise summary of key configuration and runtime settings.

        This is primarily a convenience method for sanity checking during interactive
        work and batch runs. It writes configuration metadata to `self.logger`.

        Notes
        -----
        - No values are returned.
        - Assumes `self.logger` has already been configured and that core attributes
        (sim_name, date bounds, thresholds, toggles) are present.
        """
        self.logger.info("--- SeaIceToolbox Summary ---")
        self.logger.info(f"Simulation Name     : {self.sim_name}")
        self.logger.info(f"Analysis Start Date : {self.dt0_str}")
        self.logger.info(f"Analysis End Date   : {self.dtN_str}")
        self.logger.info(f"Speed Threshold     : {self.ispd_thresh:.1e} m/s")
        self.logger.info(f"BorC-regrid Type(s) : {self.BorC2T_type}")
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

    def __init__(self,
                 P_json                      = None,# the configuration file for which there are many dependencies
                                                    # that this toolbox relies upon
                 P_CICE_grid                 = None,# name of the CICE grid file to be used; default is in JSON file
                 sim_name                    = None,# valid name of a model simulation; essentially 'valid' means
                                                    # any name given underneath the directory in the config file
                                                    # named 'AFIM_out'; the sea_ice_model class underneath the hood
                                                    # of this super-class relies on this name for processing and
                                                    # loading of simulation data
                 dt0_str                     = None,# the start period over which many methods underneath use
                                                    # format is YYYY-MM-DD; default is 1993-01-01
                 dtN_str                     = None,# the end period over which many methods underneath use
                                                    # format is YYYY-MM-DD; default is 1999-12-31
                 list_of_BorC2T              = None,# select list of ["B", "Ta", "Tb", "Tc", "Tx"];  must be a list;
                                                    # default ["Tb"]
                 iceh_frequency              = None,# 'hourly', 'daily', 'monthly', 'yearly'
                                                    # defines the history files that will be used by this toolbox
	             ice_concentration_threshold = None,# almost should never be changed from default value of
                                                    # 0.15 (grid cell concentration)
	             ice_speed_threshold         = None,# a significantly important value in the determination
                                                    # and classification (masking) of fast ice; defualt value
                                                    # 5e-4 m/s
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
                 overwrite_saved_figs        = None,# whether or not to overwite saved figures; default is false
                 save_new_figs               = None,# whether or not to write new figures to disk; default is true
                 show_figs                   = None,# whether or not to show/print figures to screen; default is false
                 delete_original_cice_iceh_nc= None,# whether or not to delete the original CICE ice history
                 client                      = None,# dask distributed client, can be externally passed here
                 force_recompile_ice_in      = False,# reinitialise ice_in JSON file; see help self.parse_simulation_metadata()
                 **kwargs):
        """
        Initialise the unified AFIM sea-ice toolbox for a given simulation/config.

        The constructor loads the AFIM JSON configuration file, normalises user-provided
        overrides (simulation name, date bounds, thresholds, hemisphere), defines key
        directories, parses simulation metadata, and initialises the composed mixins.

        Parameters
        ----------
        P_json : str or pathlib.Path, optional
            Path to the AFIM configuration JSON. If None, uses a project default path.
        P_CICE_grid : str or pathlib.Path, optional
            Optional override for the CICE grid file path used by `load_cice_grid`.
        sim_name : str, optional
            Simulation name used to resolve directories under AFIM output root.
        dt0_str, dtN_str : str, optional
            Inclusive analysis window bounds in ``YYYY-MM-DD``.
        list_of_BorC2T : list[str], optional
            Speed/vector product identifiers (e.g., ["Tb"], ["Ta","Tx"], ["B"]).
        iceh_frequency : {"hourly","daily","monthly","yearly"}, optional
            CICE history cadence to use for loading inputs and/or locating Zarr stores.
        ice_concentration_threshold : float, optional
            Concentration threshold used in masking/metrics (often 0.15).
        ice_speed_threshold : float, optional
            Speed threshold (m/s) used to classify fast ice (e.g., 5e-4).
        ice_type : str or list[str], optional
            Ice classification product(s) to work with ("FI", "PI", "SI", etc.).
        mean_period : int, optional
            Rolling mean window length (days) used for smoothed classification products.
        bin_win_days, bin_min_days : int, optional
            Binary-days window length and minimum count for persistence-based masks.
        extra_cice_vars : list[str] or bool, optional
            Additional CICE variables to include beyond `cice_vars_reqd`. If True, use
            config `cice_vars_ext`. If a list, extend required list with that list.
        hemisphere : str, optional
            Hemisphere selector (accepts common aliases).
        P_log : str or pathlib.Path, optional
            Log file path to attach to the toolbox logger.
        log_level : int or str, optional
            Logging verbosity level.
        dask_memory_limit : str, optional
            Informational; the client is expected to be created externally.
        overwrite_zarr : bool, optional
            Overwrite Zarr groups when writing outputs.
        overwrite_saved_figs : bool, optional
            Overwrite existing figure files.
        save_new_figs, show_figs : bool, optional
            Figure saving and interactive display toggles.
        delete_original_cice_iceh_nc : bool, optional
            Delete original NetCDF after conversion to Zarr (where implemented).
        client : dask.distributed.Client, optional
            Dask client for computation. This implementation requires an existing client.
        force_recompile_ice_in : bool, default False
            Force regeneration/reparse of simulation metadata.
        **kwargs
            Additional keyword arguments are attached to the instance and available to mixins.

        Raises
        ------
        FileNotFoundError
            If `P_json` is provided but cannot be opened.
        ValueError
            If a Dask client is not provided (and no manager-created client is present),
            or if hemisphere is invalid.
        KeyError
            If required keys are missing from the JSON configuration.

        Side Effects
        ------------
        - Reads JSON configuration and may open/create the log file path.
        - Configures `self.logger` (stream handler + optional file handler).
        - Defines many path attributes and state flags.
        - Initialises mixin classes.

        Notes
        -----
        - Heavyweight initialisation is intentional; it ensures all mixins share a
        consistent configuration namespace.
        - Grids and regridders are not loaded until requested by other methods.
        """
        import json
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
        self.metrics_name       = self.config.get("metrics_name"      , "mets")
        self.FI_class_types     = self.config.get("FI_class_types"    , {})
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
        self.dt0_str              = dt0_str                      if dt0_str                      is not None else self.config.get('dt0_str', '1993-01-01')
        self.dtN_str              = dtN_str                      if dtN_str                      is not None else self.config.get('dtN_str', '1999-12-31')
        self.ispd_thresh          = ice_speed_threshold          if ice_speed_threshold          is not None else self.config.get('ice_speed_thresh_hi', 5.0e-4)
        self.BorC2T_type          = list_of_BorC2T               if list_of_BorC2T               is not None else self.config.get('BorC2T_type', ['Tb'])
        self.ice_type             = ice_type                     if ice_type                     is not None else self.config.get('ice_type', 'FI')
        self.iceh_freq            = iceh_frequency               if iceh_frequency               is not None else self.config.get('iceh_freq', 'daily')
        self.mean_period          = mean_period                  if mean_period                  is not None else self.config.get('mean_period', 15)
        self.bin_win_days         = bin_win_days                 if bin_win_days                 is not None else self.config.get('bin_win_days', 11)
        self.bin_min_days         = bin_min_days                 if bin_min_days                 is not None else self.config.get('bin_min_days', 9)
        self.icon_thresh          = ice_concentration_threshold  if ice_concentration_threshold  is not None else self.config.get('ice_conc_thresh', 0.15)
        self.overwrite_zarr_group = overwrite_zarr               if overwrite_zarr               is not None else False
        self.ow_fig               = overwrite_saved_figs         if overwrite_saved_figs         is not None else False
        self.save_fig             = save_new_figs                if save_new_figs                is not None else True
        self.show_fig             = show_figs                    if show_figs                    is not None else False
        self.del_org_cice_iceh_nc = delete_original_cice_iceh_nc if delete_original_cice_iceh_nc is not None else False
        hemisphere                = hemisphere                   if hemisphere                   is not None else self.config.get('hemisphere', 'south')
        self.CICE_dict['P_G']     = P_CICE_grid                  if P_CICE_grid                  is not None else self.CICE_dict['P_G']
        self.define_hemisphere(hemisphere)
        self.ispd_thresh_str     = f"{self.ispd_thresh:.1e}".replace("e-0", "e-")
        self.D_sim               = Path(self.D_dict['AFIM_out'], sim_name)
        self.D_iceh              = Path(self.D_sim , 'history', 'daily')
        self.D_zarr              = Path(self.D_sim , 'zarr')
        self.D_graph             = Path(self.config['D_dict']['graph'], 'AFIM')
        self.D_tmp               = Path(self.config['D_dict']['tmp'])
        self.D_metrics           = Path(self.D_zarr, f"ispd_thresh_{self.ispd_thresh_str}", "metrics")
        self.sim_config          = self.parse_simulation_metadata(force_recompile=force_recompile_ice_in)   
        self.valid_BorC2T_types    = self.config.get("valid_BorC2T_types", [])
        self.valid_ice_types     = self.config.get("valid_ice_types", [])
        self.cice_vars_reqd      = self.CICE_dict["cice_vars_reqd"]
        self.spatial_dims        = self.CICE_dict["spatial_dims"]
        self.define_ispd_thresh_dir()
        self.define_ice_class_name()
        self.define_fast_ice_class_name()
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
            if self.GI_thin is not None and self.GI_thin>0 and self.GI_version>0:
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
        self.grid_loaded               = False
        SeaIceClassification.__init__(self, sim_name, **kwargs)
        SeaIceMetrics.__init__(self, **kwargs)
        SeaIcePlotter.__init__(self, **kwargs)
        SeaIceIcebergs.__init__(self, **kwargs)
        SeaIceObservations.__init__(self, **kwargs)
        SeaIceACCESS.__init__(self, **kwargs)
        self.summary()

    def setup_logging(self, logfile=None, log_level=logging.INFO):
        """
        Configure the toolbox logger with a stream handler and optional file handler.

        This method creates (or reuses) a named logger, sets the logging level, and
        ensures that (a) a console StreamHandler exists and (b) exactly one FileHandler
        is attached when `logfile` is provided.

        Parameters
        ----------
        logfile : str or pathlib.Path, optional
            File to write log messages to. If None, only console logging is configured.
        log_level : int, default logging.INFO
            Logging level for the logger and its handlers.

        Notes
        -----
        - Existing FileHandlers are removed and closed before attaching a new one.
        - `self.logger.propagate` is set False to avoid duplicate logs when parent
        loggers are configured elsewhere.
        """
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
        import json, re
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
        Initialise hemisphere configuration used across grid slicing and plotting.

        Accepts common hemisphere aliases (e.g., "south", "SH", "nh") and maps them to
        the internal `self.hemisphere_dict` configuration (slicing indices and labels).
        Also converts `nj_slice` from a (start, stop) tuple into a Python `slice`.

        Parameters
        ----------
        hemisphere : str
            Hemisphere selector. Accepted aliases include:
            - North: "north", "northern", "nh", "n"
            - South: "south", "southern", "sh", "s"

        Sets
        ----
        hemisphere_dict : dict
            Hemisphere metadata dictionary loaded from the JSON config, with `nj_slice`
            stored as a Python `slice`.
        hemisphere : str
            Canonical lowercase hemisphere string used internally.

        Raises
        ------
        ValueError
            If `hemisphere` does not match any accepted alias.

        Notes
        -----
        Downstream methods rely on `self.hemisphere_dict['nj_slice']` for slicing and
        `self.hemisphere_dict['abbreviation']` for path construction.
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
        """
        Define the NetCDF and Zarr directories for CICE ice history ("iceh") inputs.

        Parameters
        ----------
        D_sim : str or pathlib.Path, optional
            Base simulation directory. Defaults to `self.D_sim`.
        iceh_freq : str, optional
            Ice history cadence (e.g., "daily", "monthly"). Defaults to `self.iceh_freq`.

        Sets
        ----
        D_iceh_netcdf : pathlib.Path
            Directory containing NetCDF ice history files.
        D_iceh_zarr : pathlib.Path
            Zarr store path for ice history (e.g., ".../zarr/iceh_daily.zarr").

        Notes
        -----
        This method only defines paths; it does not create directories or validate
        existence.
        """
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
        self.D_ispd_thresh = Path(D_zarr, self.hemisphere_dict['abbreviation'], f"ispd_thresh_{ispd_thresh_str}")
        
    def define_classification_dir(self, ice_type=None, D_zarr=None, ispd_thresh=None):
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
        ice_type        = ice_type    or self.ice_type
        D_zarr          = D_zarr      or self.D_zarr
        ispd_thresh     = ispd_thresh or self.ispd_thresh
        ispd_thresh_str = f"{ispd_thresh:.1e}".replace("e-0", "e-").replace("e+0", "e+")
        if ice_type == "SI":
            self.D_class = Path(D_zarr, self.hemisphere_dict["abbreviation"])
        else:
            self.D_class = Path(D_zarr, self.hemisphere_dict["abbreviation"], f"ispd_thresh_{ispd_thresh_str}")

    def _check_BorC2T_type(self,BorC2T_type):
        if isinstance(BorC2T_type, str):
            BorC2T_type = [BorC2T_type]
        assert all(v in self.valid_BorC2T_types for v in BorC2T_type), f"Invalid BorC2T_type: {BorC2T_type}"

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
    
    def _get_first(self, ds: xr.Dataset, names) -> Optional[xr.DataArray]:
        """
        Return the first matching variable or coordinate from a dataset.

        Parameters
        ----------
        ds : xarray.Dataset
            Dataset to search.
        names : iterable of str
            Candidate names to check in order.

        Returns
        -------
        xarray.DataArray or None
            The first match found in ds.variables or ds.coords, else None.
        """
        for n in names:
            if n in ds.variables:
                return ds[n]
            if n in ds.coords:
                return ds.coords[n]
        return None

    def _norm_list(self, x: Optional[Iterable[str]]) -> Optional[List[str]]:
        """
        Normalize an optional iterable of strings into a clean list.

        Parameters
        ----------
        x : iterable of str or None
            Input strings (e.g., sensors, levels, versions). Elements that are None,
            empty, or whitespace-only are removed.

        Returns
        -------
        list[str] or None
            Cleaned list of non-empty strings, or None if the result is empty.
        """
        if x is None: return None
        out = [str(s).strip() for s in x if s and str(s).strip()]
        return out or None

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
        
    def define_ice_class_name(self, ice_type=None , BorC2T_type=None ):
        """
        Define the classification name string for ice type and vector component type.

        Parameters
        ----------
        ice_type : str or list of str, optional
            Type(s) of ice classification (e.g., 'fast', 'drift'). Defaults to `self.ice_type`.
        BorC2T_type : str or list of str, optional
            Type(s) of ice velocity vector used (e.g., 'B', 'Ta', 'Tb', 'Tx'). Defaults to `self.BorC2T_type`.

        Raises
        ------
        AssertionError
            If any provided `ice_type` or `BorC2T_type` is not in the list of valid types.

        Sets
        ----
        self.ice_class : str
            Combined classification string in the format "{ice_type}_{BorC2T_type}".
        """
        ice_type    = ice_type    or self.ice_type
        BorC2T_type = BorC2T_type or self.BorC2T_type
        self._check_BorC2T_type(BorC2T_type) 
        self._check_ice_type(ice_type)  
        self.ice_class = f"{ice_type}_{BorC2T_type}"
        self.logger.info(f" self.ice_class defined as {self.ice_class}")
        
    def define_fast_ice_class_name(self, BorC2T_type = None , fast_ice_class_method = 'binary-days'):
        """
        Define the canonical fast-ice classification product name.

        This method standardises the naming of fast-ice classification products based
        on the velocity product(s) used (e.g., "Tb") and the classification method
        (e.g., "binary-days"). The resulting name is stored in `self.FI_class`.

        Parameters
        ----------
        BorC2T_type : str or list[str], optional
            Velocity product identifier(s). Defaults to `self.BorC2T_type`.
        fast_ice_class_method : str, default "binary-days"
            Classification method key used to select a suffix from `self.FI_class_types`.

        Sets
        ----
        FI_class : str
            Name string such as "FI_Tb_BIN" (exact suffix depends on config mapping).

        Raises
        ------
        KeyError
            If `fast_ice_class_method` is not present in `self.FI_class_types`.
        ValueError
            If `BorC2T_type` fails validation via `_check_BorC2T_type`.

        Notes
        -----
        The class name is used in output paths and Zarr store naming.
        """
        BorC2T_type = BorC2T_type or self.BorC2T_type
        self._check_BorC2T_type(BorC2T_type)
        FI_class_name = self.FI_class_types[fast_ice_class_method]
        self.FI_class = f"FI_{BorC2T_type}_{FI_class_name}"
        self.logger.info(f" self.FI_class defined as {self.FI_class}")

    def define_ice_speed_name(self, BorC2T_type=None):
        """
        Set the canonical name for the selected ice-speed vector type.

        Parameters
        ----------
        BorC2T_type : str, optional
            One of the valid vector types (e.g., ``"B"``, ``"Ta"``, ``"Tx"``, or
            composites depending on your config). Defaults to ``self.BorC2T_type``.

        Sets
        ----
        self.ispd_name : str
            Name used throughout outputs/paths, formatted as ``f"ispd_{BorC2T_type}"``.

        Raises
        ------
        ValueError
            If `BorC2T_type` is invalid (validated by `_check_BorC2T_type`).
        """
        BorC2T_type = BorC2T_type or self.BorC2T_type
        self._check_BorC2T_type(BorC2T_type) 
        self.ispd_name = f"ispd_{BorC2T_type}"

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

    def _days_in_month(self, y: int, m: int) -> int:
        """
        Return the number of days in a given month of a given year (Gregorian calendar).

        Parameters
        ----------
        y : int
            Year (e.g., 2002).
        m : int
            Month number in [1..12].

        Returns
        -------
        int
            Number of days in the month, accounting for leap years.
        """
        if m in (1,3,5,7,8,10,12): return 31
        if m in (4,6,9,11): return 30
        return 29 if ((y % 4 == 0 and y % 100 != 0) or (y % 400 == 0)) else 28

    def _month_overlap(self, y: int, m: int, t0: pd.Timestamp, t1: pd.Timestamp) -> bool:
        """
        Check whether a given (year, month) overlaps an inclusive UTC time window.

        Parameters
        ----------
        y, m : int
            Year and month to test.
        t0, t1 : pandas.Timestamp
            Start and end timestamps of the desired window (inclusive), interpreted in UTC.

        Returns
        -------
        bool
            True if the month interval intersects [t0, t1], otherwise False.
        """
        first = pd.Timestamp(year=y, month=m, day=1, tz="UTC")
        last  = pd.Timestamp(year=y, month=m, day=self._days_in_month(y,m), tz="UTC")
        return not (last < t0 or first > t1)

    def _parse_yyyymmdd(self, name: str) -> Optional[pd.Timestamp]:
        """
        Parse a YYYYMMDD date token from a filename-like string.

        This looks for an 8-digit date token in the range 2000–2099 and returns it as a
        UTC pandas.Timestamp.

        Parameters
        ----------
        name : str
            Filename (or other string) to search.

        Returns
        -------
        pandas.Timestamp or None
            Parsed UTC timestamp for the detected date, or None if no valid token exists.

        Notes
        -----
        - The regex is constrained to 20xx years by design.
        - Invalid day/month combinations return None.
        """
        m = re.search(r"(?<!\d)(20\d{2})(0[1-9]|1[0-2])(0[1-9]|[12]\d|3[01])(?!\d)", name)
        if not m: return None
        y, mo, d = int(m.group(1)), int(m.group(2)), int(m.group(3))
        try:
            return pd.Timestamp(year=y, month=mo, day=d, tz="UTC")
        except Exception:
            return None

    def _parse_yyyymm(self, name: str) -> Optional[Tuple[int,int]]:
        """
        Parse a YYYYMM month token from a filename-like string.

        Parameters
        ----------
        name : str
            Filename (or other string) to search.

        Returns
        -------
        tuple[int, int] or None
            (year, month) if a YYYYMM token is found (2000–2099), else None.
        """
        m = re.search(r"(?<!\d)(20\d{2})(0[1-9]|1[0-2])(?!\d)", name)
        if not m: return None
        return int(m.group(1)), int(m.group(2))

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
        """
        Compute and log the total on-disk size of a directory tree.

        Parameters
        ----------
        path : str or pathlib.Path
            Directory to scan recursively.

        Notes
        -----
        - Size is computed from file sizes returned by `Path.rglob` and may be slow on
        large filesystems.
        - The result is logged in GiB (1024**3).
        """
        size_gb = sum(f.stat().st_size for f in path.rglob("*") if f.is_file()) / (1024**3)
        self.logger.info(f"Disk-usage (size) of directory {path}: {size_gb:.2f} GB")

    def create_monthly_strings(self, dt0_str=None, dtN_str=None):
        """
        Return sorted unique YYYY-MM strings between two dates (inclusive).

        Parameters
        ----------
        dt0_str, dtN_str : str, optional
            Start/end dates in ``YYYY-MM-DD``. Defaults to `self.dt0_str` and `self.dtN_str`.

        Returns
        -------
        list[str]
            Sorted list of unique month strings, e.g., ["1993-01", "1993-02", ...].

        Notes
        -----
        Uses a daily pandas date range and then de-duplicates by month.
        """
        dt0_str = dt0_str if dt0_str is not None else self.dt0_str
        dtN_str = dtN_str if dtN_str is not None else self.dtN_str
        dts     = pd.date_range(dt0_str, dtN_str, freq="D")
        return sorted(set(dt.strftime("%Y-%m") for dt in dts))

    def create_empty_valid_DS_dictionary(self, valid_zarr_DS_list=None):
        """
        Create a nested dictionary template for collecting datasets by category.

        Parameters
        ----------
        valid_zarr_DS_list : list[str], optional
            List of dataset keys to initialise per outer key. Defaults to `self.valid_ice_types`.

        Returns
        -------
        collections.defaultdict
            A defaultdict where each new outer key maps to a dict of empty lists:
            { <outer_key>: {<ds_key_1>: [], <ds_key_2>: [], ...} }.

        Notes
        -----
        This is a lightweight utility for accumulating per-year/per-month objects before
        concatenation.
        """
        from collections import defaultdict
        valid_DS_list = valid_zarr_DS_list if valid_zarr_DS_list is not None else self.valid_ice_types
        return defaultdict(lambda: {k: [] for k in valid_DS_list})

    def count_zarr_files(self, path):
        """
        Count and log the number of files under a directory tree (e.g., a Zarr store).

        Parameters
        ----------
        path : str or pathlib.Path
            Directory to count files under.

        Notes
        -----
        This counts filesystem entries, not Zarr logical arrays. It is useful for quick
        sanity checks and estimating metadata overhead.
        """
        total_files = sum(len(files) for _, _, files in os.walk(path))
        self.logger.info(f"{path} contains {total_files} files")

    def radians_to_degrees(self, da):
        """
        Convert radians to degrees.

        Parameters
        ----------
        da : array-like
            Values in radians (NumPy array, xarray DataArray, or scalar).

        Returns
        -------
        same type as `da`
            Values converted to degrees.
        """
        return (da * 180) / np.pi

    def get_static_grid(self, var):
        """
        Return a static (time-independent) slice of a grid-like DataArray.

        If the input has a "time" dimension, returns `isel(time=0)`, otherwise returns
        the object unchanged.

        Parameters
        ----------
        var : xarray.DataArray
            DataArray that may or may not have a "time" dimension.

        Returns
        -------
        xarray.DataArray
            Time-sliced or original DataArray.
        """
        if "time" in var.dims:
            return var.isel(time=0)
        return var

    def dict_to_ds(self, data_dict):
        """
        Convert a dictionary of DataArrays into an xarray.Dataset.

        Parameters
        ----------
        data_dict : dict
            Mapping of variable name -> xarray.DataArray (or array-like compatible).

        Returns
        -------
        xarray.Dataset
            Dataset with keys from `data_dict` as variables.
        """
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
        from datetime import datetime, timedelta
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
                                             engine   = netcdf_engine,
                                             parallel = True,
                                             combine  = "by_coords",
                                             cache    = True,
                                             chunks   = {})
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

        groups = sorted(
            [p.name for p in zarr_root.glob("????-??")
            if (zarr_root / p.name / ".zgroup").exists()]
        )
        if not groups:
            raise FileNotFoundError(f"No Zarr groups found in {zarr_root}.")

        # Use consolidated metadata if present; else fall back
        consolidated = (zarr_root / ".zmetadata").exists()

        ds0 = xr.open_zarr(zarr_root, group=groups[0], consolidated=consolidated)
        dsN = xr.open_zarr(zarr_root, group=groups[-1], consolidated=consolidated)

        available_dt0 = pd.to_datetime(ds0.time.values[0])
        available_dtN = pd.to_datetime(dsN.time.values[-1])

        self._zarr_group_index[key] = (groups, available_dt0, available_dtN, consolidated)
        self.logger.info(
            f"Indexed Zarr root {zarr_root} (groups={len(groups)}, "
            f"bounds={available_dt0.date()}..{available_dtN.date()}, consolidated={consolidated})"
        )
        return self._zarr_group_index[key]

    def load_cice_zarr(self,
                    sim_name=None,
                    dt0_str=None,
                    dtN_str=None,
                    D_alt=None,
                    variables=None,
                    slice_hem=False):
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

        required_groups = [
            g for g in available_groups
            if (pd.to_datetime(f"{g}-01") <= user_dtN)
            and (pd.to_datetime(f"{g}-01") + pd.offsets.MonthEnd(1) >= user_dt0)
        ]

        self.logger.info(f"Loading Zarr groups between {user_dt0.date()} and {user_dtN.date()} ({len(required_groups)} groups)")

        with dask.config.set({
            "array.slicing.split_large_chunks": True,
            "array.chunk-size": "256MiB",
            # Consider turning fuse back on unless you have a specific reason to disable it:
            # "optimization.fuse.active": True,
            "optimization.fuse.active": False,
        }):
            ds_list = []
            for g in required_groups:
                self.logger.debug(f"  - opening group {g}")
                ds = xr.open_zarr(zarr_root, group=g, consolidated=consolidated)

                if variables:
                    missing = [v for v in variables if v not in ds]
                    if missing:
                        self.logger.warning(f"  > Skipping {g}, missing: {missing}")
                        continue

                    # Much clearer logging:
                    self.logger.info(f"  [{g}] selecting variables: {variables}")
                    ds = ds[variables]

                # Optional but often beneficial: crop each group before concat
                ds = ds.sel(time=slice(user_dt0, user_dtN))
                ds_list.append(ds)

            if not ds_list:
                raise ValueError("No datasets to concatenate after filtering.")

            ds_all = xr.concat(ds_list, dim="time", coords="minimal", compat="override")

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
                            BorC2T_type    : str   = None,
                            variables   : list  = None,
                            dt0_str     : str   = None,
                            dtN_str     : str   = None,
                            D_zarr      : str   = None,
                            chunks      : dict  = None,
                            persist     : bool  = False):
        """
        Load classified ice products from yearly-grouped Zarr stores.

        Loads fast-ice/pack-ice/sea-ice classification products written as Zarr stores
        containing per-year groups. Supports daily, rolling-mean, and binary-days
        products via store naming conventions.

        Parameters
        ----------
        sim_name : str, optional
            Simulation name. Defaults to `self.sim_name`.
        bin_days : bool, default True
            If True, load the binary-days product (`*_bin.zarr`).
        roll_mean : bool, default False
            If True (and `bin_days` is False), load the rolling-mean product (`*_roll.zarr`).
        ispd_thresh : float, optional
            Speed threshold embedded in the Zarr path; defaults to `self.ispd_thresh`.
        ice_type : str, optional
            Ice type identifier ('FI', 'PI', or 'SI'). Defaults to `self.ice_type`.
        BorC2T_type : str, optional
            Velocity product identifier used in the store name (e.g., 'Ta', 'Tb', 'Tx', 'BT').
        variables : list[str], optional
            Variable subset to select after concatenation.
        dt0_str, dtN_str : str, optional
            Requested date window. Used here only to select which years to attempt to open.
        D_zarr : str or pathlib.Path, optional
            Override base Zarr directory.
        chunks : dict, optional
            Chunking mapping applied after concatenation (implementation-dependent).
        persist : bool, default False
            If True, Dask persist the concatenated dataset in memory.

        Returns
        -------
        xarray.Dataset
            Concatenated dataset across requested years.

        Raises
        ------
        FileNotFoundError
            If no valid year groups could be opened from the target store.

        Notes
        -----
        - Uses `_drop_duplicate_coords` (class-level helper) to avoid coordinate collisions
        during concat.
        - This method does not currently crop the returned dataset to dt0/dtN; caller may
        apply `sel(time=slice(...))` if needed.
        """
        from datetime import datetime
        sim_name     = sim_name     or self.sim_name
        ispd_thresh  = ispd_thresh  or self.ispd_thresh
        ice_type     = ice_type     or self.ice_type
        BorC2T_type  = BorC2T_type  or f"{''.join(self.BorC2T_type)}"
        dt0_str      = dt0_str      or self.dt0_str
        dtN_str      = dtN_str      or self.dtN_str
        chunks       = chunks       or self.CICE_dict["FI_chunks"]
        ispd_thresh_str = f"{ispd_thresh:.1e}".replace("e-0", "e-")
        if ice_type=='SI':
            D_class    = Path(D_zarr or Path(self.D_zarr, self.hemisphere_dict['abbreviation']))
            zarr_store = D_class / "SI.zarr"
        else:
            D_class = Path(D_zarr or Path(self.D_zarr, self.hemisphere_dict['abbreviation'], f"ispd_thresh_{ispd_thresh_str}"))
            if bin_days:
                zarr_store = D_class / f"{ice_type}_{BorC2T_type}_bin.zarr"
            elif roll_mean:
                zarr_store = D_class / f"{ice_type}_{BorC2T_type}_roll.zarr"
            else:
                zarr_store = D_class / f"{ice_type}_{BorC2T_type}.zarr"
        # === Loop over years, not months ===
        dt0      = datetime.strptime(dt0_str, "%Y-%m-%d")
        dtN      = datetime.strptime(dtN_str, "%Y-%m-%d")
        years    = list(range(dt0.year, dtN.year + 1))
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
        ds = xr.concat(datasets, dim="time", coords='minimal')
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
    
    def _region_mask(self, lon : xr.DataArray, lat : xr.DataArray, geo_reg : tuple[float, float, float, float], *,
                     right_open: bool = True) -> xr.DataArray:
        """
        Geographic mask for [lon_min, lon_max, lat_min, lat_max].

        right_open=True makes lon_max exclusive (helps avoid double-counting at boundaries).
        """
        lon_min, lon_max, lat_min, lat_max = geo_reg
        m = (lon >= lon_min) & (lat >= lat_min) & (lat <= lat_max)
        if right_open:
            m = m & (lon < lon_max)
        else:
            m = m & (lon <= lon_max)
        return m

    ##########################################################################################
    #                              FAST ICE SPECIFIC FUNCTIONS                               #
    ##########################################################################################
    def coarsen_and_align_simulated_FI_to_observed_FI(self, sim_ds, obs_ds, doy_vals=None, method="mean", obs_time_coord="time"):
        """
        Aggregate daily simulated fields onto observational fast-ice time windows.

        This method coarsens daily model output into the same temporal windows used by
        the AF2020 observational product (typically 24 bins defined by day-of-year start
        values). For each year and each DOY window, it aggregates model data over the
        corresponding date range and writes one output timestep per window.

        Parameters
        ----------
        sim_ds : xarray.Dataset
            Daily model output with a "time" coordinate.
        obs_ds : xarray.Dataset
            Observational dataset providing the target window definitions via its time
            coordinate (default name is `obs_time_coord`).
        doy_vals : list[int], optional
            Day-of-year start values defining observation windows. If None, uses
            `self.AF_FI_dict["DOY_vals"]`.
        method : {"mean","median"}, default "mean"
            Aggregation method to apply within each window.
        obs_time_coord : str, default "time"
            Name of the time coordinate in `obs_ds` used as the output time dimension.

        Returns
        -------
        xarray.Dataset
            Aggregated dataset with one timestep per observation window (per year),
            concatenated along `obs_time_coord`.

        Raises
        ------
        ValueError
            If no windows intersect the simulation data or `method` is unsupported.

        Notes
        -----
        - Window end dates are computed from successive DOY boundaries. The final window
        extends to the end of the year (accounting for leap years).
        - The output time stamps correspond to the start date of each window.
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

    def compute_fipdiff_stats_weighted(self, FIP: xr.Dataset,
                                       PIXEL                 : float = 5000.0,
                                       regs_dict             : dict | None = None,
                                       gi_mask_da            : xr.DataArray | None = None,
                                       clip01                : bool = True,
                                       threecat_percent_only : bool = True):
        regs_dict = regs_dict if regs_dict is not None else self.Ant_8sectors  # adjust name if needed
        for k in ("mod", "obs", "lon", "lat"):
            if k not in FIP:
                raise ValueError("FIP must contain 'mod', 'obs', 'lon', 'lat'.")
        mod = FIP["mod"].astype("float64")
        obs = FIP["obs"].astype("float64")
        if clip01:
            mod = mod.clip(0.0, 1.0)
            obs = obs.clip(0.0, 1.0)
        overlap      = xr.apply_ufunc(np.minimum, mod, obs)
        model_excess = xr.apply_ufunc(np.maximum, mod - obs, 0.0)
        obs_excess   = xr.apply_ufunc(np.maximum, obs - mod, 0.0)
        lon_m = self.normalise_longitudes(FIP["lon"], to="-180-180")
        lat   = FIP["lat"]
        cell_area_km2 = (float(PIXEL) * float(PIXEL)) / 1e6
        if gi_mask_da is not None:
            valid = (~gi_mask_da.astype(bool).broadcast_like(mod))
        else:
            valid = xr.ones_like(mod, dtype=bool)
        # loop over regional dictionary terms to create rows 
        rows = []
        for rname, reg_dict in regs_dict.items():
            geo_reg = reg_dict["plot_region"]  # must be (lon_min, lon_max, lat_min, lat_max) in -180..180
            R = self._region_mask(lon_m, lat, geo_reg) & valid
            A_km2 = float((overlap      * R).sum()) * cell_area_km2
            M_km2 = float((model_excess * R).sum()) * cell_area_km2
            O_km2 = float((obs_excess   * R).sum()) * cell_area_km2
            tot   = A_km2 + M_km2 + O_km2
            pct_A = 100.0 * A_km2 / tot if tot > 0 else np.nan
            pct_M = 100.0 * M_km2 / tot if tot > 0 else np.nan
            pct_O = 100.0 * O_km2 / tot if tot > 0 else np.nan
            MOD_km2 = float((mod * R).sum()) * cell_area_km2
            OBS_km2 = float((obs * R).sum()) * cell_area_km2
            pct_A_of_MOD = 100.0 * A_km2 / MOD_km2 if MOD_km2 > 0 else np.nan
            pct_M_of_MOD = 100.0 * M_km2 / MOD_km2 if MOD_km2 > 0 else np.nan
            pct_A_of_OBS = 100.0 * A_km2 / OBS_km2 if OBS_km2 > 0 else np.nan
            pct_O_of_OBS = 100.0 * O_km2 / OBS_km2 if OBS_km2 > 0 else np.nan
            rows.append(dict(region                    = rname,
                             pct_agreement             = pct_A,
                             pct_model_dominant        = pct_M,
                             pct_observation_dominant  = pct_O,
                             model_FIA_km2             = MOD_km2,
                             obs_FIA_km2               = OBS_km2,
                             agreement_overlap_km2     = A_km2,
                             model_excess_km2          = M_km2,
                             obs_excess_km2            = O_km2,
                             agreement_pct_of_model    = pct_A_of_MOD,
                             model_excess_pct_of_model = pct_M_of_MOD,
                             agreement_pct_of_obs      = pct_A_of_OBS,
                             obs_excess_pct_of_obs     = pct_O_of_OBS))
        # create dataframe from dictionary
        df      = pd.DataFrame(rows).set_index("region")
        A_sum   = df["agreement_overlap_km2"].sum()
        M_sum   = df["model_excess_km2"].sum()
        O_sum   = df["obs_excess_km2"].sum()
        tot_sum = A_sum + M_sum + O_sum
        # now sum together for pan-Ant. ... but what about NH??
        ant_row = dict(pct_agreement             = 100.0 * A_sum / tot_sum if tot_sum > 0 else np.nan,
                       pct_model_dominant        = 100.0 * M_sum / tot_sum if tot_sum > 0 else np.nan,
                       pct_observation_dominant  = 100.0 * O_sum / tot_sum if tot_sum > 0 else np.nan,
                       model_FIA_km2             = df["model_FIA_km2"].sum(),
                       obs_FIA_km2               = df["obs_FIA_km2"].sum(),
                       agreement_overlap_km2     = A_sum,
                       model_excess_km2          = M_sum,
                       obs_excess_km2            = O_sum,
                       agreement_pct_of_model    = np.nan,
                       model_excess_pct_of_model = np.nan,
                       agreement_pct_of_obs      = np.nan,
                       obs_excess_pct_of_obs     = np.nan)
        df_all = pd.concat([df, pd.DataFrame(ant_row, index=["ANT"])])
        if threecat_percent_only:
            return df_all[["pct_agreement", "pct_model_dominant", "pct_observation_dominant"]]
        else:
            return df_all[["model_FIA_km2", "agreement_overlap_km2", "model_excess_km2",
                           "agreement_pct_of_model", "model_excess_pct_of_model",
                           "obs_FIA_km2", "agreement_pct_of_obs", "obs_excess_km2", "obs_excess_pct_of_obs",
                           "pct_agreement", "pct_model_dominant", "pct_observation_dominant"]]
            
    def define_FIP_diff_names(self, 
                     AF2020            = False,
                     variable_name     = None,
                     lon_coord_name    = None,
                     lat_coord_name    = None,
                     time_coord_name   = None,
                     spatial_dim_names = None,
                     thresh_val        = None):
        """
        Standardise coordinate/variable naming for fast-ice persistence diagnostics.

        This is a lightweight helper that returns a consistent naming dictionary used by
        fast-ice mask and coordinate extraction routines for either:
          - AF2020 observational datasets (AF2020=True), or
          - model/CICE datasets (AF2020=False).

        Parameters
        ----------
        AF2020 : bool, default False
            If True, derive names and threshold from `self.AF_FI_dict`; otherwise use
            CICE/model conventions.
        variable_name, lon_coord_name, lat_coord_name, time_coord_name : str, optional
            Override defaults for variable and coordinate names.
        spatial_dim_names : tuple[str, str], optional
            Override spatial dim names. Defaults to `self.CICE_dict["spatial_dims"]`.
        thresh_val : float, optional
            Override the threshold value used to define "fast ice" presence.

        Returns
        -------
        dict
            Keys:
            - spatial_dims, time_dim, variable, latitude, longitude, thresh_val
        """
        name_space        = {}
        spatial_dim_names = spatial_dim_names if spatial_dim_names is not None else self.CICE_dict['spatial_dims']
        if AF2020:
            name_var   = variable_name   if variable_name   is not None else self.AF_FI_dict["variable_name"]
            name_lat   = lat_coord_name  if lat_coord_name  is not None else self.AF_FI_dict["lat_coord_name"]
            name_lon   = lon_coord_name  if lon_coord_name  is not None else self.AF_FI_dict["lon_coord_name"]
            name_dt    = time_coord_name if time_coord_name is not None else self.AF_FI_dict["time_coord_name"]
            thresh_val = thresh_val      if thresh_val      is not None else self.AF_FI_dict["threshold_value"]
        else:    
            name_var   = variable_name   if variable_name   is not None else "aice"
            name_lat   = lat_coord_name  if lat_coord_name  is not None else self.CICE_dict["tcoord_names"][1]
            name_lon   = lon_coord_name  if lon_coord_name  is not None else self.CICE_dict["tcoord_names"][0]
            name_dt    = time_coord_name if time_coord_name is not None else self.CICE_dict["time_dim"]
            thresh_val = thresh_val      if thresh_val      is not None else self.icon_thresh
        return {'spatial_dims': spatial_dim_names,
                'time_dim'    : name_dt,
                'variable'    : name_var,
                'latitude'    : name_lat,
                'longitude'   : name_lon,
                'thresh_val'  : thresh_val}

    def define_fast_ice_coordinates(self, ds,
                                    AF2020            = False,
                                    variable_name     = None,
                                    lon_coord_name    = None,
                                    lat_coord_name    = None,
                                    time_coord_name   = None,
                                    spatial_dim_names = None,
                                    thresh_val        = None):
        """
        Extract lat/lon/time arrays and standardised naming for fast-ice mask construction.

        Parameters
        ----------
        ds : xr.Dataset
            Dataset containing the required variable and coordinate fields.
        AF2020 : bool, default False
            Switch between AF2020 and CICE/model conventions.
        (other parameters)
            Passed through to `define_FIP_diff_names(...)`.

        Returns
        -------
        dict
            Keys:
            - datetimes : 1D datetime array
            - latitudes : 2D latitude array
            - longitudes : 2D longitude array
            - names : the naming dictionary returned by `define_FIP_diff_names`
        """
        names = self.define_FIP_diff_names(AF2020            = AF2020,
                                           variable_name     = variable_name,
                                           lon_coord_name    = lon_coord_name,
                                           lat_coord_name    = lat_coord_name,
                                           time_coord_name   = time_coord_name,
                                           spatial_dim_names = spatial_dim_names,
                                           thresh_val        = thresh_val)
        if "time" in ds[names['latitude']].dims:
            lat = ds[names['latitude']].isel(time=0).values
        else:
            lat = ds[names['latitude']].values
        if "time" in ds[names['longitude']].dims:
            lon = ds[names['longitude']].isel(time=0).values
        else:
            lon = ds[names['longitude']].values
        dt = ds[names['time_dim']].values
        return {'datetimes'  : dt,
                'latitudes'  : lat,
                'longitudes' : lon,
                'names'      : names}

    def define_fast_ice_mask_da(self, ds,
                                AF2020            = False,
                                variable_name     = None,
                                lon_coord_name    = None,
                                lat_coord_name    = None,
                                time_coord_name   = None,
                                spatial_dim_names = None,
                                thresh_val        = None):
        """
        Construct a binary fast-ice mask DataArray from a dataset.

        A grid cell is classified as fast ice when:
            ds[variable] >= thresh_val
        where names/threshold are determined by `define_FIP_diff_names()`.

        Parameters
        ----------
        ds : xr.Dataset
            Input dataset.
        AF2020 : bool, default False
            If True, use AF2020 conventions; else CICE/model.
        (other parameters)
            Passed through to `define_fast_ice_coordinates(...)`.

        Returns
        -------
        xr.DataArray
            Binary mask with dims given by `self.CICE_dict['three_dims']` and coordinates
            attached using the configured lat/lon/time coordinate names.

        Notes
        -----
        - This routine assumes the mask should always be expressed on the CICE naming
          conventions (coords from `self.CICE_dict`), even when sourced from AF2020.
        """
        coords  = self.define_fast_ice_coordinates(ds,
                                                    AF2020            = AF2020,
                                                    variable_name     = variable_name,
                                                    lon_coord_name    = lon_coord_name,
                                                    lat_coord_name    = lat_coord_name,
                                                    time_coord_name   = time_coord_name,
                                                    spatial_dim_names = spatial_dim_names,
                                                    thresh_val        = thresh_val)
        FI_mask = xr.where(ds[coords['names']['variable']] >= coords['names']['thresh_val'], 1.0, 0.0)
        return xr.DataArray(FI_mask,
                            dims   = self.CICE_dict['three_dims'],
                            coords = {self.CICE_dict['lat_coord_name'] : (self.CICE_dict['spatial_dims'], coords['latitudes']),
                                      self.CICE_dict['lon_coord_name'] : (self.CICE_dict['spatial_dims'], coords['longitudes']),
                                      self.CICE_dict['time_dim']       : (self.CICE_dict['time_dim']    , coords['datetimes'])})

    def define_fast_ice_persistence_da(self, fip, reG=False):
        """
        Wrap a persistence-like array into a consistently-labelled 2D DataArray.

        Parameters
        ----------
        fip : xr.DataArray
            Persistence field on either the native grid or a regridded regular grid.
        reG : bool, default False
            If True, attach coordinates from `define_regular_G(...)` (0.15° global SH);
            otherwise use fip's existing lat/lon coordinates.

        Returns
        -------
        xr.DataArray
            2D float32 DataArray with spatial dims `self.CICE_dict['spatial_dims']` and
            lat/lon coordinates attached.
        """
        if reG:
            G_dst = self.define_regular_G(0.15, region=[0,360,-90,0], spatial_dim_names=self.CICE_dict['spatial_dims'])
            lats  = G_dst[self.CICE_dict['lat_coord_name']].values
            lons  = G_dst[self.CICE_dict['lon_coord_name']].values
        else:
            lats = fip[self.CICE_dict['lat_coord_name']].values
            lons = fip[self.CICE_dict['lon_coord_name']].values
        return xr.DataArray(fip.astype('float32').values,
                            dims   = self.CICE_dict['spatial_dims'],
                            coords = {self.CICE_dict['lat_coord_name'] : (self.CICE_dict['spatial_dims'], lats),
                                    self.CICE_dict['lon_coord_name'] : (self.CICE_dict['spatial_dims'], lons)})

    def compute_fip_weight(self,
                           FIP   : xr.Dataset, 
                           mode  : str = "max",
                           t     : float = 0.1,
                           gamma : float = 1.0) -> xr.DataArray:
        """
        Compute an opacity/weight field from fast-ice persistence coverage.

        This is used for visualisation or for weighting a misfit field, producing
        values in [0, 1] where:
          - 0 => fully transparent / no weight (below threshold)
          - 1 => fully opaque / full weight (high persistence coverage)

        Parameters
        ----------
        FIP : xr.Dataset
            Must include `FIP["obs"]` (observed persistence fraction). May include
            `FIP["mod"]` (model persistence fraction).
        mode : {"max","mean","prod","obs"}, default "max"
            Coverage operator used before thresholding:
              - "max"  : max(obs, mod)
              - "mean" : 0.5*(obs + mod)
              - "prod" : sqrt(obs*mod) (geometric mean)
              - any other value falls back to obs only
        t : float, default 0.1
            Soft threshold. Values below t map to 0; values linearly ramp to 1 by 1.0.
        gamma : float, default 1.0
            Contrast shaping exponent (w ** gamma). gamma>1 increases contrast.

        Returns
        -------
        xr.DataArray
            Weight field in [0, 1], named "diff_weight", with descriptive attributes.

        Raises
        ------
        KeyError
            If FIP does not contain an "obs" variable.

        Notes
        -----
        - This routine is agnostic to the grid; it simply transforms persistence fractions.
        """
        if "obs" not in FIP:
            raise KeyError("FIP['obs'] is required")
        obs = FIP["obs"]
        mod = FIP.get("mod", None)
        if mode == "max" and mod is not None:
            cov = xr.apply_ufunc(np.maximum, obs, mod)
        elif mode == "mean" and mod is not None:
            cov = 0.5 * (obs + mod)
        elif mode == "prod" and mod is not None:
            cov = xr.apply_ufunc(np.sqrt, obs * mod)  # geometric mean
        else:
            cov = obs  # fall back to obs only
        # Soft ramp: 0 below t, then linear to 1 at 1.0
        w = (cov - t) / (1 - t)
        w = w.clip(0, 1)
        # Optional contrast shaping
        if gamma != 1.0:
            w = w ** gamma
        w.name = "diff_weight"
        w.attrs.update(dict(long_name="opacity weight", comment=f"weight ~ {mode}(obs,mod), threshold={t}, gamma={gamma}"))
        return w
            