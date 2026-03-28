from __future__ import annotations
import os, sys, logging, warnings, shutil
import xarray as xr
import pandas as pd
import numpy  as np
from contextlib import contextmanager
from pathlib    import Path
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
from sea_ice_fast           import SeaIceFast
from sea_ice_cice           import SeaIceCICE
from sea_ice_waves          import SeaIceWaves

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
                    SeaIceIcebergs, SeaIceObservations, SeaIceACCESS,
                    SeaIceGridWork, SeaIceRegridder, SeaIceFast, SeaIceCICE,
                    SeaIceWaves):
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
        self.logger.info(f"grid file           : {self.CICE_dict['P_G']}")
        self.logger.info(f"landmask file       : {self.P_KMT_org}")
        self.logger.info(f"Using GI?           : {self.use_gi}")
        if self.use_gi:
            self.logger.info(f"modified landmask file: {self.P_KMT_mod}")
        self.logger.info(f"Speed Threshold     : {self.ispd_thresh:.1e} m/s")
        self.logger.info(f"BorC-regrid Type(s) : {self.BorC2T_type}")
        self.logger.info(f"Ice Type(s)         : {self.ice_type}")
        self.logger.info(f"Mean Period         : {self.mean_period} days")
        self.logger.info(f"Binary-days Window  : {self.bin_win_days} days")
        self.logger.info(f"Binary-days Min-Days: {self.bin_min_days}")
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
	             P_log                       = None ,# the log file to send print statements to
                 log_level                   = None ,# the logging level (see python logging doc for more info)
                 dask_memory_limit           = None ,# provide the memory limit to dask, default is 16GB
                 overwrite_zarr              = False,# whether or not to overwrite a zarr; default is false
                 overwrite_saved_figs        = False,# whether or not to overwite saved figures; default is false
                 save_new_figs               = True ,# whether or not to write new figures to disk; default is true
                 show_figs                   = False,# whether or not to show/print figures to screen; default is false
                 delete_original_cice_iceh_nc= False,# whether or not to delete the original CICE ice history
                 client                      = None ,# dask distributed client, can be externally passed here
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
        self.methods_init_executed = {}
        self.sim_name = sim_name if sim_name is not None else 'test'
        if P_json is None:
            P_json = '/home/581/da1339/AFIM/src/AFIM/src/JSONs/sea_ice_config.json'
        with open(P_json, 'r') as f:
            self.config = json.load(f)
        hemisphere  = hemisphere or self.config.get('hemisphere', 'south')
        self.D_dict = self.config.get('D_dict', {})
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
        self.class_types_dict     = self.config.get("class_types_dict"  , {})
        self.CICE_dict            = self.config.get("CICE_dict"         , {})
        self.GI_dict              = self.config.get('GI_dict'           , {})
        self.Waves_dict           = self.config.get('Waves_dict'        , {})
        self.NSIDC_dict           = self.config.get('NSIDC_dict'        , {})
        self.BAS_dict             = self.config.get('BAS_dict'          , {})
        self.AF_FI_dict           = self.config.get("AF_FI_dict"        , {})
        self.Sea_Ice_Obs_dict     = self.config.get("Sea_Ice_Obs_dict"  , {})
        self.AOM2_dict            = self.config.get("AOM2_dict"         , {})
        self.MOM_dict             = self.config.get("MOM_dict"          , {})
        self.ERA5_dict            = self.config.get("ERA5_dict"         , {})
        self.ORAS_dict            = self.config.get("ORAS_dict"         , {})
        self.plot_var_dict        = self.config.get("plot_var_dict"     , {})
        self.hemispheres_dict     = self.config.get("hemispheres_dict"  , {})
        self.Ant_8sectors         = self.config.get('Ant_8sectors'      , {})
        self.Ant_2sectors         = self.config.get('Ant_2sectors'      , {})
        self.pygmt_dict           = self.config.get("pygmt_dict"        , {})
        self.pygmt_FIA_dict       = self.config.get('pygmt_FIA_dict'    , {})
        self.pygmt_FI_panel       = self.config.get('pygmt_FI_panel'    , {})
        self.dt0_str              = dt0_str                      or self.config.get('dt0_str', '1993-01-01')
        self.dtN_str              = dtN_str                      or self.config.get('dtN_str', '1999-12-31')
        self.ispd_thresh          = ice_speed_threshold          or self.config.get('ice_speed_thresh_hi', 5.0e-4)
        self.BorC2T_type          = list_of_BorC2T               or self.config.get('BorC2T_type', ['Tb'])
        self.ice_type             = ice_type                     or self.config.get('ice_type', 'FI')
        self.iceh_freq            = iceh_frequency               or self.config.get('iceh_freq', 'daily')
        self.mean_period          = mean_period                  or self.config.get('mean_period', 15)
        self.bin_win_days         = bin_win_days                 or self.config.get('bin_win_days', 11)
        self.bin_min_days         = bin_min_days                 or self.config.get('bin_min_days', 9)
        self.icon_thresh          = ice_concentration_threshold  or self.config.get('ice_conc_thresh', 0.15)
        self.CICE_dict['P_G']     = P_CICE_grid                  or self.CICE_dict['P_G']
        self.overwrite_zarr_group = overwrite_zarr
        self.ow_fig               = overwrite_saved_figs
        self.save_fig             = save_new_figs
        self.show_fig             = show_figs
        self.del_org_cice_iceh_nc = delete_original_cice_iceh_nc
        self.leap_year            = self.config.get("leap_year"         , 1996)
        self.metrics_name         = self.config.get("metrics_name"      , "mets")
        self.valid_BorC2T_types   = self.config.get("valid_BorC2T_types", [])
        self.valid_ice_types      = self.config.get("valid_ice_types"   , [])
        self.cice_vars_reqd       = self.CICE_dict["cice_vars_reqd"]
        self.spatial_dims         = self.CICE_dict["spatial_dims"]
        self.define_hemisphere(hemisphere)
        self.define_toolbox_paths()
        self.parse_simulation_metadata(force_recompile=force_recompile_ice_in)
        if extra_cice_vars is not None:
            if extra_cice_vars:
                self.cice_var_list = self.cice_vars_reqd + self.CICE_dict["cice_vars_ext"]
            else:
                self.cice_var_list = self.cice_vars_reqd + extra_cice_vars
        else:
            self.cice_var_list = self.cice_vars_reqd
        self.FIC_scale = self.config.get('FIC_scale', 1e9)
        self.SIC_scale = self.config.get('SIC_scale', 1e12)
        if self.CICE_dict['coupled']:
            self.P_KMT_org = Path(self.CICE_dict['P_KMT'])
        else:
            self.P_KMT_org = Path(self.GI_dict["D_GI_thin"],self.GI_dict['KMT_org_fmt'])
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
                self.P_KMT_mod = self.P_KMT_org
                self.use_gi    = False
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

    #######################################################################################################
    #########################################        LOGGING          #####################################
    #######################################################################################################
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

    @contextmanager
    def _suppress_large_graph_warning(self):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore",
                                    message  = r"Sending large graph of size.*",
                                    category = UserWarning,
                                    module   = r"distributed\.client")
            yield

    def _method_name(self):
        import inspect
        return inspect.currentframe().f_back.f_code.co_name

    #######################################################################################################
    #########################################    TRIGONOMETRY       #######################################
    #######################################################################################################
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

    #######################################################################################################
    ##################################### SEA-ICE ANALYSIS CONFIG #########################################
    #######################################################################################################
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

    def _check_BorC2T_type(self,BorC2T_type):
        if isinstance(BorC2T_type, str):
            BorC2T_type = [BorC2T_type]
        assert all(v in self.valid_BorC2T_types for v in BorC2T_type), f"Invalid BorC2T_type: {BorC2T_type}"

    def _check_ice_type(self,ice_type):
        if isinstance(ice_type, str):
            ice_type = [ice_type]
        assert all(v in self.valid_ice_types for v in ice_type), f"Invalid ice_type: {ice_type}"

    ##########################################################################################################
    ################################## CICE MODEL CONFIGURATION/DIAGNOSTIRICS ################################
    ##########################################################################################################
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
        self.logger.info(f"reading {self.P_ice_diag} to construct {self.P_sim_cfg}")
        if self.P_sim_cfg.exists() and not force_recompile:
            with open(self.P_sim_cfg, "r") as f:
                self.sim_config = json.load(f)
                return True
        PARAM_KEYS = ["dt", "ndtd", "ndte", "kdyn", "revised_evp", "e_yieldcurve", "e_plasticpot",
                      "Ktens", "kstrength", "Pstar", "Cstar", "Cf", "visc_method", "kmt_file"]
        PATTERNS = {key: re.compile(rf"{key}\s*=\s*(.+?)\s*:") if key != "kmt_file"
                    else re.compile(rf"{key}\s*=\s*(.+)$")
                    for key in PARAM_KEYS}
        result = {key: "" for key in PARAM_KEYS}
        try:
            with open(self.P_ice_diag, "r", encoding="utf-8", errors="replace") as f:
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
            self.sim_config = None
            return False
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
        with open(self.P_sim_cfg, "w") as f:
            json.dump(result, f, indent=2)
        self.sim_config = result
        return True

    ##########################################################################################################
    #############################     XARRAY/NUMPY DATASET/ARRAY EXTRACTION/WORK     ##########################
    ##########################################################################################################
    def align_time_coordinate_of_three_arrays(self, ds1, ds2, ds3, time_coord="time"):
        for da in [ds1, ds2, ds3]:
            da[time_coord] = pd.to_datetime(da[time_coord].values).normalize()
        t_common = np.intersect1d(np.intersect1d(ds1[time_coord].values, ds2[time_coord].values), ds3[time_coord].values)
        return ds1.sel(time=t_common), ds2.sel(time=t_common), ds3.sel(time=t_common)

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

    def _has(self, ds, var):
        """
        Return whether the wrapped input contains a variable/key/attribute named ``var``.

        This helper provides a uniform “does it exist?” check across common container types
        encountered in the sea-ice workflow (e.g., ``xarray.Dataset`` outputs, dict-like
        objects, or lightweight objects with attributes).

        Lookup precedence
        -----------------
        1. If ``ds`` is an ``xarray.Dataset``: check ``var`` in ``ds.data_vars``.
        2. Otherwise attempt membership: ``var in ds`` (for dict-like / list-like containers).
        3. If membership is not supported (raises ``TypeError``): fall back to ``hasattr(I_data, var)``.

        Parameters
        ----------
        ds  : xr.Dataset
              dataset to test var
        var : str
              variable name to test for.

        Returns
        -------
        bool
            ``True`` if ``var`` is present under the rules above, otherwise ``False``.

        Notes
        -----
        - For ``xarray.Dataset``, this checks *data variables only* (not coordinates).
          If you want coordinates too, consider checking ``(var in I_data.variables)``.
        - The function relies on an outer-scope variable ``ds`` (closure). If you
          refactor to store data on the instance, replace ``ds`` with ``self.ds``.
        """
        if isinstance(ds, xr.Dataset):
            return var in ds.data_vars
        try:
            return var in ds
        except TypeError:
            return hasattr(ds, var)

    def _get(self, ds, var):
        """
        Retrieve a variable/key/attribute named ``var`` from the wrapped input.

        This helper provides a uniform accessor across common container types
        (e.g., ``xarray.Dataset`` outputs, dict-like objects, or objects with attributes).

        Lookup precedence
        -----------------
        1. If ``ds`` is an ``xarray.Dataset``: return ``ds[var]`` (typically an
           ``xarray.DataArray`` or variable-like object).
        2. Otherwise attempt key access: ``ds[var]``.
        3. If key access fails for any reason: fall back to attribute access
           ``getattr(ds, var)``.

        Parameters
        ----------
        var : str
            Variable name to retrieve.

        Returns
        -------
        Any
            The retrieved object. For ``xarray.Dataset`` this is usually an
            ``xarray.DataArray``.

        Raises
        ------
        KeyError
            If ``var`` is not a key in a mapping-like ``ds`` and attribute fallback
            does not exist.
        AttributeError
            If key access fails and ``ds`` does not have attribute ``var``.
        Exception
            Any exception raised by ``ds[var]`` may be swallowed and replaced by
            the attribute fallback attempt.

        Notes
        -----
        - Because the fallback catches *all* exceptions from ``ds[var]``, genuine
          indexing errors (not just missing keys) will be masked. If you only want to
          fall back on missing keys, narrow the exception to ``(KeyError, TypeError)``.
        - The function relies on an outer-scope variable ``ds`` (closure). If you
          refactor to store data on the instance, replace ``ds`` with ``self.ds``.
        """
        if isinstance(ds, xr.Dataset):
            return ds[var]
        try:
            return ds[var]
        except Exception:
            return getattr(ds, var)

    ##########################################################################################################
    #############################                  MASKING                          ##########################
    ##########################################################################################################
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

    ##########################################################################################################
    #################################### PATH/DIRECTORY DEFINITIONS #########################################
    ##########################################################################################################
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

    @staticmethod
    def _class_method_key(class_method: str) -> str:
        lut = {"raw"         : "raw",
               "daily"       : "raw",
               "binary-days" : "bin",
               "binary_days" : "bin",
               "bin"         : "bin",
               "rolling-mean": "roll",
               "rolling_mean": "roll",
               "roll"        : "roll"}
        if class_method not in lut:
            raise ValueError(f"Unsupported class_method: {class_method}")
        return lut[class_method]

    def _resolve_product_store(self, ice_type,
                               class_method = "raw",
                               product      = "data",
                               D_zarr       = None,
                               BorC2T_type  = None,
                               ispd_thresh  = None,
                               bin_win_days = None,
                               bin_min_days = None,
                               mean_period  = None):
        self.define_toolbox_paths(D_zarr       = D_zarr,
                                  ice_type     = ice_type,
                                  BorC2T_type  = BorC2T_type,
                                  ispd_thresh  = ispd_thresh,
                                  bin_win_days = bin_win_days,
                                  bin_min_days = bin_min_days,
                                  mean_period  = mean_period)
        method_key = self._class_method_key(class_method)
        if ice_type in {"SI", "MIZ"}:
            if method_key != "raw":
                raise ValueError(f"{ice_type} only supports class_method='raw'")
            return self.P_zarrs[ice_type][product]
        return self.P_zarrs[ice_type][product][method_key]

    def define_toolbox_paths(self,
                             sim_name     = None,
                             D_sim        = None,
                             D_zarr       = None,
                             D_graph      = None,
                             ice_type     = None,
                             BorC2T_type  = None,
                             ispd_thresh  = None,
                             bin_win_days = None,
                             bin_min_days = None,
                             mean_period  = None):
        """
        Define canonical paths for ice-history, classified products, and metrics.

        Layout:
          <D_zarr>/
            iceh_daily.zarr/<YYYY-MM>/
            iceh_monthly.zarr/<YYYY-MM>/
            <HEM>/
              SI/
                raw.zarr
                mets.zarr
              MIZ/
                raw.zarr
                mets.zarr
              ispd_thresh_<thr>/
                FI/
                  <BorC2T>/
                    raw.zarr
                    mets.zarr
                    bin-win-XX_bin-min-YY/
                      raw.zarr
                      mets.zarr
                    roll-days-ZZ/
                      raw.zarr
                      mets.zarr
                PI/
                  <BorC2T>/
                    raw.zarr
                    mets.zarr
                    bin-win-XX_bin-min-YY/
                      raw.zarr
                      mets.zarr
                    roll-days-ZZ/
                      raw.zarr
                      mets.zarr
        """
        sim_name     = sim_name          if sim_name     is not None else self.sim_name
        self.D_sim   = Path(D_sim        if D_sim        is not None else Path(self.D_dict['AFIM_out'], sim_name))
        self.D_zarr  = Path(D_zarr       if D_zarr       is not None else Path(self.D_sim , 'zarr'))
        self.D_graph = Path(D_graph      if D_graph      is not None else Path(self.config['D_dict']['graph'], 'AFIM'))
        ice_type     = ice_type          if ice_type     is not None else self.ice_type
        BorC2T_type  = BorC2T_type       if BorC2T_type  is not None else self.BorC2T_type
        ispd_thresh  = float(ispd_thresh if ispd_thresh  is not None else self.ispd_thresh)
        bin_win_days = int(bin_win_days  if bin_win_days is not None else self.bin_win_days)
        bin_min_days = int(bin_min_days  if bin_min_days is not None else self.bin_min_days)
        mean_period  = int(mean_period   if mean_period  is not None else self.mean_period)
        self._check_ice_type(ice_type)
        self._check_BorC2T_type(BorC2T_type)
        self.P_ice_diag = self.D_sim / "ice_diag.d"
        self.P_sim_cfg  = self.D_sim / f"ice_in_AFIM_subset_{sim_name}.json"
        if isinstance(BorC2T_type, str):
            self.reG_type = BorC2T_type
        else:
            self.reG_type = "".join(BorC2T_type)
        ispd_thr_str = f"{ispd_thresh:.1e}".replace("e-0", "e-").replace("e+0", "e+")
        self.D_hem   = self.D_zarr / self.hemisphere_dict["abbreviation"]
        self.D_ispd  = self.D_hem / f"ispd_thresh_{ispd_thr_str}"
        self.D_iceh  = {"zarr": {"daily":   self.D_zarr / "iceh_daily.zarr",
                                 "monthly": self.D_zarr / "iceh_monthly.zarr"},
                        "nc"  : {"daily":   self.D_sim  / "history" / "daily",
                                 "monthly": self.D_sim  / "history" / "monthly"}}
        self.D_class = {"SI" : self.D_hem / "SI",
                        "MIZ": self.D_hem / "MIZ",
                        "FI" : {"raw":  self.D_ispd / "FI" / self.reG_type / "raw",
                                "bin":  self.D_ispd / "FI" / self.reG_type / f"bin-win-{bin_win_days:02d}_bin-min-{bin_min_days:02d}",
                                "roll": self.D_ispd / "FI" / self.reG_type / f"roll-days-{mean_period:02d}"},
                        "PI" : {"raw":  self.D_ispd / "PI" / self.reG_type / "raw",
                                "bin":  self.D_ispd / "PI" / self.reG_type / f"bin-win-{bin_win_days:02d}_bin-min-{bin_min_days:02d}",
                                "roll": self.D_ispd / "PI" / self.reG_type / f"roll-days-{mean_period:02d}"}}
        self.P_zarrs = {"SI" : {"data": self.D_class["SI"] / "raw.zarr",
                                "mets": self.D_class["SI"] / "mets.zarr"},
                        "MIZ": {"data": self.D_class["MIZ"] / "raw.zarr",
                                "mets": self.D_class["MIZ"] / "mets.zarr"},
                        "FI" : {"data": {"raw":  self.D_class["FI"]["raw"]  / "data.zarr",
                                         "bin":  self.D_class["FI"]["bin"]  / "data.zarr",
                                         "roll": self.D_class["FI"]["roll"] / "data.zarr"},
                                "mets": {"raw":  self.D_class["FI"]["raw"]  / "mets.zarr",
                                         "bin":  self.D_class["FI"]["bin"]  / "mets.zarr",
                                         "roll": self.D_class["FI"]["roll"] / "mets.zarr"}},
                        "PI" : {"data": {"raw":  self.D_class["PI"]["raw"]  / "data.zarr",
                                         "bin":  self.D_class["PI"]["bin"]  / "data.zarr",
                                         "roll": self.D_class["PI"]["roll"] / "data.zarr"},
                                "mets": {"raw":  self.D_class["PI"]["raw"]  / "mets.zarr",
                                         "bin":  self.D_class["PI"]["bin"]  / "mets.zarr",
                                         "roll": self.D_class["PI"]["roll"] / "mets.zarr"}}}
        self.methods_init_executed['define_toolbox_paths'] = True
        return True

    def define_ice_mask_name(self, ice_type=None):
        ice_type = ice_type or self.ice_type
        self._check_ice_type(ice_type)
        self.mask_name = f"{ice_type}_mask"

    def _sanitize_for_zarr_write(self, obj: xr.Dataset | xr.DataArray,
                                 drop_coords: list | None = None,
                                 cast_float32: bool = False) -> xr.Dataset | xr.DataArray:
        """
    Remove non-essential auxiliary coordinates before writing classified products.

    This is primarily to prevent year-to-year inconsistencies in saved Zarr groups,
    e.g. one year carrying NLON/NLAT while another does not, which later breaks
    xarray.concat in load_classified_ice().

    Parameters
    ----------
    obj : xr.Dataset or xr.DataArray
        Classified product to be written.
    drop_coords : list[str], optional
        Extra coords to remove. If None, uses a conservative default plus any
        entries in self.CICE_dict["drop_coords"].
    cast_float32 : bool, default False
        If True, cast float64 data variables to float32 before write.

    Returns
    -------
    xr.Dataset or xr.DataArray
        Sanitised object ready for Zarr write.
        """
        default_drop = ["NLON", "NLAT",
                        "ULON", "ULAT",
                        "TLON", "TLAT",
                        "elon", "elat",
                        "nlon", "nlat",
                        "ulon", "ulat",
                        "tlon", "tlat"]
        cfg_drop = []
        if hasattr(self, "CICE_dict") and isinstance(self.CICE_dict, dict):
            cfg_drop = list(self.CICE_dict.get("drop_coords", []))
        all_drop = drop_coords or list(dict.fromkeys(default_drop + cfg_drop))
        if isinstance(obj, xr.DataArray):
            da = obj
            present = [c for c in all_drop if c in da.coords]
            if present:
                self.logger.debug(f"write-sanitize DataArray: dropping coords {present}")
                da = da.drop_vars(present, errors="ignore")
            if cast_float32 and np.issubdtype(da.dtype, np.floating) and da.dtype != np.float32:
                da = da.astype(np.float32)
            return da
        ds = obj
        present = [c for c in all_drop if c in ds.coords]
        if present:
            self.logger.debug(f"write-sanitize Dataset: dropping coords {present}")
            ds = ds.drop_vars(present, errors="ignore")
        # Optional float downcast for storage consistency
        ds = ds.drop_vars("time_bounds", errors="ignore")
        if "time" in ds.coords:
            ds["time"].attrs.pop("bounds", None)
        if cast_float32:
            ds = ds.copy()
            for v in ds.data_vars:
                if np.issubdtype(ds[v].dtype, np.floating) and ds[v].dtype != np.float32:
                    ds[v] = ds[v].astype(np.float32)
        # if cast_float32:
        #     cast_map = {}
        #     for v in ds.data_vars:
        #         if np.issubdtype(ds[v].dtype, np.floating) and ds[v].dtype != np.float32:
        #             cast_map[v] = np.float32
        #     if cast_map:
        #         ds = ds.astype(cast_map)
        return ds

    def _clean_var_encoding(self, da: xr.DataArray,
                            keep_chunks     : bool = None,
                            keep_compressor : bool = None) -> xr.DataArray:
        enc0 = dict(getattr(da, "encoding", {}) or {})
        enc1 = {}
        if keep_chunks and "chunks" in enc0:
            enc1["chunks"] = enc0["chunks"]
        if keep_compressor and "compressor" in enc0:
            enc1["compressor"] = enc0["compressor"]
        # keep fill value if present
        if "_FillValue" in enc0:
            enc1["_FillValue"] = enc0["_FillValue"]
        da.encoding = enc1
        return da

    def _strip_unsafe_zarr_encoding(self, obj: xr.Dataset | xr.DataArray,
                                    keep_chunks    : bool = True,
                                    keep_compressor: bool = True) -> xr.Dataset | xr.DataArray:
        """
        Remove problematic / inherited variable encodings before Zarr write.

        This is especially useful when objects inherit encodings from upstream reads
        and are later written with xarray.to_zarr(). We keep only a small safe subset
        by default.
        """
        if isinstance(obj, xr.DataArray):
            return self._clean_var_encoding(obj, keep_chunks=keep_chunks, keep_compressor=keep_compressor)
        ds = obj.copy()
        for v in ds.variables:
            ds[v].encoding = {}  # clear first
        for v in ds.data_vars:
            ds[v] = self._clean_var_encoding(ds[v], keep_chunks=keep_chunks, keep_compressor=keep_compressor)
        for c in ds.coords:
            ds[c].encoding = {}
        return ds

    def _drop_duplicate_coords(self, ds: xr.Dataset,
                               dim: str = "ni") -> xr.Dataset:
        """
        Drop duplicate coordinate values along `dim` (keeping the first occurrence).
        Useful when concatenating yearly groups that may contain duplicated x-indices.
        """
        if dim in ds.coords:
            _, index = np.unique(ds[dim], return_index=True)
            ds = ds.isel({dim: sorted(index)})
        return ds

    def _normalise_concat_coords(self, ds: xr.Dataset,
                                 dim: str = "ni",
                                 drop_coords: list | None = None) -> xr.Dataset:
        """
        Normalise non-essential spatial coords before concatenation of yearly groups.

        This protects against cases where some yearly Zarr groups carry auxiliary
        coordinates (for example NLON/NLAT) and others do not, which causes
        xarray.concat(..., coords="minimal") to raise:

            ValueError: coordinate 'NLON' not present in all datasets

        Parameters
        ----------
        ds : xr.Dataset
            Dataset for a single year/group.
        dim : str, default "ni"
            Optional dimension along which duplicate coordinate values should be
            removed, keeping first occurrence.
        drop_coords : list[str], optional
            Explicit list of coord names to drop before concat. If None, uses a
            conservative default set of known non-essential geolocation coords.

        Returns
        -------
        xr.Dataset
            Dataset with problematic auxiliary coords removed and duplicate concat
            coordinates cleaned.
        """
        # 1) drop duplicate concat coord values, if present
        if dim in ds.coords:
            _, index = np.unique(ds[dim], return_index=True)
            ds = ds.isel({dim: sorted(index)})
        # 2) drop auxiliary geolocation coords that are not needed for concat
        default_drop = ["NLON", "NLAT",
                        "ULON", "ULAT",
                        "TLON", "TLAT",
                        "elon", "elat",
                        "nlon", "nlat",
                        "ulat", "ulon",
                        "tlat", "tlon"]
        cfg_drop = []
        if hasattr(self, "CICE_dict") and isinstance(self.CICE_dict, dict):
            cfg_drop = list(self.CICE_dict.get("drop_coords", []))
        drop_now = drop_coords or list(dict.fromkeys(default_drop + cfg_drop))
        present  = [c for c in drop_now if c in ds.coords]
        if present:
            self.logger.debug(f"Dropping auxiliary coords before concat: {present}")
            ds = ds.drop_vars(present, errors="ignore")
        return ds

    def _write_grouped_zarr(self, ds, store, group, overwrite_group=False, consolidated=False):
        store = Path(store)
        P_grp = store / group
        if P_grp.exists():
            if overwrite_group:
                self.logger.info(f"group {P_grp} already exists but overwrite_group = True ... OVER-WRITTING GROUP")
                shutil.rmtree(P_grp)
            else:
                self.logger.info(f"Group already exists, skipping: {P_grp}")
                return False
        ds   = self._sanitize_for_zarr_write(ds, cast_float32=True)
        ds   = self._strip_unsafe_zarr_encoding(ds)
        mode = "a" if store.exists() else "w"
        self.logger.info(f"writing to group {P_grp}")
        self.logger.info(f"writing to group with mode '{mode}' and consolidated {consolidated}")
        ds.to_zarr(store, group=group, mode=mode, consolidated=consolidated, zarr_format=2)
        return True

    ##########################################################################################################
    ####################################      NORMALISATIONS         #########################################
    ##########################################################################################################
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

    ##########################################################################################################
    ####################################   DATE/TIME MANIPULATIONS   #########################################
    ##########################################################################################################
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

    ##########################################################################################################
    ################################### BASIC STATISTICS/CLIMATOLOGY #########################################
    ##########################################################################################################
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
