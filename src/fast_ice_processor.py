import xarray as xr
import pandas as pd
import numpy as np
from pathlib import Path
from datetime import datetime, timedelta
import xesmf as xe
from scipy.spatial import cKDTree
import json
import os
import contextlib
import sys
import time
import logging
sys.path.insert(0, '/home/581/da1339/AFIM/src/AFIM/src')
from grounded_iceberg_processor import GroundedIcebergProcessor 

@contextlib.contextmanager
def suppress_stderr():
    with open(os.devnull, 'w') as devnull:
        old_stderr = sys.stderr
        sys.stderr = devnull
        try:
            yield
        finally:
            sys.stderr = old_stderr

class FastIceProcessor:
    """
    Class to compute Fast Ice metrics from sea ice model output.

    This processor applies a rolling window to compute spatial and temporal
    characteristics of fast ice using sea ice concentration and velocity thresholds.

    Fast ice is defined as sea ice with concentration above a threshold and
    velocity below a threshold (i.e., stationary sea ice).

    Parameters
    ----------
    sim_name  : str
                Name of the simulation (used to resolve directory paths and configuration).
    json_path : str, optional
                Path to the JSON configuration file. If None, the default config is used.
    roll_win  : int, optional
                Rolling window size in days (default is 15 or value specified in config).
    P_log     : str or Path, optional
                Path to output log file. Defaults to `logs/FastIceProcessor_{sim_name}.log`.

    Attributes
    ----------
    config              : dict
                          Full configuration dictionary loaded from JSON.
    roll_win            : int
                          Length of the rolling window in days.
    var_list            : list
                          List of variable names to extract from CICE history files.
    chunk_dict          : dict
                          Dictionary of chunking strategy for xarray open_mfdataset.
    FI_thresh           : float
                          Speed threshold (m/s) below which ice is considered fast.
    SIC_thresh          : float
                          Sea ice concentration threshold for defining fast ice.
    SIC_scale           : float
                          Scaling factor to convert area units (typically to 1e6 km²).
    cm2m_fact           : float
                          Conversion factor for cm to m.
    sim_dir             : Path
                          Path to simulation output directory for daily NetCDF files.
    CICE_dict           : dict
                          Dictionary of CICE-specific variable/dimension names.
    regrid_weights_path : str or Path
                          File path for xESMF weights from u-grid to t-grid regridding.
    gi_processor        : GroundedIcebergProcessor
                          Grid and landmask handler for simulation-specific preprocessing.
    use_gi              : bool
                          If True, grounded iceberg landmask is used during processing.

    Examples
    --------
    Basic usage (single timestep):

    >>> from datetime import datetime
    >>> proc = FastIceProcessor("Rothrock")
    >>> ds = proc.process_window(datetime(1994, 9, 8), save_zarr=False)

    Batch processing over time:

    The FastIceProcessor is typically run in a loop using a driver script
    such as [`compute_fast_ice.py`](https://github.com/dpath2o/AFIM/blob/main/src/python/compute_fast_ice.py):

    python:
    from datetime import datetime, timedelta
    from fast_ice_processor import FastIceProcessor

    dt_start = datetime(1993, 1, 1)
    dt_end   = datetime(1999, 12, 31)
    proc = FastIceProcessor("Rothrock", roll_win=15)
    current_date = dt_start + timedelta(days=15)

    while current_date <= dt_end - timedelta(days=15):
        proc.process_window(current_date)
        current_date += timedelta(days=15)

    For more, see the full project repository:
    🔗 https://github.com/dpath2o/AFIM
    """

    def __init__(self, sim_name, json_path=None, roll_win=None, P_log=None):
        """
        Initialize the FastIceProcessor object.

        Loads configuration, sets simulation paths, initializes logging,
        and sets up grounded iceberg processor.

        Parameters
        ----------
        sim_name  : str
                    Name of the simulation to process.
        json_path : str or Path, optional
                    Path to the JSON configuration file.
        roll_win  : int, optional
                    Rolling window size in days.
        P_log     : str or Path, optional
                    Path to log file for output logging.
        """
        self.sim_name = sim_name
        # Load config
        if json_path is None:
            json_path = "/home/581/da1339/AFIM/src/AFIM/JSONs/afim_cice_analysis.json"
        with open(json_path, 'r') as f:
            self.config = json.load(f)
        # Rolling window
        self.roll_win = roll_win or self.config.get('roll_win', 15)
        # Logging
        if P_log is None:
            P_log = Path(self.config['D_dict']['logs'], f"FastIceProcessor_{sim_name}.log")
        self.setup_logging(logfile=P_log)
        # Core simulation config
        self.sim_config      = self.config['sim_dict'][sim_name]
        self.json_hemisphere = self.config.get('hemisphere', 'south')
        self.var_list        = self.config['CICE_dict']['FI_vars']
        self.chunk_dict      = self.config['CICE_dict']['FI_chunks']
        self.FI_thresh       = self.config.get("FI_thresh", 5e-4)
        self.FIC_scale       = self.config.get("FIC_scale", 1e9)
        self.SIC_thresh      = self.config.get("SIC_thresh", 0.15)
        self.SIC_scale       = self.config.get("SIC_scale", 1e12)
        self.cm2m_fact       = self.config.get("cm2m_fact", 0.01)
        # File system paths
        self.sim_dir              = Path(self.config['D_dict']['AFIM_out'], sim_name, "history", "daily")
        self.CICE_dict            = self.config['CICE_dict']
        self.regrid_weights_path  = self.CICE_dict['P_reG_u2t_weights']
        # Instantiate grounded iceberg/grid processor
        self.gi_processor = GroundedIcebergProcessor(self.config, sim_name)
        self.use_gi       = self.gi_processor.use_gi
        # For convenience, reference GI total area directly
        self.GI_total_area = self.gi_processor.total_area if self.use_gi else 0

    def setup_logging(self, logfile=None):
        """
        Setup logging to console and file.

        Parameters
        ----------
        logfile : str or Path, optional
                  Log file path. If None, defaults to stdout only.
        """
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

    def define_hemisphere(self, hemisphere):
        if hemisphere is None:
            hemisphere = self.json_hemisphere
        if hemisphere.lower() in ['north', 'northern', 'nh', 'n', 'no']:
            self.hemisphere_geographic_extent = [0, 360, 0, 90]
            self.hemisphere_map_extent        = [-180,180,55,90]
            self.hemisphere_projection        = "S0.0/90.0/50/15C"
            self.hemisphere_map_text_location = [-120,56]
            self.hemisphere_abbreviation      = 'NH'
            self.hemisphere_nj_slice          = slice(540,1080)
        elif hemisphere.lower() in ['south', 'southern', 'sh', 's', 'so']:
            self.hemisphere_geographic_extent = [0, 360, -90, 0]
            self.hemisphere_map_extent        = [-180,180,-90,-55]
            self.hemisphere_projection        = "S0.0/-90.0/50/15C"
            self.hemisphere_map_text_location = [0,-90]
            self.hemisphere_abbreviation      = 'SH'
            self.hemisphere_nj_slice          = slice(0,540)
        else:
            raise ValueError(f"Invalid hemisphere '{hemisphere}'. Valid options are: "
                             "['north', 'south', 'northern', 'southern', 'sh', 'nh', 'SH', 'NH']")

    def slice_hemisphere(self, ds):
        """
        Define which hemisphere slice to use.

        Parameters
        ----------
        hemisphere : str
                     Hemisphere keyword ('north' or 'south') to configure slicing.
        """
        return ds.isel(nj=self.hemisphere_nj_slice)

    def load_data_window(self):
        """
        Slice dataset for the currently defined hemisphere.

        Parameters
        ----------
        ds : xarray.Dataset
             Input dataset.

        Returns
        -------
        xarray.Dataset
            Dataset sliced along nj dimension.
        """
        dates = pd.date_range(self.dtC - pd.Timedelta(days=self.roll_win // 2),
                              self.dtC + pd.Timedelta(days=self.roll_win // 2))
        self.dt0_str = dates[0].strftime('%Y-%m-%d')
        self.dtN_str = dates[-1].strftime('%Y-%m-%d')
        files = [self.sim_dir / f"iceh.{d.strftime('%Y-%m-%d')}.nc" for d in dates if (self.sim_dir / f"iceh.{d.strftime('%Y-%m-%d')}.nc").exists()]
        if not files:
            raise FileNotFoundError(f"No files found for window around {self.dtC}")
        preprocess = lambda ds: ds[self.var_list]
        ds = xr.open_mfdataset(files, combine='by_coords', parallel=True,
                               preprocess=preprocess, chunks=self.chunk_dict)
        return ds

    def regrid_to_tgrid(self, ds):
        """
        Regrid vector components (uvel, vvel) from u-grid to t-grid.

        Parameters
        ----------
        ds : xarray.Dataset
             Dataset containing 'uvel' and 'vvel'.

        Returns
        -------
        tuple of xarray.DataArray
            Regridded u and v velocities.
        """
        self.logger.info("Regridding uvel and vvel to T-grid...")
        with suppress_stderr():
            regridder = xe.Regridder(self.gi_processor.G_u, self.gi_processor.G_t, method="bilinear", extrap_method="inverse_dist", periodic=True, weights=self.regrid_weights_path)
        U_tgrid = regridder(ds["uvel"])
        V_tgrid = regridder(ds["vvel"])
        U_tgrid = self.slice_hemisphere(U_tgrid)
        V_tgrid = self.slice_hemisphere(V_tgrid)
        return U_tgrid, V_tgrid

    def compute_rolling_averages(self, ds):
        """
        Compute rolling averages over time for each fast ice-related variable.

        Parameters
        ----------
        ds : xarray.Dataset
             Dataset of raw model fields.

        Returns
        -------
        dict of xarray.DataArray
            Dictionary with variable names as keys and rolling-averaged DataArrays as values.
        """
        U_tgrid, V_tgrid = self.regrid_to_tgrid(ds)
        ds = self.slice_hemisphere(ds)
        self.logger.info("Computing rolling averages...")
        roll = lambda v: ds[v].rolling(time=self.roll_win, center=True, min_periods=1).mean()
        return {
            'aice': roll('aice'),
            'hi': roll('hi'),
            'strength': roll('strength'),
            'shear': roll('shear'),
            'divu': roll('divu'),
            'iage': roll('iage'),
            'daidtd': roll('daidtd'),
            'daidtt': roll('daidtt'),
            'dvidtd': roll('dvidtd'),
            'dvidtt': roll('dvidtt'),
            'strint': np.sqrt(ds['strintx']**2 + ds['strinty']**2).rolling(time=self.roll_win, center=True, min_periods=1).mean(),
            'speed': np.sqrt(U_tgrid**2 + V_tgrid**2).rolling(time=self.roll_win, center=True, min_periods=1).mean()
        }

    def apply_FI_mask(self, roll_dict):
        """
        Apply fast ice masking based on thresholds for concentration and velocity.

        Parameters
        ----------
        roll_dict : dict
                    Dictionary of rolling-averaged variables.

        Returns
        -------
        tuple
            Masked variables and binary fast ice mask.
        """    
        mask = (roll_dict['aice'] > self.SIC_thresh) & (roll_dict['speed'] <= self.FI_thresh)
        return {k: v.where(mask) for k, v in roll_dict.items()}, mask

    def compute_fast_ice_outputs(self, roll_mask_vars):
        """
        Compute fast ice metrics: 1D, 2D, and 3D outputs for spatial and temporal aggregation.

        This method is the backbone of the class and is intended to be called by `process_window()`.

        Parameters
        ----------
        roll_mask_vars : dict of xarray.DataArray
                         Masked variables used for final computations.

        Returns
        -------
        xarray.Dataset
            Final dataset with all fast ice metrics.
        """
        roll_mask_vars = {k: v.compute() for k, v in roll_mask_vars.items()}
        grid_cell_area = self.gi_processor.G_t['area'].isel(nj=self.hemisphere_nj_slice)
        spatial_dims   = self.CICE_dict['spatial_dims']
        time_dim       = self.CICE_dict['time_dim']
        time_coords    = roll_mask_vars['aice'].time.values
        lon_coords     = self.gi_processor.G_t['lon'].isel(nj=self.hemisphere_nj_slice).values
        lat_coords     = self.gi_processor.G_t['lat'].isel(nj=self.hemisphere_nj_slice).values
        lon_coord_name = self.CICE_dict['FI_lon_coord']
        lat_coord_name = self.CICE_dict['FI_lat_coord']
        FI_time_dim    = self.CICE_dict['FI_time_dim']
        cm2m           = self.cm2m_fact
        roll_win       = self.roll_win
        self.logger.info("compute fast ice 1D variables:")
        t0       = time.time()
        fia_roll = ((roll_mask_vars['aice']        * grid_cell_area).sum(dim=spatial_dims, skipna=True) / self.FIC_scale)
        fiv_roll = ((roll_mask_vars['hi']          * grid_cell_area).sum(dim=spatial_dims, skipna=True) / self.FIC_scale)
        fisth_td = ((roll_mask_vars['strength']    * grid_cell_area).sum(dim=spatial_dims, skipna=True))
        fish_td  = ((roll_mask_vars['shear']       * grid_cell_area).sum(dim=spatial_dims, skipna=True))
        fidiv_td = ((roll_mask_vars['divu']        * grid_cell_area).sum(dim=spatial_dims, skipna=True))
        fiag_td  = ((grid_cell_area / roll_mask_vars['iage']).sum(dim=spatial_dims, skipna=True))
        fiad_td  = ((roll_mask_vars['daidtd']      * grid_cell_area).sum(dim=spatial_dims, skipna=True))
        fiat_td  = ((roll_mask_vars['daidtt']      * grid_cell_area).sum(dim=spatial_dims, skipna=True))
        fivd_td  = ((roll_mask_vars['dvidtd']*cm2m * grid_cell_area).sum(dim=spatial_dims, skipna=True))
        fivt_td  = ((roll_mask_vars['dvidtt']*cm2m * grid_cell_area).sum(dim=spatial_dims, skipna=True))
        fistr_td = ((roll_mask_vars['strint']      * grid_cell_area).sum(dim=spatial_dims, skipna=True))
        if self.use_gi:
            fia = fia_roll.values + self.GI_total_area
            fiv = fiv_roll.values + self.GI_total_area
        else:
            fia = fia_roll.values
            fiv = fiv_roll.values
        FIA       = (time_dim,
                    fia,
                    {'units'      : '1000-km^2',
                     'long_name'  : 'fast ice area',
                     'description': 'sea ice area summed over spatial extent and masked with fast ice criteria'})
        FIV       = (time_dim,
                    fiv,
                    {'units'      : '1000-km^3',
                     'long_name'  : 'fast ice volume',
                     'description': 'sea ice volume summed over spatial extent and masked with fast ice criteria'})
        FI_STRONG = (time_dim,
                    fisth_td.data,
                    {'units'      : 'N',
                     'long_name'  : 'fast ice strength',
                     'description': 'sea ice strength times grid cell area summed over spatial extent and masked with fast ice criteria'})
        FI_SHEAR  = (time_dim,
                    fish_td.data,
                    {'units'      : 'm^2/day',
                     'long_name'  : 'fast ice shear rate',
                     'description': 'sea ice shear times grid cell area summed over spatial extent and masked with fast ice criteria'})
        FI_DIV    = (time_dim,
                     fidiv_td.data,
                     {'units'      : 'm^2/day',
                      'long_name'  : 'fast ice divergence rate',
                      'description': 'sea ice divergence times grid cell area summed over spatial extent and masked with fast ice criteria'})
        FI_AGING  = (time_dim,
                    fiag_td.data,
                    {'units'      : 'm^2/year',
                     'long_name'  : 'fast ice ageing rate',
                     'description': 'sea ice area divided by ice age summed over spatial extent and masked with fast ice criteria'})
        FI_AGROM  = (time_dim,
                    fiad_td.data,
                    {'units'      : 'm^2/day',
                     'long_name'  : 'fast ice mechanical area growth rate',
                     'description': 'sea ice area tendency dynamics (mechanical) times grid cell area summed over spatial extent and masked with fast ice criteria'})
        FI_AGROT  = (time_dim,
                    fiat_td.data,
                    {'units'      : 'm^2/day',
                     'long_name'  : 'fast ice thermodynamic area growth rate',
                     'description': 'sea ice area tendency thermodynamics times grid cell area summed over spatial extent and masked with fast ice criteria'})
        FI_VGROM  = (time_dim,
                    fivd_td.data,
                    {'units'      : 'm^3/day',
                     'long_name'  : 'fast ice mechanical volume growth rate',
                     'description': 'sea ice volume tendency dynamics (mechanical) times grid cell area summed over spatial extent and masked with fast ice criteria'})
        FI_VGROT  = (time_dim,
                    fivt_td.data,
                    {'units'      : 'm^3/day',
                     'long_name'  : 'fast ice thermodynamic volume growth rate',
                     'description': 'sea ice volume tendency thermodynamics times grid cell area summed over spatial extent and masked with fast ice criteria'})
        FI_STRESS = (time_dim,
                     fistr_td.data,
                     {'units'      : 'N',
                      'long_name'  : 'fast ice internal stress',
                      'description': 'sea ice internal stress times grid cell area summed over spatial extent and masked with fast ice criteria'})
        self.logger.info(f"\ttime taken: {time.time()-t0} seconds")
        #################################
        #######   3D VARIABLES   ########
        ################################
        self.logger.info(f"coarsen fast ice rolling averages to {roll_win}-days--i.e. 3D variables")
        t0    = time.time()
        fic   = roll_mask_vars['aice'].coarsen(time=roll_win, boundary="trim").mean()
        fihi  = roll_mask_vars['hi'].coarsen(time=roll_win, boundary="trim").mean()
        fisth = roll_mask_vars['strength'].coarsen(time=roll_win, boundary="trim").mean()
        fish  = roll_mask_vars['shear'].coarsen(time=roll_win, boundary="trim").mean()
        fidiv = roll_mask_vars['divu'].coarsen(time=roll_win, boundary="trim").mean()
        fiag  = roll_mask_vars['iage'].coarsen(time=roll_win, boundary="trim").mean()
        fiad  = roll_mask_vars['daidtd'].coarsen(time=roll_win, boundary="trim").mean()
        fiat  = roll_mask_vars['daidtt'].coarsen(time=roll_win, boundary="trim").mean()
        fivd  = roll_mask_vars['dvidtd'].coarsen(time=roll_win, boundary="trim").mean()
        fivt  = roll_mask_vars['dvidtt'].coarsen(time=roll_win, boundary="trim").mean()
        fistr = roll_mask_vars['strint'].coarsen(time=roll_win, boundary="trim").mean()
        FI_time_coords = fic.time.values
        FIC   = (self.CICE_dict['FI_three_dims'],
                 fic.values,
                 {'units'      : '%/grid-cell',
                  'long_name'  : 'fast ice concentration',
                  'description': 'sea ice concentration per grid cell masked with fast ice criteria and coarsened to fast ice mask temporal window'})
        FIHI  = (self.CICE_dict['FI_three_dims'],
                 fihi.values,
                 {'units'      : 'm',
                  'long_name'  : 'fast ice thickness',
                  'description': 'sea ice thickness per grid cell masked with fast ice criteria and coarsened to fast ice mask temporal window'})
        FISTH = (self.CICE_dict['FI_three_dims'],
                 fisth.values,
                 {'units'      : 'N/m',
                  'long_name'  : 'fast ice strength',
                  'description': 'sea ice strength per grid cell masked with fast ice criteria and coarsened to fast ice mask temporal window'})
        FISH  = (self.CICE_dict['FI_three_dims'],
                 fish.values,
                 {'units'      : '%/day',
                  'long_name'  : 'fast ice shear',
                  'description': 'sea ice shear per grid cell masked with fast ice criteria and coarsened to fast ice mask temporal window'})
        FIDIV = (self.CICE_dict['FI_three_dims'],
                 fidiv.values,
                 {'units'      : '%/day',
                  'long_name'  : 'fast ice divergence',
                  'description': 'sea ice divergence per grid cell masked with fast ice criteria and coarsened to fast ice mask temporal window'})
        FIAG  = (self.CICE_dict['FI_three_dims'],
                 fiag.values,
                 {'units'      : 'years',
                  'long_name'  : 'fast ice age',
                  'description': 'sea ice age per grid cell masked with fast ice criteria and coarsened to fast ice mask temporal window'})
        FIAD  = (self.CICE_dict['FI_three_dims'],
                 fiad.values,
                 {'units'      : '%/day',
                  'long_name'  : 'fast ice area tendency dynamics (mechanical)',
                  'description': 'sea ice area tendency dynamics (mechanical) per grid cell masked with fast ice criteria and coarsened to fast ice mask temporal window'})
        FIAT  = (self.CICE_dict['FI_three_dims'],
                 fiat.values,
                 {'units'      : '%/day',
                  'long_name'  : 'fast ice area tendency thermodynamics',
                  'description': 'sea ice area tendency thermodynamics per grid cell masked with fast ice criteria and coarsened to fast ice mask temporal window'})
        FIVD  = (self.CICE_dict['FI_three_dims'],
                 fivd.values,
                 {'units'      : 'cm/day',
                  'long_name'  : 'fast ice volume tendency dynamics (mechanical)',
                  'description': 'sea ice volume tendency dynamics (mechanical) per grid cell masked with fast ice criteria and coarsened to fast ice mask temporal window'})
        FIVT  = (self.CICE_dict['FI_three_dims'],
                 fivt.values,
                 {'units'      : 'cm/day',
                  'long_name'  : 'fast ice volume tendency thermodynamics',
                  'description': 'sea ice volume tendency thermodynamics per grid cell masked with fast ice criteria and coarsened to fast ice mask temporal window'})
        FISTR = (self.CICE_dict['FI_three_dims'],
                 fistr.values,
                 {'units'      : 'N/m^2',
                  'long_name'  : 'fast ice internal stress',
                  'description': 'sea ice internal stress per grid cell masked with fast ice criteria and coarsened to fast ice mask temporal window'})
        self.logger.info(f"\ttime taken: {time.time()-t0} seconds")
        #################################
        #######   2D VARIABLES   ########
        ################################
        self.logger.info("compute temporal sums to give spatial distributions over time--i.e. 2D variables")
        self.logger.info(f"\ttemporal mean over the period {self.dt0_str} to {self.dtN_str}")
        t0       = time.time()
        fip      = fic.sum(dim=time_dim)   / len(fic.time)                                 #units: % grid cell covered in sea ice
        fihi_sd  = fihi.sum(dim=time_dim)  / (len(fihi.time)  * roll_mask_vars['hi'].max().values)      #units: m
        fisth_sd = fisth.sum(dim=time_dim) / (len(fisth.time) * roll_mask_vars['strength'].max().values)#units: N/m
        fish_sd  = fish.sum(dim=time_dim)  / (len(fish.time)  * roll_mask_vars['shear'].max().values)   #units: 1/day
        fidiv_sd = fidiv.sum(dim=time_dim) / (len(fidiv.time) * roll_mask_vars['divu'].max().values)    #units: 1/day
        fiag_sd  = fiag.sum(dim=time_dim)  / (len(fiag.time)  * roll_mask_vars['iage'].max().values)    #units: years
        fiad_sd  = fiad.sum(dim=time_dim)  / (len(fiad.time)  * roll_mask_vars['daidtd'].max().values)  #units: 1/day
        fiat_sd  = fiat.sum(dim=time_dim)  / (len(fiat.time)  * roll_mask_vars['daidtt'].max().values)  #units: 1/day
        fivd_sd  = fivd.sum(dim=time_dim)  / (len(fivd.time)  * roll_mask_vars['dvidtd'].max().values)  #units: cm/day
        fivt_sd  = fivt.sum(dim=time_dim)  / (len(fivt.time)  * roll_mask_vars['dvidtt'].max().values)  #units: cm/day
        fistr_sd = fistr.sum(dim=time_dim) / (len(fistr.time) * roll_mask_vars['strint'].max().values)  #units: N/m^2
        FIP      = (spatial_dims,
                    fip.values,
                    {'units'       : '% (()/)',
                     'long_name'   : 'fast ice persistence',
                     'description' : 'sum of sea ice concentration per grid cell over time masked with fast ice criteria and divide by total temporal length',
                     'start_time'  : self.dt0_str,
                     'stop_time'   : self.dtN_str})
        FIHI_SD  = (spatial_dims,
                    fihi_sd.values,
                    {'units'       : '% (m)',
                     'long_name'   : 'fast ice thickness spatial distribution over time',
                     'description' : 'sum of sea ice thickness per grid cell over time multiplied by maximum thickness and masked with fast ice criteria',
                     'start_time'  : self.dt0_str,
                     'stop_time'   : self.dtN_str})
        FISTH_SD = (spatial_dims,
                    fisth_sd.values,
                    {'units'       : '% (N/m)',
                     'long_name'   : 'fast ice strength spatial distribution over time',
                     'description' : 'sum of sea ice strength per grid cell over time divided by max strength and masked with fast ice criteria',
                     'start_time'  : self.dt0_str,
                     'stop_time'   : self.dtN_str})
        FISH_SD  = (spatial_dims,
                    fish_sd.values,
                    {'units'       : '% (1/day)',
                     'long_name'   : 'fast ice shear spatial distribution over time',
                     'description' : 'sum of sea ice shear per grid cell over time divided by max shear and masked with fast ice criteria',
                     'start_time'  : self.dt0_str,
                     'stop_time'   : self.dtN_str})
        FIDIV_SD = (spatial_dims,
                    fidiv_sd.values,
                    {'units'       : '% (1/day)',
                     'long_name'   : 'fast ice divergence spatial distribution over time',
                     'description' : 'sum of sea ice divergence per grid cell over time divided by max divergence and masked with fast ice criteria',
                     'start_time'  : self.dt0_str,
                     'stop_time'   : self.dtN_str})
        FIAG_SD  = (spatial_dims,
                    fiag_sd.values,
                    {'units'       : '% (years)',
                     'long_name'   : 'fast ice age spatial distribution over time',
                     'description' : 'sum of sea ice age per grid cell over time divided by max age and masked with fast ice criteria',
                     'start_time'  : self.dt0_str,
                     'stop_time'   : self.dtN_str})
        FIAD_SD  = (spatial_dims,
                    fiad_sd.values,
                    {'units'       : '% (1/day)',
                     'long_name'   : 'fast ice area-mechanical-tendency spatial distribution over time',
                     'description' : 'sum of ice area-mechanical-tendency per grid cell over time divided by max decrement and masked with fast ice criteria',
                     'start_time'  : self.dt0_str,
                     'stop_time'   : self.dtN_str})
        FIAT_SD  = (spatial_dims,
                    fiat_sd.values,
                    {'units'       : '% (1/day)',
                     'long_name'   : 'fast ice area-thermodynamic-tendency spatial distribution over time',
                     'description' : 'sum of ice area-thermodynamic-tendency per grid cell over time divided by max increment and masked with fast ice criteria',
                     'start_time'  : self.dt0_str,
                     'stop_time'   : self.dtN_str})
        FIVD_SD  = (spatial_dims,
                    fivd_sd.values*cm2m,
                    {'units'       : '% (m/day)',
                     'long_name'   : 'fast ice volume-mechanical-tendency spatial distribution over time',
                     'description' : 'sum of ice volume-mechanical-tendency decrement per grid cell over time divided by max decrement and masked with fast ice criteria',
                     'start_time'  : self.dt0_str,
                     'stop_time'   : self.dtN_str})
        FIVT_SD  = (spatial_dims,
                    fivt_sd.values*cm2m,
                    {'units'       : '% (m/day)',
                     'long_name'   : 'fast ice volume-thermodynamic-tendency spatial distribution over time',
                     'description' : 'sum of ice volume-thermodynamic-tendency per grid cell over time divided by max increment and masked with fast ice criteria',
                     'start_time'  : self.dt0_str,
                     'stop_time'   : self.dtN_str})
        FISTR_SD = (spatial_dims,
                    fistr_sd.values,
                    {'units'       : '% (N/m^2)',
                     'long_name'   : 'fast ice internal stress spatial distribution over time',
                     'description' : 'sum of ice internal stress per grid cell over time divided by max stress and masked with fast ice criteria',
                     'start_time'  : self.dt0_str,
                     'stop_time'   : self.dtN_str})
        self.logger.info(f"\ttime taken: {time.time()-t0} seconds")
        #######################################
        #######      OUTPUT DATASET     #######
        #######################################
        self.logger.info("create output dataset:")
        t0 = time.time()
        FI = xr.Dataset(
            {# 1D variables
                'FIA'       : FIA,
                'FIV'       : FIV,
                'FI_STRONG' : FI_STRONG,
                'FI_SHEAR'  : FI_SHEAR,
                'FI_DIV'    : FI_DIV,
                'FI_AGING'  : FI_AGING,
                'FI_AGROM'  : FI_AGROM,
                'FI_AGROT'  : FI_AGROT,
                'FI_VGROM'  : FI_VGROM,
                'FI_VGROT'  : FI_VGROT,
                'FI_STRESS' : FI_STRESS,
            # 2D variables
                'FIP'      : FIP,
                'FIHI_SD'  : FIHI_SD,
                'FISTH_SD' : FISTH_SD,
                'FISH_SD'  : FISH_SD,
                'FIDIV_SD' : FIDIV_SD,
                'FIAG_SD'  : FIAG_SD,
                'FIAD_SD'  : FIAD_SD,
                'FIAT_SD'  : FIAT_SD,
                'FIVD_SD'  : FIVD_SD,
                'FIVT_SD'  : FIVT_SD,
                'FISTR_SD' : FISTR_SD,
            # 3D variables
                'FIC'   : FIC,
                'FIHI'  : FIHI,
                'FISTH' : FISTH,
                'FISH'  : FISH,
                'FIDIV' : FIDIV,
                'FIAG'  : FIAG,
                'FIAD'  : FIAD,
                'FIAT'  : FIAT,
                'FIVD'  : FIVD,
                'FIVT'  : FIVT,
                'FISTR' : FISTR
            },
            coords = {time_dim        : (time_dim     , time_coords),
                      FI_time_dim     : (FI_time_dim  , FI_time_coords),
                      lon_coord_name  : (spatial_dims , lon_coords),
                      lat_coord_name  : (spatial_dims , lat_coords) }
        )
        FI.attrs = {"title"              : "Landfast sea ice analysed from numerical sea ice model simulations",
                    "summary"            : "This dataset includes landfast sea ice variables and derived metrics "\
                                           "using a rolling window method then masking variables for threshold "\
                                           "values of sea ice concentration and sea ice speed.",
                    "source"             : "CICE v6.4.1 model output",
                    "creator_name"       : "Daniel Patrick Atwater",
                    "creator_email"      : "daniel.atwater@utas.edu.au",
                    "institution"        : "Institute of Marine and Antarctic Studies--University of Tasmania",
                    "history"            : f"Created on {datetime.now().isoformat()}",
                    "references"         : "",
                    "conventions"        : "CF-1.8",
                    "fast_ice_criteria"  : f"aice > {self.SIC_thresh} and speed <= {self.FI_thresh} m/s",
                    "grounded_iceberg_db": self.gi_processor.GI_dataset_path,
                    "landmask_file"      : self.gi_processor.KMT_path,
                    "total_area_GI"      : self.GI_total_area,
                    "time_coverage_start": self.dt0_str,
                    "time_coverage_end"  : self.dtN_str,
                    "roll_window_days"   : roll_win,
                    "geospatial_lat_min" : float(np.min(lat_coords)),
                    "geospatial_lat_max" : float(np.max(lat_coords)),
                    "geospatial_lon_min" : float(np.min(lon_coords)),
                    "geospatial_lon_max" : float(np.max(lon_coords))}
        self.logger.info(f"\ttime taken: {time.time()-t0} seconds")
        return FI

    def process_window(self, dtC, hemisphere='south', ow_zarrs=False, save_zarrs=True):
        """
        Main entry point to process fast ice metrics for a single central date.

        Parameters
        ----------
        dtC        : datetime
                     Central date for the rolling window.
        hemisphere : str
                     Hemisphere to process ('north' or 'south').
        ow_zarrs   : bool
                     Overwrite existing Zarr file if True.
        save_zarr  : bool
                     Save result to disk as a Zarr dataset if True.

        Returns
        -------
        xarray.Dataset
            Fast ice dataset containing all computed metrics.
        """
        self.dtC = dtC
        self.logger.info(f"\n\nProcessing window centered on {self.dtC} for {hemisphere}ern hemisphere")
        ds = self.load_data_window()
        self.define_hemisphere(hemisphere)
        if self.use_gi:
            self.gi_processor.load_AFIM_GI()
        else:
            self.gi_processor.load_grid_and_landmask()
        roll_avg_dict = self.compute_rolling_averages(ds)
        FI_masked_vars, _ = self.apply_FI_mask(roll_avg_dict)
        FI = self.compute_fast_ice_outputs(FI_masked_vars)
        if save_zarrs:
            D_FI_zarr = Path(self.config['D_dict']['AFIM_out'], self.sim_name, "FI")
            F_FI_zarr = f"fast_ice_{self.dtC.strftime('%Y-%m-%d')}.zarr"
            if not os.path.exists(D_FI_zarr):
                os.makedirs(D_FI_zarr)
            P_FI_zarr = Path(D_FI_zarr,F_FI_zarr)
            if not os.path.exists(P_FI_zarr) or (os.path.exists(P_FI_zarr) and ow_zarrs):
                self.logger.info(f"*** writing FI to disk: {P_FI_zarr}")
                t0 = time.time()
                FI.to_zarr(P_FI_zarr, mode='w')
                self.logger.info(f"\ttime taken: {time.time()-t0} seconds")
            else:
                self.logger.info("FI already exists or over-writting zarrs disabled ***")
                self.logger.info(f"\tskipping: {P_FI_zarr}")
        return FI
