import xarray as xr
import pandas as pd
import numpy as np
from pathlib import Path
from datetime import datetime
import xesmf as xe
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

class PackIceProcessor:
    """
    Class to compute Pack Ice metrics from sea ice model output.

    This processor uses a rolling time window to compute spatial and temporal
    statistics of pack ice from a CICE-based sea ice model simulation.
    Pack ice is defined as sea ice with concentration above a threshold and
    speed above a threshold (i.e., not fast ice).

    Parameters
    ----------
    sim_name  : str
                Name of the simulation (used to look up paths and configs).
    json_path : str, optional
                Path to the JSON config file. Defaults to the default project config.
    roll_win  : int, optional
                Rolling window (in days) for averaging. Defaults to value in config (usually 15).
    P_log     : str or Path, optional
                File path to write logs. If None, logs are saved to `logs/PackIceProcessor_{sim_name}.log`.

    Attributes
    ----------
    config               : dict
                           JSON configuration dictionary.
    roll_win             : int
                           Rolling window size in days.
    var_list             : list
                           List of CICE variable names to load.
    chunk_dict           : dict
                           Dictionary defining chunking strategy for Dask/xarray.
    FI_thresh            : float
                           Speed threshold (m/s) to distinguish fast ice (not used in pack ice masking).
    SIC_thresh           : float
                           Concentration threshold for sea ice to be considered pack ice.
    SIC_scale            : float
                           Scaling factor for area-based computations (e.g., to convert m^2 to 10^6 km^2).
    cm2m_fact            : float
                           Conversion factor from cm to m (typically 0.01).
    sim_dir              : Path
                           Directory containing simulation daily NetCDF files.
    CICE_dict            : dict
                           Dictionary of CICE-specific configuration values.
    regrid_weights_path  : str or Path
                           Path to precomputed regridding weights file for u- to t-grid interpolation.
    gi_processor         : GroundedIcebergProcessor
                           Instance for grid and mask operations. GI is not used for pack ice.
    use_gi               : bool
                           Always False for pack ice processor.

    Example
    -------
    Simple single-date use case:

    >>> from datetime import datetime
    >>> PI_proc = PackIceProcessor("Rothrock")
    >>> PI = PI_proc.process_window(datetime(1994, 9, 8), save_zarr=False)

    Looping over a full simulation (recommended usage):

    See the external loop driver script `compute_pack_ice.py`, which uses this
    class to process a time series:

    python:
    from datetime import datetime, timedelta
    from pack_ice_processor import PackIceProcessor

    dt_start = datetime(1993, 1, 1)
    dt_end   = datetime(1999, 12, 31)
    processor = PackIceProcessor("Rothrock", roll_win=15)
    current_date = dt_start + timedelta(days=15)

    while current_date <= dt_end - timedelta(days=15):
        PI = processor.process_window(current_date)
        current_date += timedelta(days=15)
    
    See full codebase at: https://github.com/dpath2o/AFIM
    """

    def __init__(self, sim_name, json_path=None, roll_win=None, P_log=None):
        self.sim_name = sim_name
        if json_path is None:
            json_path = "/home/581/da1339/AFIM/src/AFIM/JSONs/afim_cice_analysis.json"
        with open(json_path, 'r') as f:
            self.config = json.load(f)

        self.roll_win = roll_win or self.config.get('roll_win', 15)

        if P_log is None:
            P_log = Path(self.config['D_dict']['logs'], f"PackIceProcessor_{sim_name}.log")
        self.setup_logging(logfile=P_log)

        self.sim_config      = self.config['sim_dict'][sim_name]
        self.json_hemisphere = self.config.get('hemisphere', 'south')
        self.var_list        = self.config['CICE_dict']['FI_vars']
        self.chunk_dict      = self.config['CICE_dict']['FI_chunks']
        self.FI_thresh       = self.config.get("FI_thresh", 5e-4)
        self.SIC_thresh      = self.config.get("SIC_thresh", 0.15)
        self.SIC_scale       = self.config.get("SIC_thresh", 1e12)
        self.cm2m_fact       = self.config.get("cm2m_fact", 0.01)

        self.sim_dir              = Path(self.config['D_dict']['AFIM_out'], sim_name, "history", "daily")
        self.CICE_dict            = self.config['CICE_dict']
        self.regrid_weights_path  = self.CICE_dict['P_reG_u2t_weights']

        self.gi_processor = GroundedIcebergProcessor(self.config, sim_name)
        self.use_gi       = False  # GI excluded for pack ice

    def setup_logging(self, logfile=None):
        """
        Configure logging for the PackIceProcessor.

        Parameters
        ----------
        logfile : str or Path, optional
                  File path for log output. If None, logs are not written to file.
        """
        self.logger = logging.getLogger(self.sim_name + '_PI')
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
        """
        Set the hemisphere and associated nj slice.

        Parameters
        ----------
        hemisphere : str
                     Hemisphere to define ('north' or 'south'). Determines nj slice.
        """
        if hemisphere is None:
            hemisphere = self.json_hemisphere
        if hemisphere.lower() in ['north', 'northern', 'nh', 'n', 'no']:
            self.hemisphere_nj_slice = slice(540, 1080)
        elif hemisphere.lower() in ['south', 'southern', 'sh', 's', 'so']:
            self.hemisphere_nj_slice = slice(0, 540)
        else:
            raise ValueError(f"Invalid hemisphere: {hemisphere}")

    def slice_hemisphere(self, ds):
        """
        Apply hemisphere-specific slicing to dataset.

        Parameters
        ----------
        ds : xarray.Dataset
             Dataset to slice.

        Returns
        -------
        ds_sliced : xarray.Dataset
                    Dataset restricted to nj range of chosen hemisphere.
        """
        return ds.isel(nj=self.hemisphere_nj_slice)

    def load_data_window(self):
        """
        Load NetCDF files spanning the rolling window centered on the current date.

        Returns
        -------
        ds : xarray.Dataset
             Dataset containing selected variables over rolling window.

        Raises
        ------
        FileNotFoundError
            If no NetCDF files found for the window around dtC.
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
        self.logger.info(f"Loaded dataset with shape: {ds.sizes}")
        if 'aice' not in ds or ds['aice'].sizes['time'] < self.roll_win:
            self.logger.warning("aice is missing or too short after loading files")
        return ds

    def regrid_to_tgrid(self, ds):
        """
        Regrid velocity components from u-grid to t-grid using xESMF.

        Parameters
        ----------
        ds : xarray.Dataset
             Dataset containing 'uvel' and 'vvel'.

        Returns
        -------
        U, V : xarray.DataArray
               Regridded u- and v-velocity fields on the t-grid.
        """
        self.logger.info("Regridding uvel and vvel to T-grid...")
        with suppress_stderr():
            regridder = xe.Regridder(self.gi_processor.G_u, self.gi_processor.G_t,
                                     method="bilinear", extrap_method="inverse_dist", periodic=True,
                                     weights=self.regrid_weights_path)
        U = regridder(ds['uvel'])
        V = regridder(ds['vvel'])
        return self.slice_hemisphere(U), self.slice_hemisphere(V)

    def compute_rolling_averages(self, ds):
        """
        Compute rolling averages over specified time window for relevant variables.

        Parameters
        ----------
        ds : xarray.Dataset
             Dataset to compute rolling means on.

        Returns
        -------
        roll_dict : dict of xarray.DataArray
                    Dictionary of rolled variables including derived metrics.
        """
        U, V = self.regrid_to_tgrid(ds)
        ds = self.slice_hemisphere(ds)
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
            'frazil': roll('frazil'),
            'strint': np.sqrt(ds['strintx']**2 + ds['strinty']**2).rolling(time=self.roll_win, center=True, min_periods=1).mean(),
            'speed': np.sqrt(U**2 + V**2).rolling(time=self.roll_win, center=True, min_periods=1).mean()
        }

    def apply_PI_mask(self, roll_dict):
        """
        Apply pack ice mask based on concentration and speed thresholds.

        Parameters
        ----------
        roll_dict : dict
                    Dictionary of rolling-averaged variables.

        Returns
        -------
        masked_dict : dict
                      Dictionary with pack ice mask applied.
        pi_mask     : xarray.DataArray
                      Boolean mask array where pack ice criteria are met.
        """
        fi_mask = (roll_dict['aice'] > self.SIC_thresh) & (roll_dict['speed'] <= self.FI_thresh)
        pi_mask = (roll_dict['aice'] > self.SIC_thresh) & (~fi_mask)
        return {k: v.where(pi_mask) for k, v in roll_dict.items()}, pi_mask

    def compute_pack_ice_outputs(self, roll_mask_vars):
        """
        Compute summary metrics and spatial/temporal diagnostics from masked variables.

        This is the backbone method of the class and is intended to be called
        from the `process_window()` method. It performs the heavy-lifting for
        summarising pack ice across multiple time scales and dimensions.

        Parameters
        ----------
        roll_mask_vars : dict
                         Dictionary of masked and rolled variables.

        Returns
        -------
        PI : xarray.Dataset
             Output dataset containing 1D, 2D, and 3D pack ice metrics.

        Notes
        -----
        This method computes metrics like total area, volume, mechanical/thermo growth,
        stress, shear, divergence, aging, speed, and frazil formation. Includes
        both time-series, coarsened 3D fields, and spatial frequency summaries.
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
        self.logger.info("compute pack ice 1D variables:")
        t0       = time.time()
        pia_roll = ((roll_mask_vars['aice']        * grid_cell_area).sum(dim=spatial_dims, skipna=True) / self.SIC_scale)
        piv_roll = ((roll_mask_vars['hi']          * grid_cell_area).sum(dim=spatial_dims, skipna=True) / self.SIC_scale)
        pisth_td = ((roll_mask_vars['strength']    * grid_cell_area).sum(dim=spatial_dims, skipna=True))
        pish_td  = ((roll_mask_vars['shear']       * grid_cell_area).sum(dim=spatial_dims, skipna=True))
        pidiv_td = ((roll_mask_vars['divu']        * grid_cell_area).sum(dim=spatial_dims, skipna=True))
        piag_td  = ((grid_cell_area / roll_mask_vars['iage']).sum(dim=spatial_dims, skipna=True))
        piad_td  = ((roll_mask_vars['daidtd']      * grid_cell_area).sum(dim=spatial_dims, skipna=True))
        piat_td  = ((roll_mask_vars['daidtt']      * grid_cell_area).sum(dim=spatial_dims, skipna=True))
        pivd_td  = ((roll_mask_vars['dvidtd']*cm2m * grid_cell_area).sum(dim=spatial_dims, skipna=True))
        pivt_td  = ((roll_mask_vars['dvidtt']*cm2m * grid_cell_area).sum(dim=spatial_dims, skipna=True))
        pistr_td = ((roll_mask_vars['strint']      * grid_cell_area).sum(dim=spatial_dims, skipna=True))
        pispd_td = ((roll_mask_vars['speed']       * grid_cell_area).sum(dim=spatial_dims, skipna=True))
        pifzl_td = ((roll_mask_vars['frazil']*cm2m * grid_cell_area).sum(dim=spatial_dims, skipna=True))
        PIA       = (time_dim,
                    pia_roll.data,
                    {'units'      : '1e6-km^2',
                     'long_name'  : 'pack ice area',
                     'description': 'pack ice area summed over spatial extent and masked for pack ice criteria'})
        PIV       = (time_dim,
                    piv_roll.data,
                    {'units'      : '1e6-km^3',
                     'long_name'  : 'pack ice volume',
                     'description': 'pack ice volume summed over spatial extent and masked for pack ice criteria'})
        PI_STRONG = (time_dim,
                    pisth_td.data,
                    {'units'      : 'N',
                     'long_name'  : 'pack ice strength',
                     'description': 'sea ice strength times grid cell area summed over spatial extent and masked with pack ice criteria'})
        PI_SHEAR  = (time_dim,
                    pish_td.data,
                    {'units'      : 'm^2/day',
                     'long_name'  : 'pack ice shear rate',
                     'description': 'sea ice shear times grid cell area summed over spatial extent and masked with pack ice criteria'})
        PI_DIV    = (time_dim,
                     pidiv_td.data,
                     {'units'      : 'm^2/day',
                      'long_name'  : 'pack ice divergence rate',
                      'description': 'sea ice divergence times grid cell area summed over spatial extent and masked with pack ice criteria'})
        PI_AGING  = (time_dim,
                    piag_td.data,
                    {'units'      : 'm^2/year',
                     'long_name'  : 'pack ice ageing rate',
                     'description': 'sea ice area divided by ice age summed over spatial extent and masked with pack ice criteria'})
        PI_AGROM  = (time_dim,
                    piad_td.data,
                    {'units'      : 'm^2/day',
                     'long_name'  : 'pack ice mechanical area growth rate',
                     'description': 'sea ice area tendency dynamics (mechanical) times grid cell area summed over spatial extent and masked with pack ice criteria'})
        PI_AGROT  = (time_dim,
                    piat_td.data,
                    {'units'      : 'm^2/day',
                     'long_name'  : 'pack ice thermodynamic area growth rate',
                     'description': 'sea ice area tendency thermodynamics times grid cell area summed over spatial extent and masked with pack ice criteria'})
        PI_VGROM  = (time_dim,
                    pivd_td.data,
                    {'units'      : 'm^3/day',
                     'long_name'  : 'pack ice mechanical volume growth rate',
                     'description': 'sea ice volume tendency dynamics (mechanical) times grid cell area summed over spatial extent and masked with pack ice criteria'})
        PI_VGROT  = (time_dim,
                    pivt_td.data,
                    {'units'      : 'm^3/day',
                     'long_name'  : 'pack ice thermodynamic volume growth rate',
                     'description': 'sea ice volume tendency thermodynamics times grid cell area summed over spatial extent and masked with pack ice criteria'})
        PI_STRESS = (time_dim,
                     pistr_td.data,
                     {'units'      : 'N',
                      'long_name'  : 'pack ice internal stress',
                      'description': 'sea ice internal stress times grid cell area summed over spatial extent and masked with pack ice criteria'})
        PI_SPEED  = (time_dim,
                    pispd_td.data,
                    {'units'      : 'm^3/second',
                     'long_name'  : 'pack ice speed',
                     'description': 'sea ice speed times grid cell area summed over spatial extent and masked with pack ice criteria'})
        PI_FRAZIL = (time_dim,
                    pifzl_td.data,
                    {'units'      : 'm^3/day',
                     'long_name'  : 'pack ice frazil volumetric growth',
                     'description': 'sea ice frazil growth rate times grid cell area summed over spatial extent and masked with pack ice criteria'})
        self.logger.info(f"\ttime taken: {time.time()-t0} seconds")
        self.logger.info(f"\ttime taken: {time.time()-t0} seconds")
        #################################
        #######   3D VARIABLES   ########
        ################################
        self.logger.info(f"coarsen pack ice rolling averages to {roll_win}-days--i.e. 3D variables")
        t0    = time.time()
        pic   = roll_mask_vars['aice'].coarsen(time=roll_win, boundary="trim").mean()
        pihi  = roll_mask_vars['hi'].coarsen(time=roll_win, boundary="trim").mean()
        pisth = roll_mask_vars['strength'].coarsen(time=roll_win, boundary="trim").mean()
        pish  = roll_mask_vars['shear'].coarsen(time=roll_win, boundary="trim").mean()
        pidiv = roll_mask_vars['divu'].coarsen(time=roll_win, boundary="trim").mean()
        piag  = roll_mask_vars['iage'].coarsen(time=roll_win, boundary="trim").mean()
        piad  = roll_mask_vars['daidtd'].coarsen(time=roll_win, boundary="trim").mean()
        piat  = roll_mask_vars['daidtt'].coarsen(time=roll_win, boundary="trim").mean()
        pivd  = roll_mask_vars['dvidtd'].coarsen(time=roll_win, boundary="trim").mean()
        pivt  = roll_mask_vars['dvidtt'].coarsen(time=roll_win, boundary="trim").mean()
        pistr = roll_mask_vars['strint'].coarsen(time=roll_win, boundary="trim").mean()
        pispd = roll_mask_vars['speed'].coarsen(time=roll_win, boundary="trim").mean()
        pifzl = roll_mask_vars['frazil'].coarsen(time=roll_win, boundary="trim").mean()*cm2m
        PI_time_coords = pic.time.values
        PIC   = (self.CICE_dict['FI_three_dims'],
                 pic.values,
                 {'units'      : '%/grid-cell',
                  'long_name'  : 'pack ice concentration',
                  'description': 'sea ice concentration per grid cell masked with pack ice criteria and coarsened to pack ice mask temporal window'})
        PIHI  = (self.CICE_dict['FI_three_dims'],
                 pihi.values,
                 {'units'      : 'm',
                  'long_name'  : 'pack ice thickness',
                  'description': 'sea ice thickness per grid cell masked with pack ice criteria and coarsened to pack ice mask temporal window'})
        PISTH = (self.CICE_dict['FI_three_dims'],
                 pisth.values,
                 {'units'      : 'N/m',
                  'long_name'  : 'pack ice strength',
                  'description': 'sea ice strength per grid cell masked with pack ice criteria and coarsened to pack ice mask temporal window'})
        PISH  = (self.CICE_dict['FI_three_dims'],
                 pish.values,
                 {'units'      : '%/day',
                  'long_name'  : 'pack ice shear',
                  'description': 'sea ice shear per grid cell masked with pack ice criteria and coarsened to pack ice mask temporal window'})
        PIDIV = (self.CICE_dict['FI_three_dims'],
                 pidiv.values,
                 {'units'      : '%/day',
                  'long_name'  : 'pack ice divergence',
                  'description': 'sea ice divergence per grid cell masked with pack ice criteria and coarsened to pack ice mask temporal window'})
        PIAG  = (self.CICE_dict['FI_three_dims'],
                 piag.values,
                 {'units'      : 'years',
                  'long_name'  : 'pack ice age',
                  'description': 'sea ice age per grid cell masked with pack ice criteria and coarsened to pack ice mask temporal window'})
        PIAD  = (self.CICE_dict['FI_three_dims'],
                 piad.values,
                 {'units'      : '%/day',
                  'long_name'  : 'pack ice area tendency dynamics (mechanical)',
                  'description': 'sea ice area tendency dynamics (mechanical) per grid cell masked with pack ice criteria and coarsened to pack ice mask temporal window'})
        PIAT  = (self.CICE_dict['FI_three_dims'],
                 piat.values,
                 {'units'      : '%/day',
                  'long_name'  : 'pack ice area tendency thermodynamics',
                  'description': 'sea ice area tendency thermodynamics per grid cell masked with pack ice criteria and coarsened to pack ice mask temporal window'})
        PIVD  = (self.CICE_dict['FI_three_dims'],
                 pivd.values,
                 {'units'      : 'cm/day',
                  'long_name'  : 'pack ice volume tendency dynamics (mechanical)',
                  'description': 'sea ice volume tendency dynamics (mechanical) per grid cell masked with pack ice criteria and coarsened to pack ice mask temporal window'})
        PIVT  = (self.CICE_dict['FI_three_dims'],
                 pivt.values,
                 {'units'      : 'cm/day',
                  'long_name'  : 'pack ice volume tendency thermodynamics',
                  'description': 'sea ice volume tendency thermodynamics per grid cell masked with pack ice criteria and coarsened to pack ice mask temporal window'})
        PISTR = (self.CICE_dict['FI_three_dims'],
                 pistr.values,
                 {'units'      : 'N/m^2',
                  'long_name'  : 'pack ice internal stress',
                  'description': 'sea ice internal stress per grid cell masked with pack ice criteria and coarsened to pack ice mask temporal window'})
        PISPD = (self.CICE_dict['FI_three_dims'],
                 pispd.values,
                 {'units'      : 'm/s',
                  'long_name'  : 'pack ice speed',
                  'description': 'sea ice speed per grid cell masked with pack ice criteria and coarsened to pack ice mask temporal window'})
        PIFZL = (self.CICE_dict['FI_three_dims'],
                 pifzl.values,
                 {'units'      : 'm/day',
                  'long_name'  : 'pack ice frazil growth rate',
                  'description': 'sea ice frazil growth rate per grid cell masked with pack ice criteria and coarsened to pack ice mask temporal window'})
        self.logger.info(f"\ttime taken: {time.time()-t0} seconds")
        #################################
        #######   2D VARIABLES   ########
        ################################
        self.logger.info("compute temporal sums to give spatial distributions over time--i.e. 2D variables")
        self.logger.info(f"\ttemporal mean over the period {self.dt0_str} to {self.dtN_str}")
        t0       = time.time()
        pic_sd   = pic.sum(dim=time_dim)   / len(pic.time)                                              #units: % grid cell covered in sea ice
        pihi_sd  = pihi.sum(dim=time_dim)  / (len(pihi.time)  * roll_mask_vars['hi'].max().values)      #units: m
        pisth_sd = pisth.sum(dim=time_dim) / (len(pisth.time) * roll_mask_vars['strength'].max().values)#units: N/m
        pish_sd  = pish.sum(dim=time_dim)  / (len(pish.time)  * roll_mask_vars['shear'].max().values)   #units: 1/day
        pidiv_sd = pidiv.sum(dim=time_dim) / (len(pidiv.time) * roll_mask_vars['divu'].max().values)    #units: 1/day
        piag_sd  = piag.sum(dim=time_dim)  / (len(piag.time)  * roll_mask_vars['iage'].max().values)    #units: years
        piad_sd  = piad.sum(dim=time_dim)  / (len(piad.time)  * roll_mask_vars['daidtd'].max().values)  #units: 1/day
        piat_sd  = piat.sum(dim=time_dim)  / (len(piat.time)  * roll_mask_vars['daidtt'].max().values)  #units: 1/day
        pivd_sd  = pivd.sum(dim=time_dim)  / (len(pivd.time)  * roll_mask_vars['dvidtd'].max().values)  #units: cm/day
        pivt_sd  = pivt.sum(dim=time_dim)  / (len(pivt.time)  * roll_mask_vars['dvidtt'].max().values)  #units: cm/day
        pistr_sd = pistr.sum(dim=time_dim) / (len(pistr.time) * roll_mask_vars['strint'].max().values)  #units: N/m^2
        pispd_sd = pispd.sum(dim=time_dim) / (len(pispd.time) * roll_mask_vars['speed'].max().values)   #units: m/s
        pifzl_sd = pifzl.sum(dim=time_dim) / (len(pifzl.time) * roll_mask_vars['frazil'].max().values)   #units: cm/day
        PIC_SD   = (spatial_dims,
                    pic_sd.values,
                    {'units'       : '% (1/day)',
                     'long_name'   : 'pack ice concentration frequency',
                     'description' : 'sum of sea ice concentration per grid cell over time masked with pack ice criteria and divide by total temporal length',
                     'start_time'  : self.dt0_str,
                     'stop_time'   : self.dtN_str})
        PIHI_SD  = (spatial_dims,
                    pihi_sd.values,
                    {'units'       : '% (1/day)',
                     'long_name'   : 'pack ice thickness frequency',
                     'description' : 'sum of sea ice thickness per grid cell over time multiplied by maximum thickness and masked with pack ice criteria',
                     'start_time'  : self.dt0_str,
                     'stop_time'   : self.dtN_str})
        PISTH_SD = (spatial_dims,
                    pisth_sd.values,
                    {'units'       : '% (1/day)',
                     'long_name'   : 'pack ice strength frequency',
                     'description' : 'sum of sea ice strength per grid cell over time divided by max strength and masked with pack ice criteria',
                     'start_time'  : self.dt0_str,
                     'stop_time'   : self.dtN_str})
        PISH_SD  = (spatial_dims,
                    pish_sd.values,
                    {'units'       : '% (1/day)',
                     'long_name'   : 'pack ice shear frequency',
                     'description' : 'sum of sea ice shear per grid cell over time divided by max shear and masked with pack ice criteria',
                     'start_time'  : self.dt0_str,
                     'stop_time'   : self.dtN_str})
        PIDIV_SD = (spatial_dims,
                    pidiv_sd.values,
                    {'units'       : '% (1/day)',
                     'long_name'   : 'pack ice divergence frequency',
                     'description' : 'sum of sea ice divergence per grid cell over time divided by max divergence and masked with pack ice criteria',
                     'start_time'  : self.dt0_str,
                     'stop_time'   : self.dtN_str})
        PIAG_SD  = (spatial_dims,
                    piag_sd.values,
                    {'units'       : '% (1/day)',
                     'long_name'   : 'pack ice age frequency',
                     'description' : 'sum of sea ice age per grid cell over time divided by max age and masked with pack ice criteria',
                     'start_time'  : self.dt0_str,
                     'stop_time'   : self.dtN_str})
        PIAD_SD  = (spatial_dims,
                    piad_sd.values,
                    {'units'       : '% (1/day)',
                     'long_name'   : 'pack ice area-mechanical-tendency frequency',
                     'description' : 'sum of sea ice area-mechanical-tendency per grid cell over time divided by max decrement and masked with pack ice criteria',
                     'start_time'  : self.dt0_str,
                     'stop_time'   : self.dtN_str})
        PIAT_SD  = (spatial_dims,
                    piat_sd.values,
                    {'units'       : '% (1/day)',
                     'long_name'   : 'pack ice area-thermodynamic-tendency frequency',
                     'description' : 'sum of sea ice area-thermodynamic-tendency per grid cell over time divided by max increment and masked with pack ice criteria',
                     'start_time'  : self.dt0_str,
                     'stop_time'   : self.dtN_str})
        PIVD_SD  = (spatial_dims,
                    pivd_sd.values*cm2m,
                    {'units'       : '% (1/day)',
                     'long_name'   : 'pack ice volume-mechanical-tendency frequency',
                     'description' : 'sum of sea ice volume-mechanical-tendency decrement per grid cell over time divided by max decrement and masked with pack ice criteria',
                     'start_time'  : self.dt0_str,
                     'stop_time'   : self.dtN_str})
        PIVT_SD  = (spatial_dims,
                    pivt_sd.values*cm2m,
                    {'units'       : '% (1/day)',
                     'long_name'   : 'pack ice volume-thermodynamic-tendency frequency',
                     'description' : 'sum of sea ice volume-thermodynamic-tendency per grid cell over time divided by max increment and masked with pack ice criteria',
                     'start_time'  : self.dt0_str,
                     'stop_time'   : self.dtN_str})
        PISTR_SD = (spatial_dims,
                    pistr_sd.values,
                    {'units'       : '% (1/day)',
                     'long_name'   : 'pack ice internal stress frequency',
                     'description' : 'sum of sea ice internal stress per grid cell over time divided by max stress and masked with pack ice criteria',
                     'start_time'  : self.dt0_str,
                     'stop_time'   : self.dtN_str})
        PISPD_SD = (spatial_dims,
                    pispd_sd.values,
                    {'units'       : '% (1/day)',
                     'long_name'   : 'pack ice speed frequency',
                     'description' : 'sum of sea ice speed per grid cell over time divided by max speed and masked with pack ice criteria',
                     'start_time'  : self.dt0_str,
                     'stop_time'   : self.dtN_str})
        PIFZL_SD = (spatial_dims,
                    pifzl_sd.values,
                    {'units'       : '% (1/day)',
                     'long_name'   : 'pack ice frazil rate frequency',
                     'description' : 'sum of sea ice frazil rate per grid cell over time divided by max frazil growth rate and masked with pack ice criteria',
                     'start_time'  : self.dt0_str,
                     'stop_time'   : self.dtN_str})
        self.logger.info(f"\ttime taken: {time.time()-t0} seconds")
        #######################################
        #######      OUTPUT DATASET     #######
        #######################################
        self.logger.info("create output dataset:")
        t0 = time.time()
        PI = xr.Dataset(
            {# 1D variables
                'PIA'       : PIA,
                'PIV'       : PIV,
                'PI_STRONG' : PI_STRONG,
                'PI_SHEAR'  : PI_SHEAR,
                'PI_DIV'    : PI_DIV,
                'PI_AGING'  : PI_AGING,
                'PI_AGROM'  : PI_AGROM,
                'PI_AGROT'  : PI_AGROT,
                'PI_VGROM'  : PI_VGROM,
                'PI_VGROT'  : PI_VGROT,
                'PI_STRESS' : PI_STRESS,
                'PI_SPEED'  : PI_SPEED,
                'PI_FRAZIL' : PI_FRAZIL,
            # 2D variables
                'PIC_SD'   : PIC_SD,
                'PIHI_SD'  : PIHI_SD,
                'PISTH_SD' : PISTH_SD,
                'PISH_SD'  : PISH_SD,
                'PIDIV_SD' : PIDIV_SD,
                'PIAG_SD'  : PIAG_SD,
                'PIAD_SD'  : PIAD_SD,
                'PIAT_SD'  : PIAT_SD,
                'PIVD_SD'  : PIVD_SD,
                'PIVT_SD'  : PIVT_SD,
                'PISTR_SD' : PISTR_SD,
                'PISPD_SD' : PISPD_SD,
                'PIFZL_SD' : PIFZL_SD,
            # 3D variables
                'PIC'   : PIC,
                'PIHI'  : PIHI,
                'PISTH' : PISTH,
                'PISH'  : PISH,
                'PIDIV' : PIDIV,
                'PIAG'  : PIAG,
                'PIAD'  : PIAD,
                'PIAT'  : PIAT,
                'PIVD'  : PIVD,
                'PIVT'  : PIVT,
                'PISTR' : PISTR,
                'PISPD' : PISPD,
                'PIFZL' : PIFZL
            },
            coords = {time_dim        : (time_dim     , time_coords),
                      FI_time_dim     : (FI_time_dim  , PI_time_coords),
                      lon_coord_name  : (spatial_dims , lon_coords),
                      lat_coord_name  : (spatial_dims , lat_coords) }
        )
        PI.attrs = {"title"              : "Pack ice analysed from numerical sea ice model simulations",
                    "summary"            : "This dataset includes sea ice variables and derived metrics "\
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
                    "pack_ice_criteria"  : f"aice > {self.SIC_thresh} and not fast_ice_criteria",
                    "landmask_file"      : self.gi_processor.KMT_path,
                    "time_coverage_start": self.dt0_str,
                    "time_coverage_end"  : self.dtN_str,
                    "roll_window_days"   : roll_win,
                    "geospatial_lat_min" : float(np.min(lat_coords)),
                    "geospatial_lat_max" : float(np.max(lat_coords)),
                    "geospatial_lon_min" : float(np.min(lon_coords)),
                    "geospatial_lon_max" : float(np.max(lon_coords))}
        self.logger.info(f"\ttime taken: {time.time()-t0} seconds")
        return PI

    def process_window(self, dtC, hemisphere='south', ow_zarrs=False, save_zarr=True):
        """
        High-level method to process a centered time window and save Zarr output.

        Parameters
        ----------
        dtC       : datetime
                    Center date of the rolling window.
        hemisphere: str, optional
                    Hemisphere to process ('north' or 'south'). Default is 'south'.
        ow_zarrs  : bool, optional
                    If True, overwrite existing Zarr outputs. Default is False.
        save_zarr : bool, optional
                    Whether to save output to disk as a Zarr store. Default is True.

        Returns
        -------
        PI : xarray.Dataset
             Dataset with processed pack ice metrics.
        """
        import dask
        dask.config.set(scheduler='single-threaded')
        self.dtC = dtC
        self.logger.info(f"\n\nProcessing PI window centered on {self.dtC} for {hemisphere}ern hemisphere")
        self.define_hemisphere(hemisphere)
        ds = self.load_data_window()
        self.gi_processor.load_grid_and_landmask()
        roll_avg_dict     = self.compute_rolling_averages(ds)
        PI_masked_vars, _ = self.apply_PI_mask(roll_avg_dict)
        PI                = self.compute_pack_ice_outputs(PI_masked_vars)
        if save_zarr:
            D_PI_zarr = Path(self.config['D_dict']['AFIM_out'], self.sim_name, "PI")
            F_PI_zarr = f"pack_ice_{self.dtC.strftime('%Y-%m-%d')}.zarr"
            if not os.path.exists(D_PI_zarr):
                os.makedirs(D_PI_zarr)
            P_PI_zarr = Path(D_PI_zarr,F_PI_zarr)
            if not os.path.exists(P_PI_zarr) or (os.path.exists(P_PI_zarr) and ow_zarrs):
                self.logger.info(f"*** writing PI to disk: {P_PI_zarr}")
                t0 = time.time()
                PI.to_zarr(P_PI_zarr, mode='w')
                self.logger.info(f"\ttime taken: {time.time()-t0} seconds")
            else:
                self.logger.info("PI already exists or over-writting zarrs disabled ***")
                self.logger.info(f"\tskipping: {P_PI_zarr}")
        self.logger.info("Pack ice processing complete.")
        return PI  # placeholder for final dataset creation
