import xarray as xr
import pandas as pd
import numpy as np
from pathlib import Path
from datetime import datetime, timedelta
import xesmf as xe
from scipy.spatial import cKDTree
import json
import os, shutil
import sys
import time
import logging
import dask
sys.path.insert(0, '/home/581/da1339/AFIM/src/AFIM/src')
from grounded_iceberg_processor import GroundedIcebergProcessor 

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
    def __init__(self, sim_name, extra_cice_vars=None, hemisphere=None, P_log=None, json_path=None,
                 overwrite_gridded_climatology=False):
        """
        !!! GETTING STRAIGHT FROM BEGINNING IS SUPER IMPORTANT !!!
        Intended use is for dt0_str and dtN_str to define a long period to run the analysis; then use process_window() method to
        break this up into smaller self.roll_win periods
        """
        self.sim_name               = sim_name
        self.ow_gridded_climatology = overwrite_gridded_climatology
        if json_path is None:
            json_path = '/home/581/da1339/AFIM/src/AFIM/src/JSONs/afim_cice_analysis.json'
        with open(json_path, 'r') as f:
            self.config = json.load(f)
        self.doy_vals       = [1,16,31,46,61,76,91,106,121,136,151,166,181,196,211,226,241,256,271,286,301,316,331,346]
        self.cice_vars_reqd = self.config['CICE_dict']["FI_cice_vars_reqd"]
        self.cice_vars_ext  = extra_cice_vars if extra_cice_vars is not None else self.config['CICE_dict']["FI_cice_vars_ext"]
        self.cice_var_list  = self.cice_vars_reqd + self.cice_vars_ext
        self.roll_win       = self.config.get('roll_win', 15)
        if P_log is None:
            P_log = Path(self.config['D_dict']['logs'], f'FastIceProcessor_{sim_name}.log')
        self.setup_logging(logfile=P_log)
        self.sim_config = self.config['sim_dict'][sim_name]
        self.sim_dir = Path(self.config['D_dict']['AFIM_out'], sim_name, 'history', 'daily')
        self.gi_processor = GroundedIcebergProcessor(P_json=json_path, sim_name=sim_name)
        self.gi_processor.load_grid_and_landmask()
        self.use_gi = self.gi_processor.use_gi
        if self.use_gi:
            self.gi_processor.load_AFIM_GI()
        self.GI_total_area = self.gi_processor.total_area if self.use_gi else 0
        hemisphere = hemisphere if hemisphere is not None else self.config.get('hemisphere', 'south')
        self.define_hemisphere(hemisphere)
        self.FI_thresh = self.config.get('FI_thresh', 0.0005)
        self.FIC_scale = self.config.get('FIC_scale', 1e9)
        self.SIC_thresh = self.config.get('SIC_thresh', 0.15)
        self.SIC_scale = self.config.get('SIC_scale', 1e12)
        self.cm2m_fact = self.config.get('cm2m_fact', 0.01)
        self.CICE_dict = self.config['CICE_dict']
        self.regrid_weights_path = self.CICE_dict['P_reG_u2t_weights']

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
            self.hemisphere                   = 'north'
        else:
            raise ValueError(f"Invalid hemisphere '{hemisphere}'. Valid options are: "
                             "['north', 'south', 'northern', 'southern', 'sh', 'nh', 'SH', 'NH']")

    def slice_hemisphere(self, var_dict):
        return { k: v.isel(nj=self.hemisphere_nj_slice) if k not in {'FI_OBS_CLI'} else v for k, v in var_dict.items() }

    def _extract_cice_vars_from_PI_var_dict(self, FI_var_list):
        fi_meta = self.config["FI_var_dict"]
        FI_var_set_out = set()
        for v in FI_var_list:
            meta = fi_meta.get(v, {})
            cice_var = meta.get("CICE_variable")
            vec_vars = meta.get("CICE_vector_variables", [])
            # Skip derived output variables that aren't actually in the model dataset
            if cice_var and cice_var not in ['speed', 'strint', 'strair', 'strocn', 'strtlt', 'strcor']:
                FI_var_set_out.add(cice_var)
            # Always add vector components
            if isinstance(vec_vars, list):
                FI_var_set_out.update(vec_vars)
        # Always include aice (used in masking and core metrics)
        FI_var_set_out.add("aice")
        self.logger.debug(f"🧾 CICE variables required for computation: {sorted(FI_var_set_out)}")
        return sorted(FI_var_set_out)

    def load_obs_climatology(self, doy_start):
        csv_path = self.config['sea_ice_dict']['P_AF2020_cli_csv']
        df = pd.read_csv(csv_path)
        closest_doy = min(self.doy_vals, key=lambda x: abs(x - doy_start))
        row = df[df['DOY_start'] == closest_doy].copy()
        row = row.rename(columns={'Circumpolar': 'circumpolar',
                                  'IOsector': 'IOsector',
                                  'WPOsector': 'WPOsector',
                                  'RSsector': 'RSsector',
                                  'BASsector': 'BASsector',
                                  'WSsector': 'WSsector'})
        sectors = ['circumpolar', 'IOsector', 'WPOsector', 'RSsector', 'BASsector', 'WSsector']
        return row[sectors].reset_index(drop=True)

    def _load_obs_gridded(self, P_orgs):
        FI_obs = xr.open_mfdataset(P_orgs,
                                   combine='nested',
                                   concat_dim='time',
                                   parallel=True,
                                   chunks='auto',
                                   engine='netcdf4')
        FI_reG, lon_reG, lat_reG = self.AFdb_regrid_to_tgrid(FI_obs)
        mask     = xr.where(FI_reG >= 4, 1, np.nan)
        FI       = (('t_FI_obs', 'nj', 'ni'),
                    FI_reG.where(mask).data,
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
                    lon_reG.values,
                    {'long_name': 'longitude',
                     'units'    : 'degrees_north'})
        y_coords = (('nj','ni'),
                    lat_reG.values,
                    {'long_name': 'latitude',
                     'units'    : 'degrees_east'})
        return xr.Dataset({'FI'       : FI,
                           'FI_t_alt' : t_alt },
                          coords=dict(t_FI_obs=t_coords, obs_lon=x_coords, obs_lat=y_coords))

    def AFdb_regrid_to_tgrid(self, FI_obs_native):
        from pyproj import CRS, Transformer
        self.logger.info("*** Regridding 'AF_FI_OBS_2020db' to CICE T-grid *** ")
        self.logger.info("STEP 1: convert 'AF_FI_OBS_2020db' Cartesian coordinates to spherical coordinates")
        crs_nsidc   = CRS.from_epsg(3412)
        crs_wgs84   = CRS.from_epsg(4326)
        transformer = Transformer.from_crs(crs_nsidc, crs_wgs84, always_xy=True)
        X, Y        = np.meshgrid(FI_obs_native['x'].isel(time=0).values, FI_obs_native['y'].isel(time=0).values)
        lon, lat    = transformer.transform(X,Y)
        self.logger.info("STEP 2: load regridder into memory or create if it does not exist")
        FI_OBS_GRD  = xr.Dataset({ 'lon' : (('y', 'x'), lon),
                                   'lat' : (('y', 'x'), lat)})
        t1 = time.time()
        if os.path.exists(self.config["sea_ice_dict"]["AF_reG_weights"]):
            regridder = xe.Regridder(FI_OBS_GRD, self.gi_processor.G_t, method="bilinear", periodic=True, weights=self.config["sea_ice_dict"]["AF_reG_weights"])
        else:
            regridder = xe.Regridder(FI_OBS_GRD, self.gi_processor.G_t, method="bilinear", periodic=True, filename=self.config["sea_ice_dict"]["AF_reG_weights"])
        self.logger.info(f"\ttime taken: {time.time()-t1:0.2f} seconds")
        self.logger.info("STEP 3: regrid longitudes, latitudes and 'Fast_Ice_Time_series' arrays")
        lon_reG = regridder(FI_OBS_GRD["lon"])
        lat_reG = regridder(FI_OBS_GRD["lat"])
        FI_reG  = regridder(FI_obs_native["Fast_Ice_Time_series"])
        return FI_reG, lon_reG, lat_reG

    def filter_obs_gridded_by_date(self, start_date, end_date):
        D_obs    = Path(self.config['sea_ice_dict']['D_AF2020_db_org'])
        yrs_reqd = set([start_date.year, end_date.year])
        P_orgs   = [D_obs / f"FastIce_70_{yr}.nc" for yr in yrs_reqd]
        ds       = self._load_obs_gridded(P_orgs)
        # Convert FI_t_alt (int like 20000101) to datetime
        alt_dates = pd.to_datetime(ds['FI_t_alt'].values.astype(str), format='%Y%m%d')
        # Use where logic to find which obs period includes start_date
        valid_idx = (alt_dates >= pd.to_datetime(start_date)) & (alt_dates <= pd.to_datetime(end_date))
        matched = ds.sel(t_FI_obs=valid_idx)
        if matched.dims['t_FI_obs'] == 0:
            self.logger.warning(f"No matching observational dates found between {start_date} and {end_date}")
        return matched

    def create_obs_gridded_climatology(self):
        P_zarr = Path(self.config['sea_ice_dict']["P_AF_2020db_avg"])
        if P_zarr.exists() and not self.ow_gridded_climatology:
            self.logger.info(f"Averaged observational gridded climatology already exists\n{P_zarr}")
            return xr.open_zarr(P_zarr)
        elif P_zarr.exists() and self.ow_gridded_climatology:
            self.logger.info(f"Averaged observational gridded climatology exists *but* over-writing has been requested\n\tOVER-WRITING: {P_zarr}")
            shutil.rmtree(P_zarr)
        else:
            self.logger.info(f"Averaged observational gridded climatology does *NOT* exist\n\tCREATING NEW: {P_zarr}")
        # Load all gridded obs, already regridded and with date_alt available
        D_obs = Path(self.config['sea_ice_dict']['D_AF2020_db_org'])
        P_orgs = sorted(D_obs.glob("FastIce_70_*.nc"))
        self.logger.info(f"loading all observational fast ice netcdf files:\n{P_orgs}")
        ds_all = self._load_obs_gridded(P_orgs)
        # Convert alt date to datetime
        alt_dates = pd.to_datetime(ds_all['FI_t_alt'].values.astype(str), format='%Y%m%d')
        doy_vals = alt_dates.dayofyear
        ds_all = ds_all.assign_coords(doy=("t_FI_obs", doy_vals))
        # Align to the official 24 periods only
        ds_all = ds_all.sel(t_FI_obs=[i for i, doy in enumerate(doy_vals) if doy in self.doy_vals])
        # Group by DOY and average
        grouped = ds_all.groupby("doy").mean(dim="t_FI_obs", skipna=True)
        # Rename DOY dim to something more descriptive
        grouped = grouped.rename_dims({"doy" : "t_doy"})
        grouped = grouped.rename_vars({"FI"  : "FI_OBS_GRD"})
        grouped = grouped.assign_attrs(doy_vals=self.doy_vals)
        self.logger.info(f"gridded climatology now looks like\n{grouped}")
        self.logger.info(f"writing this to {P_zarr}")
        t1 = time.time()
        grouped.to_zarr(P_zarr)
        self.logger.info(f"\tzarr written in {time.time()-t1:0.2f} seconds")
        return grouped

    def load_data_window(self, window_start):
        """
        Load 15-day window of CICE data and attach observational climatology and gridded data.
        `doy_start` is matched to the nearest 15-day or final 20-day climatology period in self.doy_vals.
        """
        if self.associate_observational:
            window_end = window_start + pd.Timedelta(days=self.roll_win - 1)
            self.logger.info(f"loading data for date period: {window_start} to {window_end}")
            # Load observational climatology for the matching DOY window
            doy_start = int(window_start.strftime('%j'))
            self.logger.info(f"loading CSV climatology")
            cli_df = self.load_obs_climatology(doy_start)
            sectors = cli_df.columns.tolist()
            # Load or compute gridded obs if within range, else use climatology
            if pd.Timestamp('2000-03-01') <= window_start <= pd.Timestamp('2018-03-31'):
                self.logger.info("using gridded observations")
                obs_gridded = self.filter_obs_gridded_by_date(window_start, window_end)
                obs_grd_data = obs_gridded['FI'].mean(dim='t_FI_obs', skipna=True)
            else:
                self.logger.info("using gridded climatology")
                clim_all = self.create_obs_gridded_climatology()
                closest_doy = min(self.doy_vals, key=lambda x: abs(x - doy_start))
                self.logger.info(f"using {closest_doy}-DOY from gridded climatology to associate with this model period")
                doy_index = np.argmin(np.abs(np.array(self.doy_vals) - closest_doy))
                obs_grd_data = clim_all['FI_OBS_GRD'].isel(t_doy=doy_index)
        self.logger.info("constructing list of CICE files to load for this period")
        P_CICE_orgs = [self.sim_dir / f"iceh.{d.strftime('%Y-%m-%d')}.nc"
                       for d in self.dt_period_list
                       if (self.sim_dir / f"iceh.{d.strftime('%Y-%m-%d')}.nc").exists()]
        if not P_CICE_orgs:
            raise FileNotFoundError(f"No CICE files found for window around {window_start}")
        self.logger.info("Loading CICE files ...")
        self.logger.debug(f"{P_CICE_orgs}")
        t1 = time.time()
        preprocess = lambda ds: ds[list(self.cice_var_list)]
        CICE = xr.open_mfdataset(P_CICE_orgs, combine='by_coords', parallel=True, preprocess=preprocess)
        self.logger.info(f"\t✅ Model dataset loaded: shape {CICE.sizes}, time: {time.time() - t1:.2f} s")
        if self.associate_observational:
            self.logger.debug(f"climatology dataframe from CSV looks like\n{cli_df}")
            if cli_df.shape[0] == 1:
                data = cli_df.to_numpy()
            else:
                raise ValueError("Expected a single row from climatology CSV.")
            CICE['FI_OBS_CLI'] = xr.DataArray( data,
                                               dims=['t_doy', 'sector'],
                                               coords={'t_doy'  : [doy_start],
                                                       'sector' : cli_df.columns.tolist()} )
            CICE['FI_OBS_CLI'].attrs['doy'] = closest_doy
            CICE['FI_OBS_GRD']              = obs_grd_data
            CICE['FI_OBS_GRD'].attrs['doy'] = closest_doy
        if 'doy' in CICE.coords:
            CICE = CICE.drop_vars('doy')
        self.logger.debug(f"returning load_data_window() method with a dataset that looks like\n{CICE}")
        return CICE

    def CICE_regrid_to_tgrid(self, ds):
        self.logger.info("Regridding CICE 'uvel' and 'vvel' (sea ice velocity components) to T-grid...")
        t1 = time.time()
        regridder = xe.Regridder(self.gi_processor.G_u, self.gi_processor.G_t, method="bilinear", extrap_method="inverse_dist", periodic=True, weights=self.regrid_weights_path)
        self.logger.info(f"\ttime taken: {time.time()-t1:0.2f} seconds")
        u_reG = regridder(ds["uvel"])
        v_reG = regridder(ds["vvel"])
        return u_reG, v_reG

    def dataset_to_dictionary(self, ds):
        fi_meta = self.config["FI_var_dict"]
        CICE_dict_unrolled = {}
        for v in self.cice_vars_reqd:
            if v in ds:
                CICE_dict_unrolled[v] = ds[v]
        uice, vice = self.CICE_regrid_to_tgrid(ds)
        CICE_dict_unrolled['speed'] = np.sqrt(uice ** 2 + vice ** 2)
        for out_var in fi_meta:
            meta = fi_meta.get(out_var, {})
            cice_var = meta.get("CICE_variable")
            vec_vars = meta.get("CICE_vec_vars", [])
            if vec_vars and cice_var:
                if all(vv in ds for vv in vec_vars):
                    # Ensure the derived result is a DataArray
                    base_da = ds[vec_vars[0]]
                    derived = np.sqrt(np.sum([ds[vv] ** 2 for vv in vec_vars], axis=0))
                    CICE_dict_unrolled[cice_var] = xr.DataArray(derived,
                                                                dims=base_da.dims,
                                                                coords=base_da.coords,
                                                                attrs={"long_name": f"derived variable from {vec_vars}"} )
                else:
                    missing = [vv for vv in vec_vars if vv not in ds]
                    self.logger.debug(f"⚠️ Skipping vector-derived var {out_var} — missing components: {missing}")
            elif cice_var and cice_var in ds:
                CICE_dict_unrolled[cice_var] = ds[cice_var]
        self.logger.info(f"compute rolling mean for {self.roll_win}-days on CICE variables:\n{list(CICE_dict_unrolled.keys())}")
        roll = lambda da: da.rolling(time=self.roll_win, center=True, min_periods=1).mean()
        CICE_dict_rolled = {k: roll(v) for k, v in CICE_dict_unrolled.items()}
        coarse = lambda da: da.coarsen(time=self.roll_win, boundary="trim").mean()
        CICE_dict_coarsened = {k: coarse(v) for k, v in CICE_dict_rolled.items()}
        if self.associate_observational:
            CICE_dict_coarsened['FI_OBS_CLI'] = ds['FI_OBS_CLI']
            CICE_dict_coarsened['FI_OBS_GRD'] = ds['FI_OBS_GRD']
        return CICE_dict_coarsened

    def apply_FI_mask(self, roll_dict):
        fi_mask      = (roll_dict['aice'] > self.SIC_thresh) & (roll_dict['speed'] <= self.FI_thresh)
        masked_dict  = { k: v.where(fi_mask) if k not in {'FI_OBS_CLI', 'FI_OBS_GRD'} else v for k, v in roll_dict.items() }
        self.fi_mask = fi_mask.isel(nj=self.hemisphere_nj_slice)
        return masked_dict

    def compute_fast_ice_outputs(self, cice_vars_dict):
        self.logger.info("\n *** COMPUTING FAST ICE OUTPUTS ***")
        self.logger.info("STEP 1:\n\tSLICE HEMISPHERE")
        t1 = time.time()
        cice_vars_dict = self.slice_hemisphere(cice_vars_dict)
        self.logger.info(f"\ttime taken: {time.time()-t1:0.2f} seconds")
        self.logger.info("STEP 2:\n\tMETHOD DEFINITIONS")
        grid_cell_area = self.gi_processor.G_t['area'].isel(nj=self.hemisphere_nj_slice)
        cm2m           = self.cm2m_fact
        fi_meta        = self.config["FI_var_dict"]
        one_d_list     = [k for k, meta in fi_meta.items() if meta.get("dimensions") == "1D"]
        two_d_list     = [k for k, meta in fi_meta.items() if meta.get("dimensions") == "2D"]
        three_d_list   = [k for k, meta in fi_meta.items() if meta.get("dimensions") == "3D"]
        CICE_time_dim  = 'time'
        t_dim_str      = 't_dim'
        x_dim_str      = 'ni'
        y_dim_str      = 'nj'
        sector_dim_str = 'sector'
        t_coord_str    = t_dim_str
        x_coord_str    = 'lon'
        y_coord_str    = 'lat'
        one_dim_tup    = (t_dim_str,)
        two_dim_tup    = (y_dim_str, x_dim_str)
        three_dim_tup  = (t_dim_str, y_dim_str, x_dim_str)
        sector_dim_tup = (t_dim_str, sector_dim_str)
        x_coord_tup    = (two_dim_tup, self.gi_processor.G_t['lon'].isel(nj=self.hemisphere_nj_slice).values)
        y_coord_tup    = (two_dim_tup, self.gi_processor.G_t['lat'].isel(nj=self.hemisphere_nj_slice).values)
        t_coord_tup    = (one_dim_tup, [self.dt0_period]) # this should be the start of the {self.roll_win}-day period
                                                        # this ensures consistency with observational fast ice
        sector_coord_tup = ((sector_dim_str), cice_vars_dict['FI_OBS_CLI'].sector.values)
        #################################
        #######   1D VARIABLES   ########
        #################################
        self.logger.info("STEP 3:\n\tCOMPUTE 1D OUTPUT VARIABLES")
        t1 = time.time()
        one_d_metrics = {}
        for v in one_d_list:
            if v.endswith("_SD") or (v in ['FIA_CLI']):
                continue
            meta     = fi_meta.get(v, {})
            cice_var = meta.get("CICE_variable", None)
            self.logger.debug(f"\t📦 Available keys in cice_vars_dict: {list(cice_vars_dict.keys())}")
            if not cice_var or cice_var not in cice_vars_dict.keys():
                self.logger.warning(f"\t⚠️ Skipping 1D metric {v} — source variable '{cice_var}' missing.")
                continue
            else:
                self.logger.debug(f"\t✅ Creating 1D metric {v} — source variable '{cice_var}'")
            if v == "FIA":
                one_d_metrics[v] = ((cice_vars_dict[cice_var] * grid_cell_area).sum(dim=two_dim_tup)) / self.FIC_scale
            elif v == "FIE" and hasattr(self, "fi_mask"):
                one_d_metrics[v] = ((self.fi_mask * grid_cell_area).sum(dim=two_dim_tup)) / self.FIC_scale
            elif "AGING" in v:
                one_d_metrics[v] = (grid_cell_area / cice_vars_dict[cice_var]).sum(dim=two_dim_tup)
            elif "VGRO" in v or "FRAZIL" in v:
                one_d_metrics[v] = (cice_vars_dict[cice_var] * cm2m * grid_cell_area).sum(dim=two_dim_tup)
            else:
                one_d_metrics[v] = (cice_vars_dict[cice_var] * grid_cell_area).sum(dim=two_dim_tup)
        one_d_vars            = {k: xr.DataArray(data   = v.data,
                                                 dims   = one_dim_tup,
                                                 coords = {t_dim_str : t_coord_tup},
                                                 attrs  = fi_meta.get(k, {}))
                                 for k, v in one_d_metrics.items()}
        if self.associate_observational:
            one_d_vars['FIA_OBS'] = xr.DataArray(data   = cice_vars_dict['FI_OBS_CLI'].data,
                                                 dims   = sector_dim_tup,
                                                 coords = {t_dim_str      : t_coord_tup,
                                                           sector_dim_str : sector_coord_tup},
                                                 attrs  = fi_meta.get('FIA_OBS', {}))
        self.logger.info(f"\ttime taken: {time.time()-t1:0.2f} seconds")
        self.logger.debug(f"\t1D vars computed: {list(one_d_vars.keys())}")
        #################################
        #######   3D VARIABLES   ########
        #################################
        self.logger.info("STEP 3:\n\tCOMPUTE 3D OUTPUT VARIABLES:")
        self.logger.info("\t very little processing done and essentially CICE variables are only filtered/masked for fast ice criteria")
        t1 = time.time()
        three_d_vars = {}
        for v in three_d_list:
            if v.endswith("_SD") or (v in ['FI_GRD']):
                continue
            meta     = fi_meta.get(v, {})
            cice_var = meta.get("CICE_variable", None)
            self.logger.debug(f"\t📦 Available keys in cice_vars_dict: {list(cice_vars_dict.keys())}")
            if not cice_var or cice_var not in cice_vars_dict.keys():
                self.logger.warning(f"\t⚠️ Skipping 3D metric {v} — source variable '{cice_var}' missing.")
                continue
            else:
                self.logger.debug(f"\t✅ Creating 3D metric {v} — source variable '{cice_var}'")
                data = cice_vars_dict[cice_var].data
                if "VGRO" in v or "FZL" in v:
                    data = data * cm2m
                if data.ndim == 2:
                    dim_name, coord_vals = t_coord_tup
                    data = data.expand_dims({dim_name: coord_vals})
                three_d_vars[v] = xr.DataArray(data   = data,
                                               dims   = three_dim_tup,
                                               coords = {t_dim_str   : t_coord_tup,
                                                         x_coord_str : x_coord_tup,
                                                         y_coord_str : y_coord_tup},
                                               attrs  = meta)
        if self.associate_observational:
            FI_OBS = cice_vars_dict['FI_OBS_GRD'].isel(nj=self.hemisphere_nj_slice).expand_dims({t_dim_str: [self.dt0_period]}).drop_vars(['TLON','TLAT','ULON','ULAT','obs_lon','obs_lat'])
            three_d_vars['FI_OBS'] = xr.DataArray(data   = FI_OBS.data,
                                                  dims   = three_dim_tup,
                                                  coords = {t_dim_str   : t_coord_tup,
                                                            x_coord_str : x_coord_tup,
                                                            y_coord_str : y_coord_tup},
                                                  attrs  = fi_meta.get('FI_OBS', {}))
        self.logger.info(f"\ttime taken: {time.time()-t1:0.2f} seconds")
        self.logger.debug(f"\t3D vars computed: {list(three_d_vars.keys())}")
        #################################
        #######   2D VARIABLES   ########
        #################################
        self.logger.info("STEP 5:\n\tCOMPUTE 2D OUTPUT VARIABLES:")
        self.logger.info("\t significant temporal averaging done to compute this portion of dataset ... can take a little bit of time")
        t1 = time.time()
        two_d_vars = {}
        for v in two_d_list:
            if v=='FIP_OBS':
                continue
            meta     = fi_meta.get(v, {})
            cice_var = meta.get("CICE_variable", None)
            self.logger.debug(f"\t📦 Available keys in cice_vars_dict: {list(cice_vars_dict.keys())}")
            if not cice_var or cice_var not in cice_vars_dict.keys():
                self.logger.warning(f"\t⚠️ Skipping 2D var {v} due to missing base variable '{cice_var}'")
                continue
            else:
                self.logger.debug(f"\t✅ Creating 2D metric {v} — source variable '{cice_var}'")
            da       = cice_vars_dict[cice_var]
            norm     = da.sizes[CICE_time_dim]
            data_sum = da.sum(dim=CICE_time_dim)
            if v=='FIP':
                data_mean = data_sum / norm
            else:
                max_val   = da.max()
                data_mean = data_sum / (norm * max_val)
            two_d_vars[v] = xr.DataArray(data   = data_mean.data,
                                         dims   = two_dim_tup,
                                         coords = {x_coord_str : x_coord_tup,
                                                   y_coord_str : y_coord_tup},
                                         attrs  = {**meta,
                                                   "start_time": self.dt0_period.isoformat() if isinstance(self.dt0_period, datetime) else self.dt0_period,
                                                   "stop_time" : self.dtN_period.isoformat() if isinstance(self.dtN_period, datetime) else self.dtN_period})
        if self.associate_observational:
            FIP_OBS = cice_vars_dict['FI_OBS_GRD'].isel(nj=self.hemisphere_nj_slice).drop_vars(['TLON','TLAT','ULON','ULAT','obs_lon','obs_lat'])
            two_d_vars['FIP_OBS'] = xr.DataArray(data   = FIP_OBS.data,
                                                 dims   = two_dim_tup,
                                                 coords = {x_coord_str : x_coord_tup,
                                                           y_coord_str : y_coord_tup},
                                                 attrs  = {**fi_meta.get('FIP_OBS', {}),
                                                           "start_time": self.dt0_period.isoformat() if isinstance(self.dt0_period, datetime) else self.dt0_period,
                                                           "stop_time" : self.dtN_period.isoformat() if isinstance(self.dtN_period, datetime) else self.dtN_period})
        self.logger.info(f"\ttime taken: {time.time()-t1:0.2f} seconds")
        self.logger.debug(f"\t2D vars computed: {list(three_d_vars.keys())}")
        #######################################
        #######      OUTPUT DATASET     #######
        #######################################
        self.logger.info("STEP 6:\n\tCREATE OUTPUT DATASET")
        t1 = time.time()
        FI = xr.Dataset({**one_d_vars, **three_d_vars, **two_d_vars})
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
                    "roll_window_days"   : self.roll_win,
                    "geospatial_lat_min" : float(np.min(self.gi_processor.G_t['lat'].isel(nj=self.hemisphere_nj_slice).values)),
                    "geospatial_lat_max" : float(np.max(self.gi_processor.G_t['lat'].isel(nj=self.hemisphere_nj_slice).values)),
                    "geospatial_lon_min" : float(np.min(self.gi_processor.G_t['lon'].isel(nj=self.hemisphere_nj_slice).values)),
                    "geospatial_lon_max" : float(np.max(self.gi_processor.G_t['lon'].isel(nj=self.hemisphere_nj_slice).values))}
        # make sure 'time' dimension is scrubbed from output dataset
        FI = FI.map(lambda da: da.squeeze(dim='time') if 'time' in da.dims else da)
        FI = FI.map(lambda da: da.drop_dims('time') if 'time' in da.dims else da)
        # Drop 'time' coordinate entirely if it exists (from coords, not just dims)
        if 'time' in FI.coords:
            FI = FI.drop_vars('time')
        # Forcefully remove 'time' from any lingering data variable dimensions
        for var in FI.data_vars:
            if 'time' in FI[var].dims:
                FI[var] = FI[var].squeeze(dim='time', drop=True)
        # (Optional sanity check)
        assert 'time' not in FI.dims
        self.logger.info(f"\ttime taken: {time.time()-t1:0.2f} seconds")
        return FI

    def process_window(self, dt0_str=None, dtN_str=None, rolling_window=None, associate_observational=True,
                       write_zarr=False, ow_zarrs=False):
        """
        Default behaviour is to assume user is interactively creating FI datasets
        for a particular simulation and time period. Allows a user to test rolling
        mean period and how this effects the dataset. However, this *TURNS OFF*
        any and all loading of the observational climatology and observational gridded
        datasets. This is because the observational datasets are somewhat 'set in stone',
        so to speak. The observational datasets are based on 15-day periods starting on
        1 January to the 346th year day and then uses a 20-day period. At present
        this rigidity in the observational temporal dimension puts constraints on
        dynamically associating CICE model data with it.

        However, for now if one uses process_window() in the following manner
        dt0_str  = "1993-09-01"
        dtN_str  = "1994-12-31"
        sim_name = 'ktens-max'
        FI_proc  = FastIceProcessor(sim_name = sim_name)
        FI       = FI_proc.process_window(dt0_str = dt0_str,
                                          dtN_str = dtN_str)
        Then one can expect to get FI with the default rolling window 15 days up to the last 20-day
        period *and* expect to get the observational data associated with output. They may
        choose to write it disk under this condition. 

        One can provide optional rolling_window integers to test the effect on their dataset. However,
        should rolling_window be set to anything other than 15 then loading associating with
        observational fast ice datasets will be switched off, for reasons given above. A user can
        write to zarr files under this condition if they so choose but the author advises against
        this usage unless the user has gone through the effort of ensuring they are
        using a non-default JSON file that clearly saves these FI output datasets in a non-default
        location. 
        """
        dask.config.set(scheduler='single-threaded')
        self.dt0_str = dt0_str if dt0_str is not None else self.config.get("dt0_str", "1993-01-01")
        self.dtN_str = dtN_str if dtN_str is not None else self.config.get("dtN_str", "1993-12-31")
        self.dt0 = pd.Timestamp(self.dt0_str)
        self.dtN = pd.Timestamp(self.dtN_str)
        self.roll_win = rolling_window if rolling_window is not None else self.config.get("roll_win", 15)
        self.associate_observational = associate_observational
        if rolling_window and rolling_window != 15:
            self.associate_observational = False

        # Define static list of DOY start values for 15-day periods plus final 20-day period
        self.doy_vals = list(range(1, 362, 15)) + [347]

        # Construct aligned dtC_list to climatology periods
        dtC_list = []
        for year in range(self.dt0.year, self.dtN.year + 1):
            for doy in self.doy_vals:
                dt = pd.Timestamp(datetime(year, 1, 1)) + pd.Timedelta(days=doy - 1)
                if self.dt0 <= dt <= self.dtN:
                    dtC_list.append(dt)

        ds_all = []
        for dtC in dtC_list:
            doy = dtC.timetuple().tm_yday
            is_leap = dtC.is_leap_year if hasattr(dtC, 'is_leap_year') else (dtC.year % 4 == 0 and (dtC.year % 100 != 0 or dtC.year % 400 == 0))
            last_doy = 366 if is_leap else 365
            roll_win = 20 if doy >= 347 else 15

            self.roll_win = roll_win
            self.dtC = dtC
            self.dt0_period = dtC - pd.Timedelta(days=self.roll_win // 2)
            self.dtN_period = dtC + pd.Timedelta(days=self.roll_win // 2)
            self.dt_period_list = pd.date_range(self.dt0_period, self.dtN_period, freq='1D')

            ds = self.load_data_window(self.dt0_period)
            ds_dict = self.dataset_to_dictionary(ds)
            masked_vars = self.apply_FI_mask(ds_dict)
            FI = self.compute_fast_ice_outputs(masked_vars)

            if write_zarr:
                D_FI_zarr = Path(self.config['D_dict']['AFIM_out'], self.sim_name, "FI")
                F_FI_zarr = f"fast_ice_{self.dt0_period.strftime('%Y-%m-%d')}.zarr"
                P_FI_zarr = Path(D_FI_zarr, F_FI_zarr)
                if not D_FI_zarr.exists():
                    os.makedirs(D_FI_zarr)
                if not P_FI_zarr.exists() or ow_zarrs:
                    self.logger.info(f"*** writing FI dataset to disk: {P_FI_zarr}")
                    t1 = time.time()
                    FI.to_zarr(P_FI_zarr, mode='w')
                    self.logger.info(f"\ttime taken: {time.time()-t1:0.2f} seconds")
                else:
                    self.logger.info(f"*** FI dataset zarr file already exists and overwriting disabled:\n\t{P_FI_zarr}")

            ds_all.append(FI)

        FI_merged = xr.concat(ds_all, dim='time')
        if 'time' in FI_merged.coords:
            FI_merged = FI_merged.drop_vars('time')
        self.logger.info("✅ Pack ice processing complete.")
        return FI_merged
