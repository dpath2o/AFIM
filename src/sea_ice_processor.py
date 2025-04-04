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
#import dask
sys.path.insert(0, '/home/581/da1339/AFIM/src/AFIM/src')
from grounded_iceberg_processor import GroundedIcebergProcessor 

class SeaIceProcessor:
    """
    Class to compute sea ice metrics from sea ice model output.

    This processor applies a rolling window to compute spatial and temporal
    characteristics of sea ice (original intent was landfast sea ice) using
    sea ice concentration and velocity thresholds.

    Default is compute to fast ice as defined as sea ice concentration above a threshold and
    sea ice velocity below a threshold (i.e., stationary sea ice).

    Optionally can compute 'pack ice' that is either inclusive/exclusive of fast ice.

    Parameters
    ----------

    Attributes
    ----------

    Examples
    --------

    The SeaIceProcessor is typically run in a loop using a driver script
    such as [`compute_fast_ice.py`](https://github.com/dpath2o/AFIM/blob/main/src/python/compute_fast_ice.py):

    For more, see the full project repository:
    🔗 https://github.com/dpath2o/AFIM
    """
    def __init__(self, sim_name,
                 sea_ice                     = False,
                 pack_ice                    = False,
                 ice_concentration_threshold = None,
                 ice_speed_threshold         = None,
                 extra_cice_vars             = None,
                 hemisphere                  = None,
                 P_log                       = None,
                 P_json                      = None,
                 overwrite_AF2020_zarr       = False,
                 zarr_directory              = None):
        """
        """
        if P_json is None:
            P_json = '/home/581/da1339/AFIM/src/AFIM/src/JSONs/afim_cice_analysis.json'
        with open(P_json, 'r') as f:
            self.config = json.load(f)
        self.sim_name       = sim_name
        self.sim_config     = self.config['sim_dict'][sim_name]
        self.sim_dir        = Path(self.config['D_dict']['AFIM_out'], sim_name, 'history', 'daily')
        self.ispd_thresh    = ice_speed_threshold if ice_speed_threshold is not None else self.config.get('ice_speed_thresh', 0.0005)
        self.icon_thresh    = ice_concentration_threshold if ice_concentration_threshold is not None else self.config.get('ice_conc_thresh', 0.15)
        self.doy_vals       = self.config.get("DOY_vals",[1,16,31,46,61,76,91,106,121,136,151,166,181,196,211,226,241,256,271,286,301,316,331,346])
        self.CICE_dict      = self.config['CICE_dict']
        self.cice_vars_reqd = self.CICE_dict["FI_cice_vars_reqd"]
        self.cice_vars_ext  = extra_cice_vars if extra_cice_vars is not None else self.CICE_dict["FI_cice_vars_ext"]
        self.cice_var_list  = self.cice_vars_reqd + self.cice_vars_ext
        self.FIC_scale      = self.config.get('FIC_scale', 1e9)
        self.SIC_scale      = self.config.get('SIC_scale', 1e12)
        self.cm2m_fact      = self.config.get('cm2m_fact', 0.01)
        if sea_ice and pack_ice:
            self.logger.info("cannot enable both sea ice and pack ice output ... select one or the other ... defaulting to fast ice")
            sea_ice  = False
            pack_ice = False
        self.sea_ice  = sea_ice
        self.pack_ice = pack_ice
        D_log         = self.config['D_dict']['logs']
        if self.sea_ice:
            self.assoc_AF2020   = False
            self.ow_AF2020_zarr = False
            self.D_zarr_out     = zarr_directory if zarr_directory is not None else Path(self.config['D_dict']['AFIM_out'], self.sim_name, "SI")
            self.F_zarr_out_fmt = "sea_ice_{date_str:s}.zarr"
            P_log               = P_log if P_log is not None else Path(D_log, f'SeaIceProcessor_SI_{sim_name}.log')
        elif self.pack_ice:
            self.assoc_AF2020   = False
            self.ow_AF2020_zarr = False
            self.D_zarr_out     = zarr_directory if zarr_directory is not None else Path(self.config['D_dict']['AFIM_out'], self.sim_name, "PI")
            self.F_zarr_out_fmt = "pack_ice_{date_str:s}.zarr"
            P_log               = P_log if P_log is not None else Path(D_log, f'SeaIceProcessor_PI_{sim_name}.log')
        else:
            self.assoc_AF2020   = True
            self.sea_ice        = False
            self.pack_ice       = False
            self.ow_AF2020_zarr = overwrite_AF2020_zarr
            self.D_zarr_out     = zarr_directory if zarr_directory is not None else Path(self.config['D_dict']['AFIM_out'], self.sim_name, "FI")
            self.F_zarr_out_fmt = "fast_ice_{date_str:s}.zarr"
            P_log               = P_log if P_log is not None else Path(D_log, f'SeaIceProcessor_FI_{sim_name}.log')
        self.setup_logging(logfile=P_log)
        if not self.D_zarr_out.exists():
            os.makedirs(self.D_zarr_out)
        self.gi_processor = GroundedIcebergProcessor(P_json=P_json, sim_name=sim_name)
        self.gi_processor.load_grid_and_landmask()
        self.use_gi = self.gi_processor.use_gi
        if self.use_gi:
            self.gi_processor.load_AFIM_GI()
        self.GI_total_area = self.gi_processor.total_area if self.use_gi else 0
        hemisphere         = hemisphere if hemisphere is not None else self.config.get('hemisphere', 'south')
        self.define_hemisphere(hemisphere)

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
            self.hemisphere                   = 'south'
        else:
            raise ValueError(f"Invalid hemisphere '{hemisphere}'. Valid options are: "
                             "['north', 'south', 'northern', 'southern', 'sh', 'nh', 'SH', 'NH']")

    def slice_hemisphere(self, var_dict):
        return { k: v.isel(nj=self.hemisphere_nj_slice) if k not in {'FI_OBS_CLI'} else v for k, v in var_dict.items() }

    def _convert_NSIDC_cartesian_to_spherical(self, ds):
        from pyproj import CRS, Transformer
        time_dim = self.config['NSIDC_dict']['time_dim']
        x_dim    = self.config['NSIDC_dict']['x_dim']
        y_dim    = self.config['NSIDC_dict']['y_dim']
        x_coord  = self.config['NSIDC_dict']['x_coord']
        y_coord  = self.config['NSIDC_dict']['y_coord']
        self.logger.info("🧭 Converting NSIDC Cartesian to spherical coordinates:")
        t1 = time.time()
        crs_proj = CRS.from_proj4(self.config["NSIDC_dict"]["projection_string"])
        crs_wgs84 = CRS.from_epsg(4326)
        transformer = Transformer.from_crs(crs_proj, crs_wgs84, always_xy=True)
        x, y = ds[x_coord].values, ds[y_coord].values
        X, Y = np.meshgrid(x, y)
        lon, lat = transformer.transform(X, Y)
        ds['lon'] = ((y_dim, x_dim), lon)
        ds['lat'] = ((y_dim, x_dim), lat)
        ds = ds.swap_dims({time_dim: 'time'})
        self.logger.info(f"\t✅ NSIDC coordinate conversion complete in {time.time()-t1:.2f} seconds")
        return ds

    def _compute_NSIDC_SIA_SIE(self, ds):
        t1 = time.time()
        SIC_name = self.config["NSIDC_dict"]["SIC_name"]
        flags = self.config["NSIDC_dict"]["cdr_seaice_conc_flags"]
        y_dim = self.config["NSIDC_dict"]["y_dim"]
        x_dim = self.config["NSIDC_dict"]["x_dim"]
        aice = ds[SIC_name]
        for flag in flags:
            aice = xr.where(aice == flag / 100, np.nan, aice)
        area = xr.open_dataset(self.config["NSIDC_dict"]["P_cell_area"]).cell_area
        mask = aice > self.icon_thresh
        SIA = (aice * area).where(mask.notnull()).sum(dim=[y_dim, x_dim], skipna=True)
        SIE = (mask * area).sum(dim=[y_dim, x_dim], skipna=True)
        aice = aice.where(mask)
        SIA = SIA.rolling(time=self.roll_win, center=True).mean()
        SIE = SIE.rolling(time=self.roll_win, center=True).mean()
        SIA = SIA.coarsen(time=self.roll_win, boundary="trim").mean()
        SIE = SIE.coarsen(time=self.roll_win, boundary="trim").mean()
        #aice = aice.rolling(time=self.roll_win, center=True).mean()
        #aice = aice.coarsen(time=self.roll_win, boundary="trim").mean()
        attrs = ds.attrs.copy()
        attrs.update({"processed_by"     : "FastIceProcessor",
                      "source"           : "NSIDC CDR",
                      "processing_notes" : f"SIA/SIE calculated with threshold {self.icon_thresh} and rolling window {self.roll_win}" })
        ds_out = xr.Dataset({'SIA': (('t_nsidc',), SIA.values/self.SIC_scale),
                             'SIE': (('t_nsidc',), SIE.values/self.SIC_scale)},
                             # 'SIC': (('time', y_dim, x_dim), aice.data)},
                            coords = { 't_nsidc' : SIA['time'].values },
                                       # 'lat'  : ds['lat'],
                                       # 'lon'  : ds['lon']},
                            attrs = attrs )
        self.logger.info(f"✅ Computed NSIDC SIA/SIE in {time.time() - t1:.2f} seconds")
        return ds_out

    def load_AF2020_FI_area_CSV(self, doy_start):
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

    def regrid_AF2020_FI_to_CICE(self, FI_obs_native):
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
        lon_reG = regridder(FI_OBS_GRD["lon"]).persist()
        lat_reG = regridder(FI_OBS_GRD["lat"]).persist()
        FI_reG  = regridder(FI_obs_native["Fast_Ice_Time_series"]).persist()
        return FI_reG, lon_reG, lat_reG

    def load_AF2020_FI_org_netcdf(self, P_orgs):
        FI_obs = xr.open_mfdataset(P_orgs,
                                   combine='nested',
                                   concat_dim='time',
                                   parallel=True,
                                   chunks='auto',
                                   engine='netcdf4')
        FI_reG, lon_reG, lat_reG = self.regrid_AF2020_FI_to_CICE(FI_obs)
        mask     = xr.where(FI_reG >= 4, 1, np.nan)
        FI       = (('t_FI_obs', 'nj', 'ni'),
                    FI_reG.where(mask).values,
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
                          coords=dict(t_FI_obs=t_coords, obs_lon=x_coords, obs_lat=y_coords)).persist()

    def filter_AF2020_FI_by_date(self, start_date, end_date):
        D_obs    = Path(self.config['sea_ice_dict']['D_AF2020_db_org'])
        yrs_reqd = set([start_date.year, end_date.year])
        P_orgs   = [D_obs / f"FastIce_70_{yr}.nc" for yr in yrs_reqd]
        ds       = self.load_AF2020_FI_org_netcdf(P_orgs)
        # Convert FI_t_alt (int like 20000101) to datetime
        alt_dates = pd.to_datetime(ds['FI_t_alt'].values.astype(str), format='%Y%m%d')
        # Use where logic to find which obs period includes start_date
        valid_idx = (alt_dates >= pd.to_datetime(start_date)) & (alt_dates <= pd.to_datetime(end_date))
        matched = ds.sel(t_FI_obs=valid_idx).persist()
        if matched.dims['t_FI_obs'] == 0:
            self.logger.warning(f"No matching observational dates found between {start_date} and {end_date}")
        return matched

    def create_AF2020_FI_zarr(self):
        P_zarr = Path(self.config['sea_ice_dict']["P_AF_2020db_avg"])
        if P_zarr.exists() and not self.ow_AF2020_zarr:
            self.logger.info(f"Averaged observational gridded climatology already exists\n{P_zarr}")
            return xr.open_zarr(P_zarr)
        elif P_zarr.exists() and self.ow_AF2020_zarr:
            self.logger.info(f"Averaged observational gridded climatology exists *but* over-writing has been requested\n\tOVER-WRITING: {P_zarr}")
            shutil.rmtree(P_zarr)
        else:
            self.logger.info(f"Averaged observational gridded climatology does *NOT* exist\n\tCREATING NEW: {P_zarr}")
        # Load all gridded obs, already regridded and with date_alt available
        D_obs = Path(self.config['sea_ice_dict']['D_AF2020_db_org'])
        P_orgs = sorted(D_obs.glob("FastIce_70_*.nc"))
        self.logger.info(f"loading all observational fast ice netcdf files:\n{P_orgs}")
        ds_all = self.load_AF2020_FI_org_netcdf(P_orgs)
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
        grouped = grouped.assign_attrs(doy_vals=self.doy_vals).persist()
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
        if self.assoc_AF2020:
            window_end = window_start + pd.Timedelta(days=self.roll_win - 1)
            self.logger.info(f"loading data for date period: {window_start} to {window_end}")
            # Load observational climatology for the matching DOY window
            doy_start = int(window_start.strftime('%j'))
            self.logger.info(f"loading CSV climatology")
            cli_df = self.load_AF2020_FI_area_CSV(doy_start)
            sectors = cli_df.columns.tolist()
            # Load or compute gridded obs if within range, else use climatology
            if pd.Timestamp('2000-03-01') <= window_start <= pd.Timestamp('2018-03-31'):
                self.logger.info("using gridded observations")
                obs_gridded = self.filter_AF2020_FI_by_date(window_start, window_end)
                obs_grd_data = obs_gridded['FI'].mean(dim='t_FI_obs', skipna=True)
            else:
                self.logger.info("using gridded climatology")
                clim_all = self.create_AF2020_FI_zarr()
                closest_doy = min(self.doy_vals, key=lambda x: abs(x - doy_start))
                self.logger.info(f"using {closest_doy}-DOY from gridded climatology to associate with this model period")
                doy_index = np.argmin(np.abs(np.array(self.doy_vals) - closest_doy))
                #obs_grd_data = clim_all['FI_OBS_GRD'].sel(t_doy=closest_doy)
                obs_grd_data = clim_all['FI_OBS_GRD'].isel(t_doy=doy_index)
        if self.sea_ice or self.pack_ice:
            # Load NSIDC observational data
            self.logger.info(f"🧭 Including NSIDC SIA/SIE into model dataset for {self.sim_name}")
            D_NSIDC_orgs = Path(self.config["NSIDC_dict"]["D_original"], self.hemisphere, "daily")
            F_vers_dict = self.config["NSIDC_dict"]["file_versions"]
            F_vers_sorted = sorted(F_vers_dict.items(), key=lambda x: datetime.strptime(x[1], "%Y-%m-%d") if x[1] else datetime.min)
            P_NSIDC_orgs = []
            for d in self.dt_period_list:
                dt_str_nohyph = d.strftime('%Y%m%d')
                fver = next((ver for ver, date_str in reversed(F_vers_sorted)
                             if date_str and d >= datetime.strptime(date_str, "%Y-%m-%d")), F_vers_sorted[0][0])
                filename = f"seaice_conc_daily_{self.hemisphere_abbreviation.lower()}_{dt_str_nohyph}_{fver}_v04r00.nc"
                P_NSIDC_orgs.append(D_NSIDC_orgs / filename)
            self.logger.debug(f"Loading NSIDC files: {P_NSIDC_orgs}")
            t1 = time.time()
            # Manually load and concat files
            ds_list = []
            for f in sorted(P_NSIDC_orgs):
                ds = xr.open_dataset(f)
                ds_list.append(ds)
            NSIDC = xr.concat(ds_list, dim=self.config["NSIDC_dict"]["time_dim"])
            self.logger.info(f"✅ NSIDC dataset loaded: shape {NSIDC.sizes}, time: {time.time() - t1:.2f} s")
            NSIDC = self._convert_NSIDC_cartesian_to_spherical(NSIDC)
            NSIDC_out = self._compute_NSIDC_SIA_SIE(NSIDC)
            self.logger.info(f"✅ NSIDC processing complete, time: {time.time() - t1:.2f} s")
            #self._NSIDC_SIC = NSIDC_out['SIC']
            self._NSIDC_SIA = NSIDC_out['SIA']
            self._NSIDC_SIE = NSIDC_out['SIE']
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
        CICE = xr.open_mfdataset(P_CICE_orgs, combine='by_coords', chunks=None) #, parallel=True, preprocess=preprocess)
        self.logger.info(f"\✅ CICE files loaded: shape {CICE.sizes}, time: {time.time() - t1:.2f} s")
        if self.assoc_AF2020:
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
        regridder = xe.Regridder(self.gi_processor.G_u, self.gi_processor.G_t,
                                 method        = "bilinear",
                                 extrap_method = "inverse_dist",
                                 periodic      = True,
                                 weights       = self.CICE_dict['P_reG_u2t_weights'])
        self.logger.info(f"\ttime taken: {time.time()-t1:0.2f} seconds")
        u_reG = regridder(ds["uvel"])
        v_reG = regridder(ds["vvel"])
        return u_reG, v_reG

    def dataset_to_dictionary(self, ds):
        self.logger.info("***        DATASET TO DICTIONARY      ***")
        self.logger.info(f"   re-ogranise, roll and coarsen data")
        if self.sea_ice:
            json_varout_meta = self.config["SI_var_dict"]
        elif self.pack_ice:
            json_varout_meta = self.config["PI_var_dict"]
        else:
            json_varout_meta = self.config["FI_var_dict"]
        CICE_dict_unrolled = {}
        self.logger.info(f"create a dictionary of CICE variables to be used later")
        for v in self.cice_vars_reqd:
            if v in ds:
                CICE_dict_unrolled[v] = ds[v]
        uice, vice = self.CICE_regrid_to_tgrid(ds)
        CICE_dict_unrolled['speed'] = np.sqrt(uice ** 2 + vice ** 2)
        if self.sea_ice:
            self.logger.info(f"mapping CICE variable names to SI output-dataset dictionary found in the JSON configuration file")
        elif self.pack_ice:
            self.logger.info(f"mapping CICE variable names to PI output-dataset dictionary found in the JSON configuration file")
        else:
            self.logger.info(f"mapping CICE variable names to FI output-dataset dictionary found in the JSON configuration file")
        for out_var in json_varout_meta:
            meta = json_varout_meta.get(out_var, {})
            cice_var = meta.get("CICE_variable")
            vec_vars = meta.get("CICE_vec_vars", [])
            if vec_vars and cice_var:
                if all(vv in ds for vv in vec_vars):
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
        t1 = time.time()
        roll = lambda da: da.rolling(time=self.roll_win, center=True, min_periods=1).mean()
        CICE_dict_rolled = {k: roll(v) for k, v in CICE_dict_unrolled.items()}
        self.logger.info(f"time taken {time.time()-t1:0.2f} seconds")
        self.logger.info("coarsen the data arrays")
        t1 = time.time()
        coarse = lambda da: da.coarsen(time=self.roll_win, boundary="trim").mean()
        CICE_dict_coarsened = {k: coarse(v) for k, v in CICE_dict_rolled.items()}
        self.logger.info(f"time taken {time.time()-t1:0.2f} seconds")
        if self.assoc_AF2020:
            CICE_dict_coarsened['FI_OBS_CLI'] = ds['FI_OBS_CLI'].compute()
            CICE_dict_coarsened['FI_OBS_GRD'] = ds['FI_OBS_GRD'].compute()
        self.logger.info("persist the coarsened data in memory")
        CICE_dict_out = {}
        for k in CICE_dict_coarsened.keys():
            CICE_dict_out[k] = CICE_dict_coarsened[k].persist()
        return CICE_dict_out

    def apply_masks(self, roll_dict):
        self.logger.info("***     APPLYING MASKS     ***")
        t1 = time.time()
        sic_mask = (roll_dict['aice'] > self.icon_thresh)
        spd_mask = (roll_dict['speed'] <= self.ispd_thresh)
        if self.sea_ice:
            masked_dict    = {k: v.where(sic_mask) for k, v in roll_dict.items()}
            self.spat_mask = sic_mask.isel(nj=self.hemisphere_nj_slice).compute()
        elif self.pack_ice:
            pi_mask        = sic_mask & ~spd_mask
            masked_dict    = {k: v.where(pi_mask) for k, v in roll_dict.items()}
            self.spat_mask = pi_mask.isel(nj=self.hemisphere_nj_slice).compute()
        else:
            fi_mask        = sic_mask & spd_mask
            masked_dict    = { k: v.where(fi_mask) if k not in {'FI_OBS_CLI', 'FI_OBS_GRD'} else v for k, v in roll_dict.items() }
            self.spat_mask = fi_mask.isel(nj=self.hemisphere_nj_slice).compute()
        self.logger.info(f"time taken {time.time()-t1:0.2f} seconds")
        self.logger.info("persist the masked data in memory")
        masked_out = {}
        for k in masked_dict.keys():
            masked_out[k] = masked_dict[k].persist()
        return masked_out

    def compute_sea_ice_outputs(self, cice_vars_dict):
        def stringify_datetime(v):
            return v.isoformat() if isinstance(v, (datetime, pd.Timestamp)) else v
        if self.sea_ice:
            self.logger.info("\n ***     COMPUTING SEA ICE OUTPUTS    ***")
            json_varout_meta = self.config["SI_var_dict"]
            area_scale       = self.SIC_scale
        elif self.pack_ice:
            self.logger.info("\n ***     COMPUTING PACK ICE OUTPUTS    ***")
            json_varout_meta = self.config["PI_var_dict"]
            area_scale       = self.SIC_scale
        else:
            self.logger.info("\n ***     COMPUTING FAST ICE OUTPUTS    ***")
            json_varout_meta = self.config["FI_var_dict"]
            area_scale       = self.FIC_scale
        self.logger.info("STEP 1: SLICE HEMISPHERE")
        cice_vars_dict = self.slice_hemisphere(cice_vars_dict)
        self.logger.info("STEP 2: METHOD DEFINITIONS")
        grid_cell_area = self.gi_processor.G_t['area'].isel(nj=self.hemisphere_nj_slice).values
        cm2m           = self.cm2m_fact
        one_d_list     = [k for k, meta in json_varout_meta.items() if meta.get("dimensions") == "1D"]
        two_d_list     = [k for k, meta in json_varout_meta.items() if meta.get("dimensions") == "2D"]
        three_d_list   = [k for k, meta in json_varout_meta.items() if meta.get("dimensions") == "3D"]
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
        if not self.sea_ice and not self.pack_ice:
            sector_coord_tup = ((sector_dim_str), cice_vars_dict['FI_OBS_CLI'].sector.values)
        #################################
        #######   1D VARIABLES   ########
        #################################
        self.logger.info("STEP 3: COMPUTE 1D OUTPUT VARIABLES")
        one_d_metrics = {}
        for v in one_d_list:
            if v.endswith("_SD") or (v in ['FIA_CLI','SIA','SIE','SIC']):
                continue
            meta     = json_varout_meta.get(v, {})
            cice_var = meta.get("CICE_variable", None)
            self.logger.debug(f"\t📦 Available keys in cice_vars_dict: {list(cice_vars_dict.keys())}")
            if not cice_var or cice_var not in cice_vars_dict.keys():
                self.logger.warning(f"\t⚠️ Skipping 1D metric {v} — source variable '{cice_var}' missing.")
                continue
            else:
                self.logger.debug(f"\t✅ Creating 1D metric {v} — source variable '{cice_var}'")
            dat = cice_vars_dict[cice_var]
            if v in ["FIA", "PIA", "SIA"] :
                one_d_metrics[v] = ((dat * grid_cell_area).sum(dim=two_dim_tup)) / area_scale
            elif v in ["FIE", "PIE", "SIE"]:
                one_d_metrics[v] = ((self.spat_mask * grid_cell_area).sum(dim=two_dim_tup)) / area_scale
            elif "AGING" in v:
                one_d_metrics[v] = (grid_cell_area / dat).sum(dim=two_dim_tup)
            elif "VGRO" in v or "FRAZIL" in v:
                one_d_metrics[v] = (dat * cm2m * grid_cell_area).sum(dim=two_dim_tup)
            else:
                one_d_metrics[v] = (dat * grid_cell_area).sum(dim=two_dim_tup)
        one_d_vars = {k: xr.DataArray(data   = v,
                                      dims   = one_dim_tup,
                                      coords = {t_dim_str : t_coord_tup},
                                      attrs  = json_varout_meta.get(k, {}))
                      for k, v in one_d_metrics.items()}
        if self.assoc_AF2020:
            one_d_vars['FIA_OBS'] = xr.DataArray(data   = cice_vars_dict['FI_OBS_CLI']/1e3,
                                                 dims   = sector_dim_tup,
                                                 coords = {t_dim_str      : t_coord_tup,
                                                           sector_dim_str : sector_coord_tup},
                                                 attrs  = json_varout_meta.get('FIA_OBS', {}))
        if self.sea_ice or self.pack_ice:
            one_d_vars['SIA_OBS'] = xr.DataArray(data   = getattr(self, f"_NSIDC_SIA", None),
                                                 dims   = one_dim_tup,
                                                 coords = {t_dim_str : t_coord_tup},
                                                 attrs  = json_varout_meta.get('SIA_OBS', {}))
            one_d_vars['SIE_OBS'] = xr.DataArray(data   = getattr(self, f"_NSIDC_SIE", None),
                                                 dims   = one_dim_tup,
                                                 coords = {t_dim_str : t_coord_tup},
                                                 attrs  = json_varout_meta.get('SIE_OBS', {}))
        self.logger.debug(f"1D vars computed: {list(one_d_vars.keys())}")
        #################################
        #######   3D VARIABLES   ########
        #################################
        self.logger.info("STEP 3:\n\tCOMPUTE 3D OUTPUT VARIABLES:")
        self.logger.info("\t very little processing done and essentially CICE variables are only filtered/masked for fast ice criteria")
        t1 = time.time()
        three_d_vars = {}
        for v in three_d_list:
            if v.endswith("_SD") or (v in ['FI_GRD','SIC']):
                continue
            meta     = json_varout_meta.get(v, {})
            cice_var = meta.get("CICE_variable", None)
            self.logger.debug(f"\t📦 Available keys in cice_vars_dict: {list(cice_vars_dict.keys())}")
            if not cice_var or cice_var not in cice_vars_dict.keys():
                self.logger.warning(f"\t⚠️ Skipping 3D metric {v} — source variable '{cice_var}' missing.")
                continue
            else:
                self.logger.debug(f"\t✅ Creating 3D metric {v} — source variable '{cice_var}'")
                data = cice_vars_dict[cice_var]
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
        if self.assoc_AF2020:
            FI_OBS                 = cice_vars_dict['FI_OBS_GRD'].isel(nj=self.hemisphere_nj_slice)
            FI_OBS                 = FI_OBS.expand_dims({t_dim_str:[self.dt0_period]})
            FI_OBS                 = FI_OBS.drop_vars(['TLON','TLAT','ULON','ULAT','obs_lon','obs_lat'])
            three_d_vars['FI_OBS'] = xr.DataArray(data   = FI_OBS,
                                                  dims   = three_dim_tup,
                                                  coords = {t_dim_str   : t_coord_tup,
                                                            x_coord_str : x_coord_tup,
                                                            y_coord_str : y_coord_tup},
                                                  attrs  = json_varout_meta.get('FI_OBS', {}))
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
            if v in ['FIP_OBS','SIP']:
                continue
            meta     = json_varout_meta.get(v, {})
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
            two_d_vars[v] = xr.DataArray(data   = data_mean,
                                         dims   = two_dim_tup,
                                         coords = {x_coord_str : x_coord_tup,
                                                   y_coord_str : y_coord_tup},
                                         attrs  = {**meta})
                                                   # "start_time": stringify_datetime(self.dt0_period),
                                                   # "stop_time" : stringify_datetime(self.dtN_period)})
        if self.assoc_AF2020:
            FIP_OBS               = cice_vars_dict['FI_OBS_GRD'].isel(nj=self.hemisphere_nj_slice)
            FIP_OBS               = FIP_OBS.drop_vars(['TLON','TLAT','ULON','ULAT','obs_lon','obs_lat'])
            two_d_vars['FIP_OBS'] = xr.DataArray(data   = FIP_OBS,
                                                 dims   = two_dim_tup,
                                                 coords = {x_coord_str : x_coord_tup,
                                                           y_coord_str : y_coord_tup},
                                                 attrs  = {**json_varout_meta.get('FIP_OBS', {})})
                                                           # "start_time": stringify_datetime(self.dt0_period),
                                                           # "stop_time" : stringify_datetime(self.dtN_period)})
        self.logger.info(f"\ttime taken: {time.time()-t1:0.2f} seconds")
        self.logger.debug(f"\t2D vars computed: {list(three_d_vars.keys())}")
        #######################################
        #######      OUTPUT DATASET     #######
        #######################################
        self.logger.info("STEP 6:\n\tCREATE OUTPUT DATASET")
        t1 = time.time()
        OUT = xr.Dataset({**one_d_vars, **three_d_vars, **two_d_vars}).persist()
        # ENSURE 'time' dimensions is *completely* scrubbed from dataset
        OUT = OUT.map(lambda da: da.squeeze(dim='time') if 'time' in da.dims else da)
        OUT = OUT.map(lambda da: da.drop_dims('time') if 'time' in da.dims else da)
        if 'time' in OUT.coords:
            OUT = OUT.drop_vars('time')
        for var in OUT.data_vars:
            if 'time' in OUT[var].dims:
                OUT[var] = OUT[var].squeeze(dim='time', drop=True)
        assert 'time' not in OUT.dims
        # don't forget about the attributes
        if self.sea_ice:
            title                      = f"Sea ice analysed by SeaIceProcessor, simulation name: {self.sim_name}"
            ice_concentration_criteria = f"sea ice concentration per grid cell > {self.icon_thresh}"
            ice_speed_criteria         = f"no masking for speed"
        elif self.pack_ice:
            title                      = f"Pack ice analysed by SeaIceProcessor, simulation name: {self.sim_name}"
            ice_concentration_criteria = f"sea ice concentration per grid cell > {self.icon_thresh}"
            ice_speed_criteria         = f"sea ice speed per grid cell > {self.ispd_thresh} m/s"
        else:
            title                      = f"Fast ice analysed by SeaIceProcessor, simulation name: {self.sim_name}"
            ice_concentration_criteria = f"sea ice concentration per grid cell > {self.icon_thresh}"
            ice_speed_criteria         = f"sea ice speed per grid cell <= {self.ispd_thresh} m/s"
        OUT.attrs = {"title"                      : title,
                     "summary"                    : "This dataset includes sea ice variables and derived metrics using a "\
                                                    "mean rolling window and temporal coarsening methods, then a data "\
                                                    "masking method based on threshold provided below.",
                     "source"                     : "CICE v6.4.1 model output",
                     "creator_name"               : "Daniel Patrick Atwater",
                     "creator_email"              : "daniel.atwater@utas.edu.au",
                     "institution"                : "Institute of Marine and Antarctic Studies--University of Tasmania",
                     "history"                    : f"Created on {datetime.now().isoformat()}",
                     "references"                 : "",
                     "conventions"                : "CF-1.8",
                     "ice_concentration_criteria" : ice_concentration_criteria,
                     "ice_speed_criteria"         : ice_speed_criteria,
                    "grounded_iceberg_db"         : self.gi_processor.GI_dataset_path,
                    "landmask_file"               : self.gi_processor.KMT_path,
                    "total_area_GI"               : self.GI_total_area,
                    "time_coverage_start"         : stringify_datetime(self.dt0_period),
                    "time_coverage_end"           : stringify_datetime(self.dtN_period),
                    "roll_window_days"            : self.roll_win,
                    "geospatial_lat_min"          : float(np.min(self.gi_processor.G_t['lat'].isel(nj=self.hemisphere_nj_slice).values)),
                    "geospatial_lat_max"          : float(np.max(self.gi_processor.G_t['lat'].isel(nj=self.hemisphere_nj_slice).values)),
                    "geospatial_lon_min"          : float(np.min(self.gi_processor.G_t['lon'].isel(nj=self.hemisphere_nj_slice).values)),
                    "geospatial_lon_max"          : float(np.max(self.gi_processor.G_t['lon'].isel(nj=self.hemisphere_nj_slice).values))}
        self.logger.info(f"\ttime taken: {time.time()-t1:0.2f} seconds")
        return OUT

    def process_window(self, dt0_str=None, dtN_str=None, rolling_window=None, write_zarr=False, ow_zarrs=False):
        """
        Default behaviour is to assume user is interactively creating output datasets
        for a particular simulation and time period. Allows a user to test rolling
        mean period and how this effects the dataset. However, this *TURNS OFF*
        any and all loading of the fast ice observational climatology and observational gridded
        datasets. This is because the fast ice observational datasets are somewhat 'set in stone',
        so to speak. The observational datasets are based on 15-day periods starting on
        1 January to the 346th year-day and then uses a 20-day period. At present
        this rigidity in the observational temporal dimension puts constraints on
        dynamically associating CICE model data with it.

        The above averaging is used consistently across fast ice outputs, pack ice outputs and
        sea ice outputs as the default. 

        However, for now if one uses process_window() in the following manner
        dt0_str  = "1993-09-01"
        dtN_str  = "1994-12-31"
        sim_name = 'ktens-max'
        SI_proc  = SeaIceProcessor(sim_name = sim_name)
        FI       = SI_proc.process_window(dt0_str = dt0_str,
                                          dtN_str = dtN_str)
        Then one can expect to get FI output dataset with the default rolling window 15 days up to the last 20-day
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
        def is_leap(date_object):
            date_object.is_leap_year if hasattr(date_object, 'is_leap_year') else (date_object.year % 4 == 0 and
                                                                                   (date_object.year % 100 != 0 or
                                                                                    date_object.year % 400 == 0))
        #dask.config.set(scheduler='single-threaded')
        self.dt0_str      = dt0_str if dt0_str is not None else self.config.get("dt0_str", "1993-01-01")
        self.dtN_str      = dtN_str if dtN_str is not None else self.config.get("dtN_str", "1993-12-31")
        self.dt0          = pd.Timestamp(self.dt0_str)
        self.dtN          = pd.Timestamp(self.dtN_str)
        if (rolling_window and rolling_window != 15) or self.sea_ice or self.pack_ice:
            self.assoc_AF2020 = False
        dt0_list = []
        for year in range(self.dt0.year, self.dtN.year + 1):
            for doy in self.doy_vals:
                dt = pd.Timestamp(datetime(year, 1, 1)) + pd.Timedelta(days=doy - 1)
                if self.dt0 <= dt <= self.dtN:
                    dt0_list.append(dt)
        DS_CAT = []
        for self.dt0_period in dt0_list:
            doy                 = self.dt0_period.timetuple().tm_yday
            last_doy            = 366 if is_leap(self.dt0_period) else 365
            roll_win            = 20 if doy >= self.doy_vals[-1] else 15
            self.roll_win       = roll_win if rolling_window is None else rolling_window
            self.dtN_period     = self.dt0_period + pd.Timedelta(days=self.roll_win)
            self.dt_period_list = pd.date_range(self.dt0_period, self.dtN_period, freq='1D')
            ds                  = self.load_data_window(self.dt0_period)
            ds_dict             = self.dataset_to_dictionary(ds)
            masked_vars         = self.apply_masks(ds_dict)
            DS                  = self.compute_sea_ice_outputs(masked_vars)
            if write_zarr:
                F_zarr_out = self.F_zarr_out_fmt.format(date_str=self.dt0_period.strftime('%Y-%m-%d'))
                P_zarr_out = Path(self.D_zarr_out, F_zarr_out)
                if not P_zarr_out.exists() or ow_zarrs:
                    self.logger.info(f"*** writing OUTPUT DATASET to disk: {P_zarr_out}")
                    t1 = time.time()
                    DS.to_zarr(P_zarr_out, mode='w')
                    self.logger.info(f"time taken: {time.time()-t1:0.2f} seconds")
                else:
                    self.logger.info(f"*** skipping write to {P_zarr_out}")
                    self.logger.info("OUTPUT zarr file already exists and overwriting disabled")
            DS_CAT.append(DS)
        DS_RETURN = xr.concat(DS_CAT, dim='t_dim').compute()
        self.logger.info("✅ sea ice processing complete.")
        return DS_RETURN
