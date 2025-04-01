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
    def __init__(self, sim_name, dt0_str=None, dtN_str=None, mask_fast_ice=False, compute_rolling_mean=True,
                 hemisphere=None, P_log=None, json_path=None):
        """
        Initialize the PackIceProcessor.
        """
        self.sim_name = sim_name
        self.mask_fi  = mask_fast_ice
        self.use_gi   = False  # GI excluded for pack ice

        if json_path is None:
            json_path = "/home/581/da1339/AFIM/src/AFIM/src/JSONs/afim_cice_analysis.json"
        with open(json_path, 'r') as f:
            self.config = json.load(f)

        self.dt0_str  = dt0_str if dt0_str is not None else self.config.get("dt0_str", "1993-01-01")
        self.dtN_str  = dtN_str if dtN_str is not None else self.config.get("dtN_str", "1999-12-31")
        self.dt_range = pd.date_range(self.dt0_str, self.dtN_str)
        self.dt0      = self.dt_range[0]
        self.dtN      = self.dt_range[-1]

        # Set rolling window behavior
        if self.mask_fi:
            self.compute_rolling_mean = True
        else:
            self.compute_rolling_mean = compute_rolling_mean

        self.roll_win = self.config.get("roll_win", 15) if self.compute_rolling_mean else 1

        # Center dates for rolling window processing
        if self.compute_rolling_mean:
            self.dtC_list = pd.date_range(
                self.dt0 + pd.Timedelta(days=self.roll_win // 2),
                self.dtN - pd.Timedelta(days=self.roll_win // 2),
                freq=f"{self.roll_win}D"
            )
        else:
            self.dtC_list = [self.dt0 + pd.Timedelta(days=len(self.dt_range) // 2)]

        # First center date for logging/debugging
        self.dtC     = self.dtC_list[0]
        self.dtC_str = self.dtC.strftime("%Y-%m-%d")

        if P_log is None:
            P_log = Path(self.config['D_dict']['logs'], f"PackIceProcessor_{sim_name}.log")
        self.setup_logging(logfile=P_log)

        # Simulation and config metadata
        self.sim_config           = self.config['sim_dict'][sim_name]
        self.sim_dir              = Path(self.config['D_dict']['AFIM_out'], sim_name, "history", "daily")
        self.gi_processor         = GroundedIcebergProcessor(self.config, sim_name)
        self.gi_processor.load_grid_and_landmask()
        hemisphere                = hemisphere if hemisphere is not None else self.config.get('hemisphere', 'south')
        self.define_hemisphere(hemisphere)
        self.FI_thresh            = self.config.get("FI_thresh", 5e-4)
        self.SIC_thresh           = self.config.get("SIC_thresh", 0.15)
        self.SIC_scale            = self.config.get("SIC_scale", 1e12)
        self.cm2m_fact            = self.config.get("cm2m_fact", 0.01)
        self.CICE_dict            = self.config['CICE_dict']
        self.regrid_weights_path = self.CICE_dict['P_reG_u2t_weights']


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
        if hemisphere.lower() in ['north', 'northern', 'nh', 'n', 'no']:
            self.hemisphere_nj_slice = slice(540, 1080)
            self.hemisphere_abbreviation = 'nh'
            self.hemisphere = 'north'
        elif hemisphere.lower() in ['south', 'southern', 'sh', 's', 'so']:
            self.hemisphere_nj_slice = slice(0, 540)
            self.hemisphere_abbreviation = 'sh'
            self.hemisphere = 'south'
        else:
            raise ValueError(f"Invalid hemisphere: {hemisphere}")

    def slice_hemisphere(self, var_dict):
        return {k: v.isel(nj=self.hemisphere_nj_slice) for k, v in var_dict.items()}

    def _extract_cice_vars_from_PI_var_dict(self, var_list):
        pi_meta = self.config["PI_var_dict"]
        var_set = set()
        for v in var_list:
            meta = pi_meta.get(v, {})
            cice_var = meta.get("CICE_variable")
            vec_vars = meta.get("CICE_vector_variables", [])
            # Skip derived output variables that aren't actually in the model dataset
            if cice_var and cice_var not in ['speed', 'strint', 'strair', 'strocn', 'strtlt', 'strcor']:
                var_set.add(cice_var)
                # Always add vector components
            if isinstance(vec_vars, list):
                var_set.update(vec_vars)
        # Always include aice (used in masking and core metrics)
        var_set.add("aice")
        self.logger.debug(f"🧾 CICE variables required for computation: {sorted(var_set)}")
        return sorted(var_set)

    def load_data_window(self, process_NSIDC=None):
        if not hasattr(self, 'var_list'):
            raise AttributeError("`self.var_list` must be defined before calling load_data_window().")
        self.process_NSIDC = process_NSIDC if process_NSIDC is not None else self.config["NSIDC_dict"].get("process_SIA", False)
        P_CICE_orgs = [
            self.sim_dir / f"iceh.{d.strftime('%Y-%m-%d')}.nc"
            for d in self.dt_range if (self.sim_dir / f"iceh.{d.strftime('%Y-%m-%d')}.nc").exists()
        ]
        if not P_CICE_orgs:
            raise FileNotFoundError(f"No CICE files found for window around {self.dtC}")
        self.logger.debug(f"Loading model files: {P_CICE_orgs}")
        t1 = time.time()
        preprocess = lambda ds: ds[list(self.var_list)]
        CICE = xr.open_mfdataset(P_CICE_orgs, combine='by_coords', parallel=True,
                                 preprocess=preprocess)#, chunks=self.chunk_dict)
        self.logger.info(f"✅ Model dataset loaded: shape {CICE.sizes}, time: {time.time() - t1:.2f} s")
        if not self.process_NSIDC:
            return CICE
        # Load NSIDC observational data
        self.logger.info(f"🧭 Including NSIDC SIA/SIE into model dataset for {self.sim_name}")
        D_NSIDC_orgs = Path(self.config["NSIDC_dict"]["D_original"], self.hemisphere, "daily")
        F_vers_dict = self.config["NSIDC_dict"]["file_versions"]
        F_vers_sorted = sorted(F_vers_dict.items(), key=lambda x: datetime.strptime(x[1], "%Y-%m-%d") if x[1] else datetime.min)
        P_NSIDC_orgs = []
        for d in self.dt_range:
            dt_str_nohyph = d.strftime('%Y%m%d')
            fver = next((ver for ver, date_str in reversed(F_vers_sorted)
                         if date_str and d >= datetime.strptime(date_str, "%Y-%m-%d")), F_vers_sorted[0][0])
            filename = f"seaice_conc_daily_{self.hemisphere_abbreviation}_{dt_str_nohyph}_{fver}_v04r00.nc"
            P_NSIDC_orgs.append(D_NSIDC_orgs / filename)
        self.logger.debug(f"Loading NSIDC files: {P_NSIDC_orgs}")
        t1 = time.time()
        # Manually load and concat files
        ds_list = []
        for f in sorted(P_NSIDC_orgs):
            ds = xr.open_dataset(f)
            ds_list.append(ds)
        NSIDC = xr.concat(ds_list, dim=self.config["NSIDC_dict"]["time_dim"])
        #NSIDC = xr.open_mfdataset(P_NSIDC_orgs, parallel=True)
        self.logger.info(f"✅ NSIDC dataset loaded: shape {NSIDC.sizes}, time: {time.time() - t1:.2f} s")
        NSIDC = self._convert_NSIDC_cartesian_to_spherical(NSIDC)
        NSIDC_out = self._compute_NSIDC_SIA_SIE(NSIDC)
        self.logger.info(f"✅ NSIDC processing complete, time: {time.time() - t1:.2f} s")
        self._NSIDC_SIC = NSIDC_out['SIC']  # Store separately so it’s not merged
        self._NSIDC_SIA = NSIDC_out['SIA']
        self._NSIDC_SIE = NSIDC_out['SIE']
        return CICE

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
        self.logger.info(f"\t✅ Conversion complete in {time.time()-t1:.2f} seconds")
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
        mask = aice > self.SIC_thresh
        SIA = (aice * area).where(mask.notnull()).sum(dim=[y_dim, x_dim], skipna=True)
        SIE = (mask * area).sum(dim=[y_dim, x_dim], skipna=True)
        aice = aice.where(mask)
        if self.compute_rolling_mean:
            SIA = SIA.rolling(time=self.roll_win, center=True).mean()
            SIE = SIE.rolling(time=self.roll_win, center=True).mean()
            aice = aice.rolling(time=self.roll_win, center=True).mean()
        attrs = ds.attrs.copy()
        if self.compute_rolling_mean:
            proc_notes = f"SIA/SIE calculated with threshold {self.SIC_thresh} and rolling window {self.roll_win}"
        else:
            proc_notes = f"SIA/SIE calculated with threshold {self.SIC_thresh}"
        attrs.update({
            "processed_by": "PackIceProcessor",
            "source": "NSIDC CDR",
            "processing_notes": proc_notes
        })
        ds_out = xr.Dataset(
            {
                'SIA': (('time',), SIA.data/self.SIC_scale),
                'SIE': (('time',), SIE.data/self.SIC_scale),
                'SIC': (('time', y_dim, x_dim), aice.data)
            },
            coords={
                'time': SIA['time'],
                'lat': ds['lat'],
                'lon': ds['lon']
            },
            attrs=attrs
        )
        self.logger.info(f"✅ Computed SIA/SIE in {time.time() - t1:.2f} seconds")
        return ds_out

    def regrid_to_tgrid(self, ds):
        self.logger.info("Regridding uvel and vvel to T-grid...")
        with suppress_stderr():
            regridder = xe.Regridder(self.gi_processor.G_u, self.gi_processor.G_t,
                                     method="bilinear", extrap_method="inverse_dist", periodic=True,
                                     weights=self.regrid_weights_path)
        U = regridder(ds['uvel'])
        V = regridder(ds['vvel'])
        #return self.slice_hemisphere(U), self.slice_hemisphere(V)
        return U,V

    def apply_PI_mask(self, var_dict):
        self.logger.info("applying mask")
        if self.mask_fi:
            fi_mask = (var_dict['aice'] > self.SIC_thresh) & (var_dict['speed'] <= self.FI_thresh)
            pi_mask = (var_dict['aice'] > self.SIC_thresh) & (~fi_mask)
        else:
            pi_mask = (var_dict['aice'] > self.SIC_thresh)
        # need to save a sliced mask for use in PIE calculation later
        self.pi_mask = pi_mask.isel(nj=self.hemisphere_nj_slice)
        return {k: v.where(pi_mask) for k, v in var_dict.items()}

    def dataset_to_dictionary(self, ds, var_list=None):
        pi_meta = self.config["PI_var_dict"]
        full_dict = {}

        # Always include base vars
        base_vars = ['aice']
        if self.process_NSIDC:
            base_vars += ['_NSIDC_SIC', '_NSIDC_SIE', '_NSIDC_SIA']
        for v in base_vars:
            if v in ds:
                full_dict[v] = ds[v]

        # Compute speed if masking fast ice
        if self.mask_fi:
            U, V = self.regrid_to_tgrid(ds)
            full_dict['speed'] = np.sqrt(U**2 + V**2)

        # Add variables needed for var_list
        for out_var in var_list:
            meta = pi_meta.get(out_var, {})
            cice_var = meta.get("CICE_variable")
            vec_vars = meta.get("CICE_vector_variables", [])

            if vec_vars and cice_var:
                if all(vv in ds for vv in vec_vars):
                    # Derived variable from vector components
                    full_dict[cice_var] = np.sqrt(sum([ds[vv] ** 2 for vv in vec_vars]))
                else:
                    missing = [vv for vv in vec_vars if vv not in ds]
                    self.logger.debug(
                        f"⚠️ Skipping vector-derived var {out_var} — missing components: {missing}"
                    )

            elif cice_var and cice_var in ds:
                # Scalar variable available directly in ds
                full_dict[cice_var] = ds[cice_var]

        self.logger.debug(f"🧾 dataset_to_dictionary keys: {list(full_dict.keys())}")

        # Apply rolling mean if enabled
        if self.compute_rolling_mean:
            self.logger.info("compute rolling mean")
            roll = lambda da: da.rolling(time=self.roll_win, center=True, min_periods=1).mean()
            return {k: roll(v) for k, v in full_dict.items()}
        else:
            return full_dict

    def compute_pack_ice_outputs(self, mask_vars, var_list=None):
        grid_cell_area = self.gi_processor.G_t['area'].isel(nj=self.hemisphere_nj_slice)
        spatial_dims   = self.CICE_dict['spatial_dims']
        time_dim       = self.CICE_dict['time_dim']
        time_coords    = mask_vars['aice'].time.values
        lon_coords     = self.gi_processor.G_t['lon'].isel(nj=self.hemisphere_nj_slice).values
        lat_coords     = self.gi_processor.G_t['lat'].isel(nj=self.hemisphere_nj_slice).values
        cm2m           = self.cm2m_fact
        roll_win       = self.roll_win if self.compute_rolling_mean else 1
        if self.mask_fi:
            lon_coord_name = self.CICE_dict['FI_lon_coord']
            lat_coord_name = self.CICE_dict['FI_lat_coord']
            FI_time_dim    = self.CICE_dict['FI_time_dim']
            three_dims     = self.CICE_dict['FI_three_dims']
        else:
            lon_coord_name = self.CICE_dict['FI_lon_coord']
            lat_coord_name = self.CICE_dict['FI_lat_coord']
            FI_time_dim    = self.CICE_dict['time_dim']
            three_dims     = self.CICE_dict['three_dims']
        mask_vars_hem = self.slice_hemisphere(mask_vars)
        pi_meta = self.config["PI_var_dict"]
        # Categorize var_list by dimension
        one_d_list = [v for v in var_list if pi_meta.get(v, {}).get("dimensions") == "1D"]
        three_d_list = [v for v in var_list if pi_meta.get(v, {}).get("dimensions") == "3D"]
        two_d_list = [v for v in var_list if pi_meta.get(v, {}).get("dimensions") == "2D"]
        #########################
        ### 1D Variables
        #########################
        self.logger.info("compute pack ice 1D variables:")
        one_d_metrics = {}
        for v in one_d_list:
            if v.endswith("_SD") or (v in ['SIA','SIE','SIC']):
                continue
            meta = pi_meta.get(v, {})
            cice_var = meta.get("CICE_variable", None)
            self.logger.debug(f"📦 Available keys in mask_vars_hem: {list(mask_vars_hem.keys())}")
            if not cice_var or cice_var not in mask_vars_hem.keys():
                self.logger.warning(f"⚠️ Skipping 1D metric {v} — source variable '{cice_var}' missing.")
                continue
            else:
                self.logger.debug(f"✅ Creating 1D metric {v} — source variable '{cice_var}'")
            if v == "PIA":
                one_d_metrics[v] = ((mask_vars_hem[cice_var] * grid_cell_area).sum(dim=spatial_dims)) / self.SIC_scale
            elif v == "PIE" and hasattr(self, "pi_mask"):
                one_d_metrics[v] = ((self.pi_mask * grid_cell_area).sum(dim=spatial_dims)) / self.SIC_scale
            elif "AGING" in v:
                one_d_metrics[v] = (grid_cell_area / mask_vars_hem[cice_var]).sum(dim=spatial_dims)
            elif "VGRO" in v or "FRAZIL" in v:
                one_d_metrics[v] = (mask_vars_hem[cice_var] * cm2m * grid_cell_area).sum(dim=spatial_dims)
            else:
                one_d_metrics[v] = (mask_vars_hem[cice_var] * grid_cell_area).sum(dim=spatial_dims)
        one_d_vars = {k: xr.DataArray(data=v.data,
                                      dims=(time_dim,),
                                      coords={time_dim: time_coords},
                                      attrs=pi_meta.get(k, {}))
                      for k, v in one_d_metrics.items()}
        self.logger.info(f"\t1D vars computed: {list(one_d_vars.keys())}")
        #########################
        ### 3D Variables
        #########################
        self.logger.info("compute pack ice 3D variables:")
        three_d_vars = {}
        for v in three_d_list:
            if v.endswith("_SD") or (v in ['SIA','SIE','SIC']):
                continue
            meta = pi_meta.get(v, {})
            cice_var = meta.get("CICE_variable", None)
            self.logger.debug(f"📦 Available keys in mask_vars_hem: {list(mask_vars_hem.keys())}")
            if not cice_var or cice_var not in mask_vars_hem.keys():
                self.logger.warning(f"⚠️ Skipping 3D metric {v} — source variable '{cice_var}' missing.")
                continue
            else:
                self.logger.debug(f"✅ Creating 3D metric {v} — source variable '{cice_var}'")
            try:
                if self.compute_rolling_mean:
                    data = mask_vars_hem[cice_var].coarsen(time=roll_win, boundary="trim").mean()
                    time_coords_3d = data.time.values
                    time_dim = FI_time_dim
                else:
                    data = mask_vars_hem[cice_var]
                    time_coords_3d = time_coords
                if "VGRO" in v or "FZL" in v:
                    data = data * cm2m
                three_d_vars[v] = xr.DataArray(data=data.data,
                                               dims=three_dims,
                                               coords={time_dim: time_coords_3d,
                                                       lon_coord_name: (spatial_dims, lon_coords),
                                                       lat_coord_name: (spatial_dims, lat_coords)},
                                               attrs=meta)
            except Exception as e:
                self.logger.warning(f"⚠️ Skipping 3D var {v} due to error: {e}")
        self.logger.info(f"\t3D vars computed: {list(three_d_vars.keys())}")
        #########################
        ### 2D Variables
        #########################
        self.logger.info("compute temporal means to give spatial distributions over time (2D)")
        two_d_vars = {}
        for v in two_d_list:
            if not v.endswith("_SD"):
                continue
            if v in ['SIA','SIE','SIC']:
                continue
            base = v.replace("_SD", "")
            meta = pi_meta.get(v, {})
            cice_var = pi_meta.get(base, {}).get("CICE_variable", None)
            self.logger.debug(f"📦 Available keys in mask_vars_hem: {list(mask_vars_hem.keys())}")
            if not cice_var or cice_var not in mask_vars_hem.keys():
                self.logger.warning(f"⚠️ Skipping 2D var {v} due to missing base variable '{cice_var}'")
                continue
            else:
                self.logger.debug(f"✅ Creating 2D metric {v} — source variable '{cice_var}'")
            try:
                da = mask_vars_hem[cice_var]
                # Determine time dimension automatically
                time_dim_auto = next((d for d in da.dims if d in ['time', 't', 't_fi']), None)
                if not time_dim_auto:
                    raise ValueError(f"No valid time dimension found for 2D variable {v}")
                norm = da.sizes[time_dim_auto]
                data_sum = da.sum(dim=time_dim_auto)
                if base in ["PIHI", "PISTH", "PISH", "PIDIV", "PIAG", "PIAD", "PIAT",
                            "PIVD", "PIVT", "PISTR", "PISPD", "PIFZL"]:
                    max_val = da.max().values
                    data_mean = data_sum / (norm * max_val)
                else:
                    data_mean = data_sum / norm
                two_d_vars[v] = xr.DataArray(data=data_mean.data,
                                             dims=spatial_dims,
                                             coords={lon_coord_name: (spatial_dims, lon_coords),
                                                     lat_coord_name: (spatial_dims, lat_coords)},
                                             attrs={**meta,
                                                    "start_time": self.dt0_str.isoformat() if isinstance(self.dt0_str, datetime) else self.dt0_str,
                                                    "stop_time": self.dtN_str.isoformat() if isinstance(self.dtN_str, datetime) else self.dtN_str})
            except Exception as e:
                self.logger.warning(f"⚠️ Skipping 2D var {v} due to error: {e}")
        self.logger.info(f"\t2D vars computed: {list(two_d_vars.keys())}")
        #########################
        ### Final Dataset
        #########################
        self.logger.info("create output dataset:")
        PI = xr.Dataset({**one_d_vars, **three_d_vars, **two_d_vars})
        if self.process_NSIDC:
            for var in ['SIA', 'SIE', 'SIC']:
                da = getattr(self, f"_NSIDC_{var}", None)
                if da is not None:
                    meta = pi_meta.get(var, {})
                    PI[var] = xr.DataArray(data=da.data, dims=da.dims, coords=da.coords, attrs=meta)
        fi_criteria = (f"aice > {self.SIC_thresh} and speed <= {self.FI_thresh} m/s"
                       if self.mask_fi else "fast ice masking not employed on this dataset")
        PI.attrs = {"title": "Pack ice analysed from numerical sea ice model simulations",
                    "summary": "Pack ice metrics computed from model output, masked by fast ice criteria where requested.",
                    "source": "CICE v6.4.1 model output",
                    "creator_name": "Daniel Patrick Atwater",
                    "creator_email": "daniel.atwater@utas.edu.au",
                    "institution": "Institute of Marine and Antarctic Studies--University of Tasmania",
                    "history": f"Created on {datetime.now().isoformat()}",
                    "conventions": "CF-1.8",
                    "fast_ice_criteria": fi_criteria,
                    "pack_ice_criteria": f"aice > {self.SIC_thresh} and not fast_ice_criteria",
                    "landmask_file": self.gi_processor.KMT_path,
                    "time_coverage_start": self.dt0_str.isoformat() if isinstance(self.dt0_str, datetime) else self.dt0_str,
                    "time_coverage_end": self.dtN_str.isoformat() if isinstance(self.dtN_str, datetime) else self.dtN_str,
                    "roll_window_days": roll_win,
                    "geospatial_lat_min": float(np.min(lat_coords)),
                    "geospatial_lat_max": float(np.max(lat_coords)),
                    "geospatial_lon_min": float(np.min(lon_coords)),
                    "geospatial_lon_max": float(np.max(lon_coords))}
        return PI

    def process_window(self, var_list=None, save_zarr=False, ow_zarrs=False):
        import dask
        dask.config.set(scheduler='single-threaded')
        pi_meta = self.config["PI_var_dict"]
        valid_keys = set(pi_meta.keys())
        # ⬇️ Force full computation if saving with fast ice masking
        if save_zarr and self.mask_fi:
            var_list = list(valid_keys)
        elif var_list is None:
            var_list = list(valid_keys) if self.mask_fi else ['PIA', 'PIE']
        else:
            var_list = [v for v in var_list if v in valid_keys]
        self.requested_PI_vars = var_list  # ✅ Used by load_data_window()
        if not var_list:
            self.logger.warning("⚠️ No valid output variables selected — nothing will be computed.")
        else:
            self.logger.debug(f"✅ Pack ice variables selected for output: {var_list}")
        # ✅ Extract required CICE variables for loading
        self.var_list = self._extract_cice_vars_from_PI_var_dict(var_list)
        if self.mask_fi and save_zarr:
            ds_all = []
            for dtC in self.dtC_list:
                self.dtC = dtC
                self.dt_range = pd.date_range(
                    dtC - pd.Timedelta(days=self.roll_win // 2),
                    dtC + pd.Timedelta(days=self.roll_win // 2)
                )
                ds = self.load_data_window()
                self.logger.debug(f"📦 Variables in loaded dataset: {list(ds.data_vars)}")
                ds_dict = self.dataset_to_dictionary(ds, self.requested_PI_vars)
                masked_vars = self.apply_PI_mask(ds_dict)
                PI = self.compute_pack_ice_outputs(masked_vars, var_list=self.requested_PI_vars)
                # Output path
                D_PI_zarr = Path(self.config['D_dict']['AFIM_out'], self.sim_name, "PI")
                F_PI_zarr = f"pack_ice_{dtC.strftime('%Y-%m-%d')}.zarr"
                P_PI_zarr = D_PI_zarr / F_PI_zarr
                if not D_PI_zarr.exists():
                    os.makedirs(D_PI_zarr)

                if not P_PI_zarr.exists() or ow_zarrs:
                    self.logger.info(f"*** writing PI to disk: {P_PI_zarr}")
                    PI.to_zarr(P_PI_zarr, mode='w')
                else:
                    self.logger.info(f"Zarr exists or overwriting disabled: {P_PI_zarr}")
                ds_all.append(PI)
            PI_merged = xr.concat(ds_all, dim='time')
            self.logger.info("✅ Pack ice processing complete.")
            return PI_merged
        else:
            # Interactive mode with first center date
            self.dtC = self.dtC_list[0]
            self.dt_range = pd.date_range(
                self.dtC - pd.Timedelta(days=self.roll_win // 2),
                self.dtC + pd.Timedelta(days=self.roll_win // 2)
            )
            ds = self.load_data_window()
            ds_dict = self.dataset_to_dictionary(ds, self.requested_PI_vars)
            masked_vars = self.apply_PI_mask(ds_dict)
            PI = self.compute_pack_ice_outputs(masked_vars, var_list=self.requested_PI_vars)
            self.logger.info("✅ Pack ice processing complete.")
            return PI
