import json, os, shutil, time, logging
import xarray as xr
import pandas as pd
import numpy  as np
import xesmf  as xe
from pathlib import Path
from datetime import datetime, timedelta

class SeaIceObservations:
    def __init__(self, **kwargs):
        return
    
    ####################################################################################################
    ##                                              NSIDC         
    ####################################################################################################
    def load_local_NSIDC(self, dt0_str=None, dtN_str=None, local_directory=None, monthly_files=False):
        """
        Load NSIDC sea ice concentration data from local NetCDF files over a date range.

        INPUTS:
           dt0_str         : str, optional; start date in 'YYYY-MM-DD' format. Defaults to self.dt0_str.
           dtN_str         : str, optional; end date in 'YYYY-MM-DD' format. Defaults to self.dtN_str.
           local_directory : str or Path, optional; directory containing the NSIDC NetCDF files.
                             If None, uses self.NSIDC_dict["D_original"].
           monthly_files   : bool, optional; iff True, loads monthly files; otherwise loads daily files.

        OUTPUTS
           xarray.Dataset; combined dataset loaded with xarray.open_mfdataset.
        """
        def promote_time(ds):
            ds = ds[["cdr_seaice_conc"]] if "cdr_seaice_conc" in ds else ds
            if "time" in ds.variables and "tdim" in ds.dims:
                ds = ds.swap_dims({"tdim": "time"})
                ds = ds.set_coords("time")
            return ds
        f_fmt    = self.NSIDC_dict["G02202_v4_file_format"]
        dt0_str  = dt0_str if dt0_str is not None else self.dt0_str
        dtN_str  = dtN_str if dtN_str is not None else self.dtN_str
        freq_str = 'monthly'             if monthly_files  else 'daily'
        D_local  = Path(local_directory) if local_directory else Path(self.NSIDC_dict["D_original"], self.hemisphere, freq_str)
        dt0      = datetime.strptime(dt0_str, '%Y-%m-%d')
        dtN      = datetime.strptime(dtN_str, '%Y-%m-%d')
        f_vers_parsed = [ (ver, datetime.strptime(date_str, "%Y-%m-%d") if date_str else datetime.min)
                          for ver, date_str in self.NSIDC_dict["file_versions"].items() ]
        f_vers_parsed.sort(key=lambda x: x[1])  # Sort by date
        date_range = pd.date_range(start=dt0, end=dtN, freq='MS' if monthly_files else 'D')
        P_NSIDCs   = []
        for dt_ in date_range:
            f_version = next( (ver for ver, d in reversed(f_vers_parsed) if dt_ >= d), f_vers_parsed[0][0] )
            F_        = f_fmt.format(freq = freq_str,
                                     hem  = self.hemisphere_abbreviation.lower(),
                                     date = dt_.strftime('%Y%m%d'),
                                     ver  = f_version)
            P_        = Path(D_local,F_)
            if P_.exists():
                P_NSIDCs.append(P_)
            else:
                self.logger.warning(f"Missing file: {P_}")
        if not P_NSIDCs:
            raise FileNotFoundError(f"No NSIDC files found in {D_local} between {dt0_str} and {dtN_str}")
        self.logger.info(f"Using xarray to load {len(P_NSIDCs)} files in {D_local}. Last file: {P_NSIDCs[-1].name}")
        return xr.open_mfdataset(P_NSIDCs,
                                 combine    = "by_coords",
                                 parallel   = True,
                                 preprocess = promote_time)

    def NSIDC_coordinate_transformation(self, ds):
        from pyproj import CRS, Transformer
        """
        Convert NSIDC dataset Cartesian coordinates (xgrid, ygrid) into geographic coordinates (longitude, latitude).

        This method transforms the Cartesian coordinates used in NSIDC datasets into
        standard geographic coordinates (longitude, latitude) based on the Polar
        Stereographic projection.

        INPUTS
            ds (xarray.Dataset): The NSIDC dataset containing Cartesian coordinates.

        OUTPUTS
            xarray.Dataset: The dataset with added `lon` and `lat` variables representing
                            geographic coordinates.

        NOTES:
            - The conversion is done using the Proj4 string specified in `self.NSIDC_dict["projection_string"]`.
            - The transformed geographic coordinates are added to the dataset as new variables
              named `lon` and `lat`.
            - The method assumes the dataset contains `xgrid` and `ygrid` variables, which
              represent the Cartesian coordinates in meters.
        """
        x = ds['xgrid'].values
        y = ds['ygrid'].values
        assert x.shape == y.shape, "xgrid and ygrid must have the same shape"
        crs_nsidc = CRS.from_proj4(self.NSIDC_dict["projection_string"])
        crs_wgs84 = CRS.from_epsg(self.sea_ice_dic["projection_wgs84"])
        trans     = Transformer.from_crs(crs_proj, crs_wgs84, always_xy=True)
        lon, lat  = transformer.transform(x, y)
        ds['lat'] = (('y', 'x'), lat)
        ds['lon'] = (('y', 'x'), lon)
        return ds

    def compute_NSIDC_metrics(self,
                              dt0_str         = None,
                              dtN_str         = None,
                              local_load_dir  = None,
                              monthly_files   = False,
                              P_zarr          = None,
                              overwrite       = False):
        flags    = self.NSIDC_dict["cdr_seaice_conc_flags"]
        SIC_name = self.NSIDC_dict["SIC_name"]
        dt0_str  = dt0_str if dt0_str is not None else self.dt0_str
        dtN_str  = dtN_str if dtN_str is not None else self.dtN_str
        if monthly_files:
            freq_str = 'monthly'
        else:
            freq_str = 'daily'
        D_local = local_load_dir if local_load_dir is not None else Path(self.NSIDC_dict["D_original"], self.hemisphere, freq_str)
        P_zarr  = P_zarr         if P_zarr         is not None else Path(D_local, "zarr" )
        if not os.path.exists(P_zarr) or overwrite:
            NSIDC = self.load_local_NSIDC(dt0_str         = dt0_str,
                                          dtN_str         = dtN_str,
                                          local_directory = D_local)
            aice  = NSIDC[SIC_name]
            for flag in flags:
                aice = xr.where(aice==flag/100, np.nan, aice)
            area     = xr.open_dataset(self.NSIDC_dict["P_cell_area"]).cell_area / self.SIC_scale
            mask     = aice > self.icon_thresh
            SIA      = (aice * area).where(mask.notnull()).sum(dim=['y', 'x'], skipna=True)
            SIE      = (mask * area).sum(dim=['y', 'x'], skipna=True)
            NSIDC_ts = xr.Dataset({'SIA' : (('time'),SIA.values),
                                   'SIE' : (('time'),SIE.values)},
                                  coords = {'time' : (('time'), SIA.time.values)})
            self.logger.info(f"**writing** NSIDC time series zarr file {P_zarr}")
            NSIDC_ts.to_zarr(P_zarr)
            return NSIDC_ts
        else:
            self.logger.info(f"**loading** previously created NSIDC time series zarr file {P_zarr}")
            return xr.open_zarr(P_zarr, consolidated=True)

    #-------------------------------------------------------------------------------------------
    #                                 OBSERVATIONAL DATA
    #-------------------------------------------------------------------------------------------
    def load_AF2020_FI_area_timeseries(self, year: int = None) -> pd.DataFrame:
        """
        Load the full AF2020 observed landfast sea ice area timeseries CSV
        and optionally filter by a specific year.

        Parameters
        ----------
        year : int, optional
            If provided, filters the CSV to return only rows matching this year.

        Returns
        -------
        df : pd.DataFrame
            A DataFrame containing Circumpolar and sector ice area columns,
            indexed by (year, DOY_start).
        """
        csv_path = self.AF_FI_dict['P_AF2020_csv']
        df       = pd.read_csv(csv_path)
        # Standardize column names
        df = df.rename(columns={
            'Circumpolar': 'circumpolar',
            'IOsector': 'IOsector',
            'WPOsector': 'WPOsector',
            'RSsector': 'RSsector',
            'BASsector': 'BASsector',
            'WSsector': 'WSsector'
        })
        # Add multi-index and optional filter
        df = df.set_index(['Year', 'DOY_start'])
        if year is not None:
            if year not in df.index.get_level_values('Year'):
                raise ValueError(f"Year {year} not found in AF2020 CSV.")
            df = df.loc[year]
        self.logger.info(f"Loaded AF2020 FI timeseries from: {csv_path}")
        return df.reset_index()

    def load_AF2020_FI_area_CSV(self, doy_start):
        csv_path    = self.AF_FI_dict['P_AF2020_cli_csv']
        df          = pd.read_csv(csv_path)
        closest_doy = min(self.AF_FI_dict['DOY_vals'], key=lambda x: abs(x - doy_start))
        row         = df[df['DOY_start'] == closest_doy].copy()
        row         = row.rename(columns={'Circumpolar': 'circumpolar',
                                          'IOsector'   : 'IOsector',
                                          'WPOsector'  : 'WPOsector',
                                          'RSsector'   : 'RSsector',
                                          'BASsector'  : 'BASsector',
                                          'WSsector'   : 'WSsector'})
        sectors     = ['Circumpolar', 'IOsector', 'WPOsector', 'RSsector', 'BASsector', 'WSsector']
        self.logger.info(f"loaded: {csv_path}")
        return row[sectors].reset_index(drop=True)

    def interpolate_obs_fia(self, csv_df, full_doy=np.arange(1, 366)):
        from scipy.interpolate import interp1d
        df = csv_df.copy()
        if 'Circumpolar' not in df.columns or 'DOY_start' not in df.columns:
            raise ValueError("Expected columns 'DOY_start' and 'circumpolar' in CSV.")
        x = df['DOY_start'].values
        y = df['Circumpolar'].values/1e3
        if len(x) != len(y):
            raise ValueError("x and y arrays must be equal in length.")
        interp_func = interp1d(x, y, kind='linear', fill_value="extrapolate")
        return xr.DataArray(interp_func(full_doy), coords=[("doy", full_doy)])

    def load_AF2020db(self, dt0_str="2000-01-01", dtN_str="2018-12-31", overwrite_zarr=False):
        t0 = time.time()
        D_zarr = Path(self.Sea_Ice_Obs_dict['D_local'],"FI_obs")
        D_zarr.mkdir(parents=True, exist_ok=True)
        P_zarr = Path(D_zarr,"AF-FI-2020db_org-timestep_reG.zarr")
        if Path(P_zarr).exists() and any(Path(P_zarr).iterdir()) and not overwrite_zarr:
            self.logger.info(f"Loading existing AF2020 observational Zarr file: {P_zarr}")
            ds = xr.open_zarr(P_zarr, consolidated=Path(P_zarr, ".zmetadata").exists())
        else:
            self.logger.info(f"{P_zarr} not found *or* overwrite requested. Generating from NetCDF...")
            ds = self._load_AF2020_FI_org_and_save_zarr(P_zarr)
        self.logger.info(f"Zarr load/generation completed in {time.time()-t0:.2f} seconds")
        dt0 = datetime.strptime(dt0_str, "%Y-%m-%d")
        dtN = datetime.strptime(dtN_str, "%Y-%m-%d")
        ds = ds.assign_coords(obs_date=("t_FI_obs", pd.to_datetime(ds['FI_t_alt'].values.astype(str), format="%Y%m%d")))
        ds = ds.set_index(t_FI_obs="obs_date")
        ds = ds.sel(t_FI_obs=slice(dt0, dtN))
        if ds.dims.get("t_FI_obs", 0) == 0:
            self.logger.warning(f"No matching observational dates found between {dt0_str} and {dtN_str}")
        return ds

    def _load_AF2020_FI_org_and_save_zarr(self, P_zarr):
        from pyproj import CRS, Transformer
        D_obs = Path(self.AF_FI_dict['D_AF2020_db_org'])
        P_orgs = sorted(D_obs.glob("FastIce_70_*.nc"))
        self.logger.info(f"Loading original AF2020 NetCDF files: {P_orgs}")
        FI_obs = xr.open_mfdataset(P_orgs, engine='netcdf4', combine='by_coords')
        FI_obs = FI_obs.chunk({'time': 32}) 
        FI_OBS = xr.where(FI_obs['Fast_Ice_Time_series'] >= 4, 1.0, np.nan)
        FI_OBS.name = "FI"
        FI_OBS.attrs["long_name"] = FI_obs['Fast_Ice_Time_series'].attrs['long_name']
        crs_obs = CRS.from_epsg(self.AF_FI_dict["projection_FI_obs"])
        crs_sph = CRS.from_epsg(self.AF_FI_dict["projection_wgs84"])
        transformer = Transformer.from_crs(crs_obs, crs_sph, always_xy=True)
        X, Y = np.meshgrid(FI_obs['x'].isel(time=0).values, FI_obs['y'].isel(time=0).values)
        lon_obs, lat_obs = transformer.transform(X, Y)
        FI_OBS = FI_OBS.assign_coords(lon=("nj", "ni", lon_obs), lat=("nj", "ni", lat_obs))
        G_t = self.define_cice_grid(grid_type='t', mask=False, build_grid_corners=False)
        P_weights = self.AF_FI_dict["AF_reG_weights"]
        reuse_weights = os.path.exists(P_weights)
        self.logger.info(f"Regridding observational data using weights: {P_weights} (reuse={reuse_weights})")
        reG_obs = xe.Regridder(FI_OBS.isel(time=0), G_t,
                            method="bilinear", periodic=False,
                            ignore_degenerate=True,
                            reuse_weights=reuse_weights,
                            filename=P_weights)
        FI_OBS_reG = reG_obs(FI_OBS)
        FI_OBS_reG.name = "FI"
        FI_OBS_reG.attrs["long_name"] = FI_obs['Fast_Ice_Time_series'].attrs['long_name']
        # Construct Dataset
        ds_out = xr.Dataset(data_vars={
            "FI": FI_OBS_reG,
            "FI_t_alt": ("t_FI_obs", FI_obs.date_alt,
                        {"long_name": FI_obs.date_alt.attrs.get("long_name", ""),
                        "description": FI_obs.date_alt.attrs.get("description", "")})
        }, coords={
            "t_FI_obs": ("t_FI_obs", FI_obs.time,
                        {"description": "Start date of 15- or 20-day image mosaic window.",
                        "units": "days since 2000-1-1 0:0:0"}),
            "lon": (("nj", "ni"), G_t['lon'].values, {"units": "degrees_east"}),
            "lat": (("nj", "ni"), G_t['lat'].values, {"units": "degrees_north"})
        })
        self.logger.info(f"Saving regridded dataset to Zarr: {P_zarr}")
        ds_out.to_zarr(P_zarr, mode="w", consolidated=True)
        return ds_out

    # def define_AF2020_reG_weights(self, FI_obs_da):
    #     from pyproj import CRS, Transformer
    #     self.logger.info("define model grid")
    #     G_t = self.define_cice_grid( grid_type='t' , mask=False , build_grid_corners=False )
    #     G_t.to_netcdf("/g/data/gv90/da1339/grids/CICE_Tgrid_for_AF2020db_regrid.nc")
    #     self.logger.info("defining AF2020 regridder weights")
    #     F_weights     = self.AF_FI_dict["AF_reG_weights"]
    #     weights_exist = os.path.exists(F_weights)
    #     self.logger.info("convert 'AF_FI_OBS_2020db' Cartesian coordinates to spherical coordindates")
    #     crs_obs          = CRS.from_epsg(self.AF_FI_dict["projection_FI_obs"]) #unique to observations
    #     crs_spherical    = CRS.from_epsg(self.AF_FI_dict["projection_wgs84"])  #spherical coordinates
    #     transformer      = Transformer.from_crs(crs_obs, crs_spherical, always_xy=True)
    #     X, Y             = np.meshgrid(FI_obs_native['x'].isel(time=0).values, FI_obs_native['y'].isel(time=0).values)
    #     lon_obs, lat_obs = transformer.transform(X,Y)
    #     da_obs           = xr.DataArray(data   = FI_obs_da,
    #                                     dims   = ["nj", "ni"],
    #                                     coords = dict(lon = (["nj", "ni"], lon_obs),
    #                                                   lat = (["nj", "ni"], lat_obs)))
    #     self.logger.info(f"Model lon: {G_t['lon'].values.min()} to {G_t['lon'].values.max()}")
    #     self.logger.info(f"Obs lon:   {G_obs['lon'].min()} to {G_obs['lon'].max()}")
    #     self.logger.info(f"{'🔁 Reusing' if weights_exist else '⚙️ Creating'} regrid weights: {F_weights}")
    #     self.reG_AF2020 = xe.Regridder(G_obs, G_t,
    #                                    method            = "bilinear",
    #                                    periodic          = False,
    #                                    ignore_degenerate = True,
    #                                    reuse_weights     = weights_exist,
    #                                    filename          = F_weights)

    # def load_AF2020_FI_org_netcdf(self, P_orgs):
    #     from pyproj import CRS, Transformer
    #     import xesmf as xe
    #     # Setup
    #     P_weights     = self.AF_FI_dict["AF_reG_weights"]
    #     reuse_weights = os.path.exists(P_weights)
    #     self.logger.info("loading FI observations with xarray mfdataset")
    #     self.logger.info(f"loading these files:\n{P_orgs}")
    #     FI_obs = xr.open_mfdataset(P_orgs, engine='netcdf4', chunks={'time': 1})
    #     self.logger.info("masking FI observations for values greater than 4")
    #     FI_OBS                    = xr.where(FI_obs['Fast_Ice_Time_series'] >= 4, 1.0, np.nan)
    #     FI_OBS.name               = "FI"
    #     FI_OBS.attrs["long_name"] = FI_obs['Fast_Ice_Time_series'].attrs['long_name']
    #     self.logger.info("define model grid")
    #     G_t = self.define_cice_grid(grid_type='t', mask=False, build_grid_corners=False)
    #     self.logger.info("convert 'AF_FI_OBS_2020db' Cartesian coordinates to spherical coordindates")
    #     crs_obs          = CRS.from_epsg(self.AF_FI_dict["projection_FI_obs"])
    #     crs_sph          = CRS.from_epsg(self.AF_FI_dict["projection_wgs84"])
    #     transformer      = Transformer.from_crs(crs_obs, crs_sph, always_xy=True)
    #     X, Y             = np.meshgrid(FI_obs['x'].isel(time=0).values, FI_obs['y'].isel(time=0).values)
    #     lon_obs, lat_obs = transformer.transform(X, Y)
    #     FI_OBS = FI_OBS.assign_coords(lon=(("nj", "ni"), lon_obs),
    #                                   lat=(("nj", "ni"), lat_obs))
    #     self.logger.info("*** Regridding 'AF_FI_OBS_2020db' to CICE T-grid ***")
    #     self.logger.info(f"\tDefine regridder; reusing weights: {reuse_weights}")
    #     reG_obs = xe.Regridder(FI_OBS.isel(time=0),  # representative grid for weights
    #                            G_t,
    #                            method="bilinear",
    #                            periodic=False,
    #                            ignore_degenerate=True,
    #                            reuse_weights=reuse_weights,
    #                          filename=P_weights)
    #     self.logger.info("\tRegridding masked observational dataset")
    #     FI_OBS_reG = reG_obs(FI_OBS)  # shape: (time, nj, ni)
    #     FI_OBS_reG.name = "FI"
    #     FI_OBS_reG.attrs["long_name"] = FI_obs['Fast_Ice_Time_series'].attrs['long_name']
    #     ds_out = xr.Dataset(data_vars = {"FI"       : FI_OBS_reG,
    #                                      "FI_t_alt" : (("t_FI_obs",), FI_obs.date_alt, {"long_name"   : FI_obs.date_alt.attrs.get("long_name", ""),
    #                                                                                     "description" : FI_obs.date_alt.attrs.get("description", "")})},
    #                         coords    = {"t_FI_obs" : (("t_FI_obs",), FI_obs.time, {"description": "Start date of 15- or 20-day image mosaic window.",
    #                                                                                 "units": "days since 2000-1-1 0:0:0"}),
    #                                      "lon"      : (("nj", "ni"), G_t['lon'].values, {"units": "degrees_east"}),
    #                                      "lat"      : (("nj", "ni"), G_t['lat'].values, {"units": "degrees_north"})})
    #     self.logger.info("converted AF2020 database for use with SeaIceProcessor (Dask-friendly)")
    #     return ds_out

    # def filter_AF2020_FI_by_date(self, dt0_str=None, dtN_str=None):
    #     dt0_str   = dt0_str or self.dt0_str
    #     dtN_str   = dtN_str or self.dtN_str 
    #     dt0       = datetime.strptime(dt0_str, "%Y-%m-%d")
    #     dtN       = datetime.strptime(dtN_str, "%Y-%m-%d")
    #     D_obs     = Path(self.AF_FI_dict['D_AF2020_db_org'])
    #     yrs_reqd  = list(range(dt0.year, dtN.year + 1))
    #     P_orgs    = [D_obs / f"FastIce_70_{yr}.nc" for yr in yrs_reqd]
    #     ds        = self.load_AF2020_FI_org_netcdf(P_orgs)
    #     alt_dates = pd.to_datetime(ds['FI_t_alt'].values.astype(str), format='%Y%m%d')
    #     ds        = ds.assign_coords(obs_date=("t_FI_obs", alt_dates))
    #     ds        = ds.set_index(t_FI_obs="obs_date")
    #     matched   = ds.sel(t_FI_obs=slice(dt0, dtN))
    #     if matched.dims['t_FI_obs'] == 0:
    #         self.logger.warning(f"No matching observational dates found between {dt0_str} and {dtN_str}")
    #     return matched

    # def create_AF2020_FI_zarr(self, overwrite=False):
    #     P_zarr = Path(self.AF_FI_dict["P_AF_2020db_avg"])
    #     if P_zarr.exists():
    #         if overwrite:
    #             self.logger.info(f"Averaged observational gridded climatology exists *but* over-writing has been requested\n\tOVER-WRITING: {P_zarr}")
    #             shutil.rmtree(P_zarr)
    #         else:
    #             self.logger.info(f"Averaged observational gridded climatology already exists\n{P_zarr}")
    #             return xr.open_zarr(P_zarr)
    #     self.logger.info(f"Averaged observational gridded climatology does *NOT* exist\n\tCREATING NEW: {P_zarr}")
    #     D_obs  = Path(self.AF_FI_dict['D_AF2020_db_org'])
    #     P_orgs = sorted(D_obs.glob("FastIce_70_*.nc"))
    #     self.logger.info(f"loading all observational fast ice netcdf files:\n{P_orgs}")
    #     ds_all    = self.load_AF2020_FI_org_netcdf(P_orgs)
    #     alt_dates = pd.to_datetime(ds_all['FI_t_alt'].astype(str), format='%Y%m%d')
    #     doy_vals  = alt_dates.dayofyear
    #     ds_all    = ds_all.assign_coords(doy=("t_FI_obs", doy_vals))
    #     ds_all    = ds_all.where(ds_all.doy.isin(self.AF_FI_dict['DOY_vals']), drop=True)
    #     grouped   = ds_all.groupby("doy").mean(dim="t_FI_obs", skipna=True)
    #     grouped   = grouped.rename_dims({"doy": "t_doy"})
    #     grouped   = grouped.rename_vars({"FI": "FI_OBS_GRD"})
    #     grouped   = grouped.assign_attrs(doy_vals=self.AF_FI_dict['DOY_vals']).persist()
    #     self.logger.info(f"gridded climatology now looks like\n{grouped}")
    #     self.logger.info(f"writing this to {P_zarr}")
    #     t1 = time.time()
    #     grouped.to_zarr(P_zarr)
    #     self.logger.info(f"\tzarr written in {time.time()-t1:0.2f} seconds")
    #     return grouped

