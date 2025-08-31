from __future__ import annotations
import json, os, shutil, time, logging, re
import xarray   as xr
import pandas   as pd
import numpy    as np
import xesmf    as xe
from pathlib    import Path
from datetime   import datetime, timedelta
from typing     import Iterable, List, Optional, Tuple

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
                                     hem  = self.hemisphere_dict['abbreviation'].lower(),
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

    ####################################################################################################
    ##                                              AF2020         
    ####################################################################################################
    def load_AF2020_FIA_summary(self, 
                                repeat_climatology: bool = True,
                                start             : str  = "1994-01-01",
                                end               : str  = "1999-12-31") -> xr.Dataset:
        """
        Load AF2020 landfast sea ice area time series and compute day-of-year climatology.

        Parameters
        ----------
        repeat_climatology : bool, optional
            If True, returns climatology repeated over the date range. If False, only returns 1D raw data.
        start : str, optional
            Start date for repeated climatology (default: "1994-01-01").
        end : str, optional
            End date for repeated climatology (default: "1999-12-31").
        smooth : int, optional
            Rolling window size for smoothing the day-of-year climatology (default: 5 days).

        Returns
        -------
        xr.Dataset
            Dataset with:
                - 'FIA_obs':           Daily FIA time series [time, region]
                - 'FIA_clim':          Smoothed climatology [doy, region]
                - 'FIA_clim_repeat':   Climatology repeated over [start, end] (optional)
        """
        from scipy.interpolate import interp1d
        csv_path = self.AF_FI_dict['P_AF2020_csv']
        df = pd.read_csv(csv_path)
        regions = ["circumpolar", "indian_ocean", "western_pacific", "ross_sea", "amundsen_sea", "weddell_sea"]
        df = df.rename(columns={'Circumpolar': 'circumpolar',
                                'IOsector'   : 'indian_ocean',
                                'WPOsector'  : 'western_pacific',
                                'RSsector'   : 'ross_sea',
                                'BASsector'  : 'amundsen_sea',
                                'WSsector'  : 'weddell_sea'})
        df = df.dropna(subset=["Year", "Month_start", "Day_of_month_start"])
        df["date"] = pd.to_datetime({"year": df["Year"].astype(int),
                                    "month": df["Month_start"].astype(int),
                                    "day": df["Day_of_month_start"].astype(int)}, errors="coerce")
        df = df.dropna(subset=["date"])
        df["doy"] = df["date"].dt.dayofyear
        df = df.set_index("date")
        df = df[regions + ["doy"]]
        fia_obs = df[regions] / 1000.0  # km² → 1000 km²
        fia_obs_xr = fia_obs.to_xarray().to_array("region").transpose("date", "region")
        if 'doy' not in df.columns or 'circumpolar' not in df.columns:
            raise ValueError("Expected 'doy' and 'circumpolar' columns for climatology generation.")
        doy_vals = np.arange(1, 366)
        x = df.groupby("doy")["circumpolar"].mean().dropna().index.values
        y = df.groupby("doy")["circumpolar"].mean().dropna().values / 1000.0
        interp_func = interp1d(x, y, kind='linear', fill_value="extrapolate")
        clim_interp = pd.DataFrame({"circumpolar": interp_func(doy_vals)}, index=doy_vals)
        fia_clim_doy = xr.DataArray(clim_interp["circumpolar"].values[:, None],  # <-- Add [:, None] to make it 2D
                                    coords={"doy": doy_vals, "region": ["circumpolar"]},
                                    dims=["doy", "region"])

        if repeat_climatology:
            model_dates = pd.date_range(start, end, freq="D")
            doy_map = model_dates.dayofyear
            valid_mask = doy_map < 366
            model_dates = model_dates[valid_mask]
            doy_map = doy_map[valid_mask]
            repeated = fia_clim_doy.sel(doy=xr.DataArray(doy_map, dims="time")).values
            fia_clim_re = xr.DataArray(repeated,
                                       coords = {"time": model_dates, "region": ["circumpolar"]},
                                       dims   = ["time", "region"])
        else:
            fia_clim_re = None
        ds = xr.Dataset({"FIA_obs": fia_obs_xr.sel(region="circumpolar"),
                        "FIA_clim": fia_clim_doy})
        if fia_clim_re is not None:
            ds["FIA_clim_repeat"] = fia_clim_re
        self.logger.info(f"AF2020 FIA data loaded from {csv_path}")
        return ds

    def AF2020_interpolate_FIA_clim(self, csv_df, full_doy=np.arange(1, 366)):
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

    def AF2020_clim_to_model_time(self, model: xr.DataArray, clim: xr.DataArray) -> xr.DataArray:
        """
        Expand a DOY-based climatology to match model time steps.

        Parameters
        ----------
        model : xr.DataArray
            Model time series with `time` coordinate.
        clim : xr.DataArray
            Climatology indexed by DOY and region.

        Returns
        -------
        xr.DataArray
            Climatology repeated over `model.time`, matched by DOY.
        """
        # Compute day-of-year for each model date (drop Feb 29 if not in obs)
        doy = xr.DataArray(model['time'].dt.dayofyear, coords={"time": model.time}, dims="time")
        # Match each model day to a climatological value
        matched = clim.sel(doy=doy, method="nearest")  # or method="pad"
        matched.coords["time"] = model.time  # Align with model time
        return matched

    def build_AF2020_grid_corners(self, lat_c, lon_c):
        """
        Construct (ny+1, nx+1) corner arrays from center lat/lon arrays.
        Assumes a regular 2D rectilinear grid in EPSG:3412 projection.

        Parameters:
        -----------
        lat_c, lon_c : 2D arrays (ny, nx)
            Latitude and longitude of pixel centers.

        Returns:
        --------
        lat_b, lon_b : 2D arrays (ny+1, nx+1)
            Latitude and longitude of grid cell corners.
        """
        ny, nx = lat_c.shape
        lat_padded = np.pad(lat_c, ((1, 1), (1, 1)), mode='edge')
        lon_padded = np.pad(lon_c, ((1, 1), (1, 1)), mode='edge')
        lat_b = 0.25 * (lat_padded[0:-1, 0:-1] + lat_padded[0:-1, 1:] +
                        lat_padded[1:, 0:-1] + lat_padded[1:, 1:])
        lon_b = 0.25 * (lon_padded[0:-1, 0:-1] + lon_padded[0:-1, 1:] +
                        lon_padded[1:, 0:-1] + lon_padded[1:, 1:])
        lon_b = ((lon_b + 360) % 360)
        return lat_b, lon_b

    def _load_AF2020_FI_org_and_save_zarr(self, P_zarr):
        D_obs  = Path(self.AF_FI_dict['D_AF2020_db_org'])
        P_orgs = sorted(D_obs.glob("FastIce_70_*.nc"))
        self.logger.info(f"Loading original AF2020 NetCDF files: {P_orgs}")
        FI_obs        = xr.open_mfdataset(P_orgs, engine='netcdf4', combine='by_coords')
        FI_obs        = FI_obs.chunk({'time': 32})
        lat_c         = FI_obs.latitude
        lon_c         = FI_obs.longitude
        lat_b, lon_b  = self.build_AF2020_grid_corners(lat_c.values, lon_c.values)
        G_obs         = xr.Dataset({"lat"   : (("Y", "X"), lat_c),
                                    "lon"   : (("Y", "X"), lon_c),
                                    "lat_b" : (("Y_b", "X_b"), lat_b),
                                    "lon_b" : (("Y_b", "X_b"), lon_b)})
        G_t           = self.define_cice_grid(grid_type='t', mask=False, build_grid_corners=True)
        P_weights     = self.AF_FI_dict["AF_reG_weights"]
        reuse_weights = os.path.exists(P_weights)
        self.logger.info(f"Regridding observational data using weights: {P_weights} (reuse={reuse_weights})")
        reG_obs = xe.Regridder(G_obs, G_t,
                               method            = "conservative",
                               periodic          = False,
                               ignore_degenerate = True,
                               reuse_weights     = reuse_weights,
                               filename          = P_weights)
        fast_ice_mask = xr.where(FI_obs['Fast_Ice_Time_series'] >= 4, 1.0, 0.0)
        FI_frac       = fast_ice_mask.sum(dim='time') / fast_ice_mask['time'].size
        FI_frac       = FI_frac.where(FI_frac > 0)  # set 0s to NaN (optional: mask ocean)
        FI_frac_reG   = reG_obs(FI_frac)
        ds_out = xr.Dataset(data_vars = {"FI"       : FI_frac_reG,
                                         "FI_t_alt" : ("t_FI_obs", FI_obs.date_alt, {"long_name"   : FI_obs.date_alt.attrs.get("long_name", ""),
                                                                                     "description" : FI_obs.date_alt.attrs.get("description", "")})},
                            coords = {"t_FI_obs" : ("t_FI_obs", FI_obs.time, {"description": "Start date of 15- or 20-day image mosaic window.",
                                                                              "units": "days since 2000-1-1 0:0:0"}),
                                      "lon"      : (("nj", "ni"), G_t['lon'].values, {"units": "degrees_east"}),
                                      "lat"      : (("nj", "ni"), G_t['lat'].values, {"units": "degrees_north"})})
        self.logger.info(f"Saving regridded dataset to Zarr: {P_zarr}")
        ds_out.to_zarr(P_zarr, mode="w", consolidated=True)
        return ds_out

    ####################################################################################################
    ##                                          CCI SIT         
    ####################################################################################################
    def _norm_list(self, x: Optional[Iterable[str]]) -> Optional[List[str]]:
        if x is None: return None
        out = [str(s).strip() for s in x if s and str(s).strip()]
        return out or None

    def _days_in_month(self, y: int, m: int) -> int:
        if m in (1,3,5,7,8,10,12): return 31
        if m in (4,6,9,11): return 30
        return 29 if ((y % 4 == 0 and y % 100 != 0) or (y % 400 == 0)) else 28

    def _month_overlap(self, y: int, m: int, t0: pd.Timestamp, t1: pd.Timestamp) -> bool:
        first = pd.Timestamp(year=y, month=m, day=1, tz="UTC")
        last  = pd.Timestamp(year=y, month=m, day=self._days_in_month(y,m), tz="UTC")
        return not (last < t0 or first > t1)

    def _parse_yyyymmdd(self, name: str) -> Optional[pd.Timestamp]:
        m = re.search(r"(?<!\d)(20\d{2})(0[1-9]|1[0-2])(0[1-9]|[12]\d|3[01])(?!\d)", name)
        if not m: return None
        y, mo, d = int(m.group(1)), int(m.group(2)), int(m.group(3))
        try:
            return pd.Timestamp(year=y, month=mo, day=d, tz="UTC")
        except Exception:
            return None

    def _parse_yyyymm(self, name: str) -> Optional[Tuple[int,int]]:
        m = re.search(r"(?<!\d)(20\d{2})(0[1-9]|1[0-2])(?!\d)", name)
        if not m: return None
        return int(m.group(1)), int(m.group(2))

    def build_sea_ice_satellite_paths(self,
                                      dt0_str      : str,
                                      dtN_str      : str,
                                      hemisphere   : Optional[str] = None,
                                      institutions : Optional[Iterable[str]] = ("ESA","AWI"),
                                      sensors      : Optional[Iterable[str]] = None,          # ["cryosat2","envisat","sentinel3a","sentinel3b"]
                                      levels       : Optional[Iterable[str]] = None,          # ESA: ["L2P","L3C"]; AWI: ["l2p_release","l3cp_release"]
                                      versions     : Optional[Iterable[str]] = None,          # ESA only; and only "v2.0" or "v3.0"
                                      root_esa     : Optional[str] = None, 
                                      root_awi     : Optional[str] = None) -> List[Path]:
        """
        Build a de-duplicated, time-filtered list of NetCDF files for satellite sea-ice
        products across ESA and/or AWI directory layouts.

        Time filtering:
            - Returns files whose date/month lies within the closed interval
            [dt0_str, dtN_str], interpreted in UTC.
            - Daily products (e.g., ESA L2P) are filtered by date parsed from the filename
            (YYYYMMDD). Monthly products (e.g., ESA L3C or AWI monthly layout) are
            filtered by overlap between the (year, month) and [dt0, dtN].

        Supported directory patterns (walked automatically):
            ESA “flat” layout:
                <root_esa>/<L2P|L3C>/<sensor>/<hemisphere>/*.nc
                (hemisphere directory may be 'NH'/'SH' or lowercase)
                - L2P: filenames contain YYYYMMDD
                - L3C: filenames often contain YYYYMM

            ESA versioned layout:
                <root_esa>/<L2P|L3C>/<sensor>/<version>/<HEMISPHERE>/<YYYY>/<MM>/*.nc
                (HEMISPHERE can be 'NH'/'SH' or lowercase; <MM> subdirs may be absent)

            AWI release layout:
                <root_awi>/<l2p_release|l3cp_release>/<hemisphere>/<sensor>/<YYYY>[/<MM>]/*.nc

        Parameters
        ----------
        dt0_str, dtN_str
            Start and end of the desired time window (inclusive). Any pandas-parsable
            string is accepted (e.g., "2002-01-01", "2012-12-31"). Interpreted as UTC.
        hemisphere
            Hemisphere selector. If None, uses `self.hemisphere_dict['abbreviation']`.
            Case is normalized; both 'sh'/'nh' and 'SH'/'NH' are supported in dir scans.
        institutions
            Iterable of institutions to include (case-insensitive). Valid values include
            "ESA" and/or "AWI". If None or empty, both are searched.
        sensors
            Optional iterable of sensor names to filter (case-insensitive, matched to
            directory names under the institution layout).
        levels
            Optional iterable of level strings. ESA: "L2P", "L3C".
            AWI: "l2p_release", "l3cp_release". If "l2p"/"l3c" are provided for AWI
            they are mapped to "l2p_release"/"l3cp_release".
        versions
            ESA only: one or more version directory names (e.g., "v2.0", "v3.0").
            Case-insensitive. Ignored for AWI.
        root_esa, root_awi
            Root directories for ESA and AWI datasets.

        Returns
        -------
        List[Path]
            Sorted unique list of file paths within the requested time window that
            satisfy the filters.

        Notes
        -----
        - Filename parsing assumes dates follow YYYYMMDD (daily) or YYYYMM (monthly)
        within 2000–2099. Files without such tokens are kept only when they live in
        a year/month directory that overlaps [dt0, dtN].
        - If `institutions` is None or empty, both ESA and AWI are scanned.
        - All timestamps are treated as UTC to avoid tz-comparison issues.

        Examples
        --------
        >>> paths = self.build_sea_ice_satellite_paths(
        ...     "2002-01-01", "2002-12-31",
        ...     hemisphere="sh",
        ...     institutions=["ESA", "AWI"],
        ...     sensors=["cryosat2", "envisat"],
        ...     levels=["L2P","L3C","l2p_release"]
        ... )
        >>> len(paths), paths[:3]
        """
        dt0_str = dt0_str if dt0_str is not None else self.dt0_str
        dtN_str = dtN_str if dtN_str is not None else self.dtN_str
        t0 = pd.Timestamp(dt0_str, tz="UTC")
        t1 = pd.Timestamp(dtN_str, tz="UTC")
        if t1 < t0:
            raise ValueError("dtN_str date must be >= dt0_str date")
        hem = hemisphere if hemisphere is not None else self.hemisphere_dict["abbreviation"]
        hem_upper = str(hem).upper()
        hem_lower = str(hem).lower()
        hem_dirs_esa = [hem_upper, hem_lower]  # try both cases for ESA
        root_esa = root_esa if root_esa is not None else self.Sea_Ice_Obs_dict['ESA-CCI']
        root_awi = root_awi if root_awi is not None else self.Sea_Ice_Obs_dict['AWI']
        root_esa = Path(root_esa)
        root_awi = Path(root_awi)
        # filters (case-normalized)
        insts = set(s.upper() for s in (self._norm_list(institutions) or []))
        sensor_set = set(s.lower() for s in (self._norm_list(sensors) or [])) if sensors else None
        version_set = set(v.lower() for v in (self._norm_list(versions) or [])) if versions else None
        level_set = set(l.lower() for l in (self._norm_list(levels) or [])) if levels else None
        out: List[Path] = []
        # ---------------- ESA ----------------
        if not insts or "ESA" in insts:
            esa_levels = ["L2P", "L3C"] if not level_set else ([L.upper() for L in level_set if L.upper() in ("L2P", "L3C")] or ["L2P", "L3C"])
            for L in esa_levels:
                Ldir = Path(root_esa,L)
                if not Ldir.exists():
                    continue
                # <root>/<L>/<sensor>/...
                for sdir in sorted(p for p in Ldir.iterdir() if p.is_dir()):
                    sensor = sdir.name.lower()
                    if sensor_set and sensor not in sensor_set:
                        continue
                    # A) flat: .../L2P|L3C/<sensor>/<hem>/*.nc   (hem 'nh'/'sh' any case)
                    for hem_dir in hem_dirs_esa:
                        flat_dir = Path(sdir,hem_dir)
                        if flat_dir.exists():
                            for fp in flat_dir.glob("*.nc"):
                                if L == "L2P":
                                    ts = self._parse_yyyymmdd(fp.name)
                                    if ts is None or ts < t0 or ts > t1:
                                        continue
                                else:  # L3C
                                    ym = self._parse_yyyymm(fp.name)
                                    if ym and not self._month_overlap(ym[0], ym[1], t0, t1):
                                        continue
                                out.append(fp)
                    # B) versioned: .../L2P|L3C/<sensor>/<ver>/<HEM>/YYYY/MM/*.nc
                    for vdir in sorted(p for p in sdir.iterdir() if p.is_dir()):
                        vname = vdir.name.lower()
                        # skip if this 'vdir' is actually the hemisphere dir from flat layout
                        if vname in {"nh", "sh"}:
                            continue
                        if version_set and vname not in version_set:
                            continue
                        for hem_dir in hem_dirs_esa:
                            hdir = Path(vdir,hem_dir)
                            if not hdir.exists():
                                continue
                            for ydir in sorted(hdir.glob("[12][0-9][0-9][0-9]")):
                                y = int(ydir.name)
                                mdirs = list(sorted(ydir.glob("[01][0-9]")))
                                if mdirs:
                                    for mdir in mdirs:
                                        m = int(mdir.name)
                                        if not self._month_overlap(y, m, t0, t1):
                                            continue
                                        for fp in sorted(mdir.glob("*.nc")):
                                            if L == "L2P":
                                                ts = self._parse_yyyymmdd(fp.name)
                                                if ts is None or ts < t0 or ts > t1:
                                                    continue
                                            out.append(fp)
                                else:
                                    for fp in sorted(ydir.glob("*.nc")):
                                        if L == "L2P":
                                            ts = self._parse_yyyymmdd(fp.name)
                                            if ts is None or ts < t0 or ts > t1:
                                                continue
                                        else:
                                            ym = self._parse_yyyymm(fp.name)
                                            if ym and not self._month_overlap(ym[0], ym[1], t0, t1):
                                                continue
                                        out.append(fp)
        # ---------------- AWI ----------------
        if not insts or "AWI" in insts:
            if not level_set:
                awi_levels: List[str] = ["l2p_release", "l3cp_release"]
            else:
                # map friendly ESA-like inputs to AWI release names
                canonical = []
                for L in level_set:
                    l = L.lower()
                    if l in ("l2p_release", "l3cp_release"):
                        canonical.append(l)
                    elif l == "l2p":
                        canonical.append("l2p_release")
                    elif l == "l3c":
                        canonical.append("l3cp_release")
                awi_levels = sorted(set(canonical)) or ["l2p_release", "l3cp_release"]
            LHEM = hem_lower  # AWI uses lower-case hemisphere dirs
            for Lraw in awi_levels:
                Ldir = Path(root_awi,Lraw,LHEM)
                if not Ldir.exists():
                    continue
                for sdir in sorted(p for p in Ldir.iterdir() if p.is_dir()):
                    sensor = sdir.name.lower()
                    if sensor_set and sensor not in sensor_set:
                        continue
                    for ydir in sorted(sdir.glob("[12][0-9][0-9][0-9]")):
                        y = int(ydir.name)
                        mdirs = list(sorted(ydir.glob("[01][0-9]")))
                        if mdirs:
                            for mdir in mdirs:
                                m = int(mdir.name)
                                if not self._month_overlap(y, m, t0, t1):
                                    continue
                                out.extend(sorted(mdir.glob("*.nc")))
                        else:
                            # year-only dir; filter by yyyymm in filenames when available
                            for fp in sorted(ydir.glob("*.nc")):
                                ym = self._parse_yyyymm(fp.name)
                                if ym and not self._month_overlap(ym[0], ym[1], t0, t1):
                                    continue
                                out.append(fp)
        # de-dup & sort
        return sorted(set(out))
        
    def bin_esa_thickness_to_cice_grid(self, ESA_CCI, region=None):
        from scipy.stats import binned_statistic_2d
        from tqdm import tqdm
        x_edges = np.linspace(CICE_SO['TLON'].min().item(), CICE_SO['TLON'].max().item(), CICE_SO['TLON'].shape[1] + 1)
        y_edges = np.linspace(CICE_SO['TLAT'].min().item(), CICE_SO['TLAT'].max().item(), CICE_SO['TLAT'].shape[0] + 1)
        unique_days = np.unique(ESA_CCI['time'].dt.floor('D').values)
        sit_list = []
        sigma_list = []
        time_list = []
        for day in tqdm(unique_days, desc="Binning ESA CCI by day"):
            ESA_day = ESA_CCI.where(ESA_CCI['time'].dt.floor('D') == day, drop=True)
            lon = ESA_day.lon.values
            lat = ESA_day.lat.values
            sit = ESA_day.sea_ice_thickness.values
            sit_sigma = ESA_day.sea_ice_thickness_uncertainty.values
            mask = np.isfinite(lon) & np.isfinite(lat) & np.isfinite(sit) & np.isfinite(sit_sigma)
            lon = lon[mask]
            lat = lat[mask]
            sit = sit[mask]
            sit_sigma = sit_sigma[mask]
            binned_sit, _, _, _ = binned_statistic_2d(lon, lat, sit, statistic="mean", bins=[x_edges, y_edges])
            binned_sigma, _, _, _ = binned_statistic_2d(lon, lat, sit_sigma, statistic="mean", bins=[x_edges, y_edges])
            sit_list.append(binned_sit.T[None, ...])
            sigma_list.append(binned_sigma.T[None, ...])
            time_list.append(np.datetime64(day))
        ESA_sit_reG = xr.DataArray(data   = np.concatenate(sit_list, axis=0),
                                   dims   = ("time", "nj", "ni"),
                                   coords = {"time": time_list,
                                             "TLON": (("nj", "ni"), CICE_SO.TLON.values),
                                             "TLAT": (("nj", "ni"), CICE_SO.TLAT.values)},
                                   name   = "ESA_sit")
        ESA_sit_sigma_reG = xr.DataArray(data   = np.concatenate(sigma_list, axis=0),
                                         dims   = ("time", "nj", "ni"),
                                         coords = {"time": time_list,
                                                   "TLON": (("nj", "ni"), CICE_SO.TLON.values),
                                                   "TLAT": (("nj", "ni"), CICE_SO.TLAT.values)},
                                        name    = "ESA_sit_sigma")
        return xr.Dataset({"ESA_sit": ESA_sit_reG, "ESA_sit_sigma": ESA_sit_sigma_reG})