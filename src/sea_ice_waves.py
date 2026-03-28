from __future__ import annotations
import os
import calendar
from pathlib import Path
from typing import Optional, Tuple
import numpy as np
import pandas as pd
import xarray as xr

__all__ = ["SeaIceWaves"]

class SeaIceWaves:
    """
    Wave-spectrum preprocessing utilities for CICE6/Icepack wave forcing.

    This mixin downloads or opens monthly CAWCR station spectral files, collapses
    directional spectra to 1D frequency spectra E(f), remaps station spectra onto
    the native CICE T grid, and writes NetCDF files in the exact layout expected by
    `ice_read_nc_xyf()` in the current CICE fork:

        efreq(ni, nj, nfreq, time)

    Notes
    -----
    1. This class assumes your toolbox provides:
       - self.logger
       - self.Waves_dict
       - self.CICE_dict
       - self.load_cice_grid()
       - the sparse regridding helpers added to SeaIceRegridder

    2. The current Fortran `get_wave_spec()` in your CICE fork reads only record 1.
       These files are therefore "correctly formatted for CICE", but a small Fortran
       follow-up is still needed for time-varying hourly wave forcing.
    """

    def __init__(self, **kwargs):
        return

    # ------------------------------------------------------------------
    # basic path / metadata helpers
    # ------------------------------------------------------------------
    def _require_waves_dict(self):
        if not hasattr(self, "Waves_dict"):
            raise AttributeError("Missing Waves_dict in configuration JSON.")
        required = [
            "cawcr_spec_url_fmt",
            "D_processed",
            "wave_file_fmt",
            "P_reG_cawcr2cice_weights",
            "time_dim",
            "station_dim",
            "freq_dim",
            "dir_dim",
            "station_lon_name",
            "station_lat_name",
            "variance_density_name",
            "k_nearest",
            "idw_power",
            "radius_km",
        ]
        missing = [k for k in required if k not in self.Waves_dict]
        if missing:
            raise KeyError(f"Waves_dict missing required keys: {missing}")

    def _month_bounds(self, year: int, month: int) -> Tuple[pd.Timestamp, pd.Timestamp]:
        ndays = calendar.monthrange(year, month)[1]
        t0 = pd.Timestamp(year=year, month=month, day=1)
        t1 = pd.Timestamp(year=year, month=month, day=ndays, hour=23, minute=59, second=59)
        return t0, t1

    def _cawcr_url(self, year: int, month: int) -> str:
        self._require_waves_dict()
        yyyymm = f"{year:04d}{month:02d}"
        return self.Waves_dict["cawcr_spec_url_fmt"].format(
            year=year,
            month=f"{month:02d}",
            yyyymm=yyyymm,
        )

    def _processed_wave_path(self, year: int, month: int) -> Path:
        self._require_waves_dict()
        yyyymm = f"{year:04d}{month:02d}"
        d_out = Path(self.Waves_dict["D_processed"])
        d_out.mkdir(parents=True, exist_ok=True)
        fname = self.Waves_dict["wave_file_fmt"].format(
            year=year,
            month=f"{month:02d}",
            yyyymm=yyyymm,
        )
        return d_out / fname

    # ------------------------------------------------------------------
    # CAWCR input
    # ------------------------------------------------------------------
    def open_cawcr_month(self, year: int, month: int, chunks: Optional[dict] = None) -> xr.Dataset:
        """
        Open a monthly CAWCR station-spectrum file via OPeNDAP/THREDDS.
        """
        url = self._cawcr_url(year, month)
        chunks = chunks if chunks is not None else self.Waves_dict.get("chunks", {})
        engine = self.Waves_dict.get("thredds_engine", "netcdf4")

        self.logger.info(f"Opening CAWCR spectra: {url}")
        ds = xr.open_dataset(url, engine=engine, chunks=chunks)
        return ds

    def _find_first_present(self, ds: xr.Dataset, candidates, kind: str) -> str:
        for name in candidates:
            if name in ds.variables or name in ds.coords or name in ds.dims:
                return name
        raise KeyError(f"Could not find {kind}; tried {candidates}")

    def normalise_cawcr_names(self, ds: xr.Dataset) -> xr.Dataset:
        """
        Standardise the expected CAWCR station-spectrum variable and dimension names.
        """
        wd = self.Waves_dict

        station_dim = self._find_first_present(ds, [wd["station_dim"], "station"], "station dim")
        time_dim    = self._find_first_present(ds, [wd["time_dim"], "time"], "time dim")
        freq_dim    = self._find_first_present(ds, [wd["freq_dim"], "frequency", "freq"], "frequency dim")
        dir_dim     = self._find_first_present(ds, [wd["dir_dim"], "direction", "dir"], "direction dim")

        efth_name = self._find_first_present(
            ds,
            [wd["variance_density_name"], "efth", "Efth"],
            "directional spectral density variable",
        )
        lon_name = self._find_first_present(
            ds,
            [wd["station_lon_name"], "longitude", "lon"],
            "station longitude variable",
        )
        lat_name = self._find_first_present(
            ds,
            [wd["station_lat_name"], "latitude", "lat"],
            "station latitude variable",
        )

        rename_map = {}
        if station_dim != "station":
            rename_map[station_dim] = "station"
        if time_dim != "time":
            rename_map[time_dim] = "time"
        if freq_dim != "frequency":
            rename_map[freq_dim] = "frequency"
        if dir_dim != "direction":
            rename_map[dir_dim] = "direction"
        if lon_name != "station_lon":
            rename_map[lon_name] = "station_lon"
        if lat_name != "station_lat":
            rename_map[lat_name] = "station_lat"
        if efth_name != "efth":
            rename_map[efth_name] = "efth"

        ds = ds.rename(rename_map)

        # Make station lon/lat proper coordinates if needed
        if "station_lon" in ds and "station_lon" not in ds.coords:
            ds = ds.set_coords("station_lon")
        if "station_lat" in ds and "station_lat" not in ds.coords:
            ds = ds.set_coords("station_lat")

        return ds

    # ------------------------------------------------------------------
    # directional -> 1D frequency spectrum
    # ------------------------------------------------------------------
    def _as_station_1d_coord(self, da: xr.DataArray, coord_name: str) -> xr.DataArray:
        """
        Coerce a station coordinate to a 1D DataArray with dims ('station',).

        Handles cases where CAWCR exposes lon/lat with an extra singleton dim or
        with a time dimension.
        """
        coord = da[coord_name]
        # already good
        if coord.dims == ("station",):
            return coord
        # strip any singleton dimensions
        squeeze_dims = [d for d in coord.dims if d != "station" and coord.sizes[d] == 1]
        if squeeze_dims:
            coord = coord.squeeze(squeeze_dims, drop=True)
        # common case: coordinate varies only trivially over time
        if "time" in coord.dims and "station" in coord.dims:
            coord = coord.isel(time=0, drop=True)
        # if station is not the only remaining dimension, fail loudly
        other_dims = [d for d in coord.dims if d != "station"]
        if other_dims:
            raise ValueError(f"{coord_name} must be 1D over station after cleanup; got dims {coord.dims}")
        return coord

    def collapse_directional_spectrum(self, ds: xr.Dataset) -> xr.DataArray:
        """
        Convert directional variance density Efth(f,theta) into 1D E(f) by integrating
        over direction.

        Returns
        -------
        xr.DataArray
            efreq_station(time, station, frequency)
        """
        ds = self.normalise_cawcr_names(ds)
        if "efth" not in ds:
            raise KeyError("Expected directional spectrum variable 'efth' after renaming.")
        theta = ds["direction"]
        theta_vals = theta.values.astype(float)
        # CAWCR directions are commonly degrees; convert to radians for integration.
        if np.nanmax(np.abs(theta_vals)) > 2 * np.pi + 1e-6:
            theta_rad = np.deg2rad(theta_vals)
        else:
            theta_rad = theta_vals
        dtheta = np.gradient(theta_rad)
        dtheta_da = xr.DataArray(dtheta,
                                 dims=("direction",),
                                 coords={"direction": ds["direction"]},
                                 name="dtheta")
        efreq = (ds["efth"] * dtheta_da).sum("direction", skipna=True)
        efreq = efreq.transpose("time", "station", "frequency")
        efreq.name = "efreq_station"
        efreq.attrs.update({"long_name": "direction-integrated wave energy spectrum",
                            "units": ds["efth"].attrs.get("units", "m2/Hz")})
        efreq = efreq.assign_coords(station_lon=ds["station_lon"],
                                    station_lat=ds["station_lat"])
        efreq = efreq.transpose("time", "station", "frequency")
        if hasattr(efreq, "chunk"):
            efreq = efreq.chunk({"time"      : int(self.Waves_dict.get("time_chunk", 24)),
                                 "station"   : int(self.Waves_dict.get("station_chunk", 256)),
                                 "frequency" : -1})
        return efreq

    # ------------------------------------------------------------------
    # target grid helpers
    # ------------------------------------------------------------------
    def get_cice_tgrid_lonlat(self) -> Tuple[xr.DataArray, xr.DataArray]:
        """
        Get native CICE T-grid longitudes and latitudes as 2D DataArrays on (nj, ni).
        """
        if hasattr(self, "load_cice_grid"):
            self.load_cice_grid()

        # Preferred path: already-populated grid objects
        if hasattr(self, "G_t") and self.G_t is not None:
            if isinstance(self.G_t, xr.Dataset):
                if "lon" in self.G_t and "lat" in self.G_t:
                    lon = self.G_t["lon"]
                    lat = self.G_t["lat"]
                    return lon, lat

        # Fallback: read raw CICE grid file directly
        P_G = self.CICE_dict["P_G"]
        tlon_name, tlat_name = self.CICE_dict.get("tcoord_names", ["TLON", "TLAT"])
        dsG = xr.open_dataset(P_G)
        if tlon_name not in dsG or tlat_name not in dsG:
            raise KeyError(f"Could not find {tlon_name}/{tlat_name} in {P_G}")
        lon = dsG[tlon_name]
        lat = dsG[tlat_name]
        return lon, lat

    def get_cice_ocean_mask(self) -> Optional[xr.DataArray]:
        """
        Return a native-grid ocean mask on (nj, ni) if available.
        """
        P_KMT = self.CICE_dict.get("P_KMT")
        if P_KMT and os.path.exists(P_KMT):
            ds = xr.open_dataset(P_KMT)
            if "kmt_org" in ds:
                return xr.where(ds["kmt_org"] > 0, 1, 0)
            if "kmt" in ds:
                return xr.where(ds["kmt"] > 0, 1, 0)
        return None

    # ------------------------------------------------------------------
    # main regridding workflow
    # ------------------------------------------------------------------
    def build_or_load_cawcr_to_cice_weights(self, src_lon: np.ndarray, src_lat: np.ndarray, tgt_lon: xr.DataArray, tgt_lat: xr.DataArray,
                                            target_mask : Optional[xr.DataArray] = None,
                                            overwrite   : bool                   = False):
        """
        Build or load sparse IDW weights from CAWCR stations to the native CICE T grid.
        """
        p_weights = self.Waves_dict["P_reG_cawcr2cice_weights"]
        k         = int(self.Waves_dict["k_nearest"])
        power     = float(self.Waves_dict["idw_power"])
        radius_km = float(self.Waves_dict["radius_km"])
        return self.build_or_load_station_to_curvilinear_sparse_weights(src_lon=src_lon, src_lat=src_lat, tgt_lon=tgt_lon.values, tgt_lat=tgt_lat.values,
                                                                        p_weights   = p_weights,
                                                                        target_mask = None if target_mask is None else target_mask.values,
                                                                        k           = k,
                                                                        power       = power,
                                                                        radius_km   = radius_km,
                                                                        overwrite   = overwrite)

    def regrid_station_spectra_to_cice(self, efreq_station: xr.DataArray, overwrite_weights: bool = False) -> xr.DataArray:
        """
        Regrid station-based 1D spectra to the native CICE T grid.

        Input
        -----
        efreq_station(time, station, frequency)

        Output
        ------
        efreq_grid(time, nj, ni, frequency)
        """
        if tuple(efreq_station.dims) != ("time", "station", "frequency"):
            efreq_station = efreq_station.transpose("time", "station", "frequency")
        tgt_lon, tgt_lat = self.get_cice_tgrid_lonlat()
        tgt_mask = self.get_cice_ocean_mask()
        # active target cells
        tgt_active = np.isfinite(tgt_lon.values) & np.isfinite(tgt_lat.values)
        if tgt_mask is not None:
            tgt_active &= tgt_mask.values.astype(bool)
        if getattr(self, "hemisphere", "south").lower().startswith("s"):
            tgt_active &= tgt_lat.values <= float(self.Waves_dict.get("target_lat_max", -35.0))
        else:
            tgt_active &= tgt_lat.values >= float(self.Waves_dict.get("target_lat_min", 35.0))
        station_lon = self._as_station_1d_coord(efreq_station, "station_lon")
        station_lat = self._as_station_1d_coord(efreq_station, "station_lat")
        src_lon     = station_lon.values.astype(float)
        src_lat     = station_lat.values.astype(float)
        # build a labeled 1D station mask
        src_keep    = np.isfinite(src_lon) & np.isfinite(src_lat)
        if getattr(self, "hemisphere", "south").lower().startswith("s"):
            src_keep &= src_lat <= float(self.Waves_dict.get("source_lat_max", -20.0))
        else:
            src_keep &= src_lat >= float(self.Waves_dict.get("source_lat_min", 20.0))
        src_keep      = xr.DataArray(src_keep, dims=("station",), coords={"station": efreq_station["station"]})
        efreq_station = efreq_station.isel(station=src_keep)
        src_lon       = self._as_station_1d_coord(efreq_station, "station_lon").values.astype(float)
        src_lat       = self._as_station_1d_coord(efreq_station, "station_lat").values.astype(float)
        self.logger.info(f"station_lon dims: {efreq_station['station_lon'].dims}, {efreq_station['station_lon'].shape}")
        self.logger.info(f"station_lat dims: {efreq_station['station_lat'].dims}, {efreq_station['station_lat'].shape}")
        W = self.build_or_load_cawcr_to_cice_weights(src_lon=src_lon, src_lat=src_lat, tgt_lon=tgt_lon, tgt_lat=tgt_lat,
                                                     target_mask = xr.DataArray(tgt_active, dims=tgt_lon.dims, coords=tgt_lon.coords),
                                                     overwrite   = overwrite_weights)
        self.logger.info(f"efreq_station shape  = {efreq_station.shape}")
        self.logger.info(f"efreq_station chunks = {getattr(efreq_station.data, 'chunks', None)}")
        out = self.apply_sparse_station_regridder(values       = efreq_station,                                   # (time, station, frequency)
                                                  weights      = W,                                              # sparse (n_target, n_station)
                                                  target_shape = tgt_lon.shape,
                                                  fill_value   = float(self.Waves_dict.get("fill_value", 0.0)),
                                                  time_chunk   = int(self.Waves_dict.get("time_chunk", 24)))
        da = xr.DataArray(out,
                          dims   = ("time", "nj", "ni", "frequency"),
                          coords = {"time": efreq_station["time"],
                                    "frequency": efreq_station["frequency"],
                                    "TLON": (("nj", "ni"), tgt_lon.values),
                                    "TLAT": (("nj", "ni"), tgt_lat.values)},
                          name   = "efreq_grid",
                          attrs  = {"long_name": "CAWCR station spectra regridded to native CICE T grid",
                                    "units": efreq_station.attrs.get("units", "m2/Hz")})
        return da

    # ------------------------------------------------------------------
    # output writer for CICE
    # ------------------------------------------------------------------
    def build_cice_wave_dataset(self, efreq_grid: xr.DataArray) -> xr.Dataset:
        """
        Convert efreq_grid(time,nj,ni,frequency) into a dataset with

            efreq(time,nfreq,nj,ni)

        This matches the dimension order that the current CICE/NetCDF Fortran
        read path is effectively seeing for successful slicing via ice_read_nc_xyf().
        """
        if tuple(efreq_grid.dims) != ("time", "nj", "ni", "frequency"):
            efreq_grid = efreq_grid.transpose("time", "nj", "ni", "frequency")
        ni_name    = self.CICE_dict.get("x_dim", "ni")
        nj_name    = self.CICE_dict.get("y_dim", "nj")
        nfreq_name = "nfreq"
        time_name  = "time"
        ni         = np.arange(efreq_grid.sizes["ni"], dtype=np.int32)
        nj         = np.arange(efreq_grid.sizes["nj"], dtype=np.int32)
        # coordinate values for frequency bins
        wavefreq  = efreq_grid["frequency"].values.astype(np.float32)
        # rename frequency dim and reorder to time,nfreq,nj,ni
        efreq_out = (efreq_grid.rename({"frequency": nfreq_name}).transpose(time_name, nfreq_name, "nj", "ni").astype(np.float32))
        ds_out    = xr.Dataset(data_vars = {"efreq"   : ((time_name, nfreq_name, nj_name, ni_name), efreq_out.values),
                                            "TLON"    : ((nj_name, ni_name)                       , efreq_grid["TLON"].values.astype(np.float32)),
                                            "TLAT"    : ((nj_name, ni_name)                       , efreq_grid["TLAT"].values.astype(np.float32)),
                                            "wavefreq": ((nfreq_name,)                            , wavefreq)},
                               coords    = {time_name : efreq_grid["time"].values,
                                            nfreq_name: wavefreq,
                                            nj_name   : nj,
                                            ni_name   : ni},
                               attrs     = {"title"   : "CAWCR wave spectra remapped to native CICE grid",
                                            "source"  : "CAWCR monthly station spectra, directionally integrated and regridded",
                                            "note"    : "efreq written with dims (time,nfreq,nj,ni) for compatibility with current ice_read_nc_xyf() behavior"})
        # Optional dwavefreq
        fq = ds_out["wavefreq"].values.astype(np.float64)
        if fq.size > 1:
            ds_out["dwavefreq"] = xr.DataArray(np.gradient(fq).astype(np.float32), dims=(nfreq_name,))
        # Apply fill over non-finite / negative values
        fill_value           = np.float32(self.Waves_dict.get("fill_value", 0.0))
        ef                   = ds_out["efreq"].values
        ef[~np.isfinite(ef)] = fill_value
        ef[ef < 0.0]         = 0.0
        ds_out["efreq"][:]   = ef.astype(np.float32, copy=False)
        ds_out["efreq"].attrs.update({"long_name": "wave energy spectrum", "units": "m2/Hz"})
        # xarray wants _FillValue in encoding, not attrs
        for vname in ds_out.variables:
            ds_out[vname].attrs.pop("_FillValue", None)
            ds_out[vname].attrs.pop("missing_value", None)
        return ds_out

    def write_cice_wave_netcdf(self, ds_out, year, month, overwrite=False):
        """
        Write the monthly CICE-ready wave NetCDF.
        """
        p_out = self._processed_wave_path(year, month)
        if p_out.exists() and not overwrite:
            self.logger.info(f"Wave file exists; reusing {p_out}")
            return p_out
        ds_out = ds_out.copy()
        for vname in ds_out.variables:
            ds_out[vname].attrs.pop("_FillValue", None)
            ds_out[vname].attrs.pop("missing_value", None)
        fill_value = np.float32(self.Waves_dict.get("fill_value", 0.0))
        nfreq      = ds_out.sizes["nfreq"]
        encoding   = {"efreq"   : {"dtype"      : "float32",
                                   "zlib"       : True,
                                   "complevel"  : 1,
                                   "_FillValue" : fill_value,
                                   "chunksizes" : (1, nfreq, 180, 180)},
                      "TLON"    : {"dtype"       : "float32",
                                   "zlib"        : True,
                                   "complevel"   : 1,
                                   "chunksizes"  : (180, 180)},
                      "TLAT"    : {"dtype"       : "float32",
                                   "zlib"        : True,
                                   "complevel"   : 1,
                                   "chunksizes"  : (180, 180)},
                      "wavefreq": {"dtype"       : "float32"}}
        if "dwavefreq" in ds_out:
            encoding["dwavefreq"] = {"dtype": "float32"}
        self.logger.info(f"Writing CICE wave forcing: {p_out}")
        ds_out.to_netcdf(p_out,
                         format         = "NETCDF4",
                         engine         = "netcdf4",
                         encoding       = encoding,
                         unlimited_dims = ["time"])
        return p_out

    def reorder_cawcr_efreq_file(self, p_in,
                                 p_out      = None,
                                 efreq_var  = "efreq",
                                 time_chunk = 24,
                                 zlib       = True,
                                 complevel  = 1,
                                 overwrite  = False):
        """
        Reorder a CAWCR-for-CICE NetCDF from:

            efreq(ni, nj, nfreq, time)

        to:

            efreq(time, nfreq, nj, ni)

        using chunked reads/writes so the whole array is never loaded at once.

        Parameters
        ----------
        p_in : str or Path
            Input NetCDF path.
        p_out : str or Path, optional
            Output NetCDF path. Defaults to '<stem>_reordered.nc'.
        efreq_var : str, default 'efreq'
            Name of the spectral variable.
        time_chunk : int, default 24
            Number of time records to process per chunk.
        zlib : bool, default True
            Whether to compress output efreq.
        complevel : int, default 1
            NetCDF compression level.
        overwrite : bool, default False
            Overwrite output if it exists.

        Returns
        -------
        Path
            Output file path.
        """
        from netCDF4 import Dataset
        p_in = Path(p_in)
        if p_out is None:
            p_out = p_in.with_name(f"{p_in.stem}_reordered.nc")
        p_out = Path(p_out)
        if p_out.exists():
            if overwrite:
                p_out.unlink()
            else:
                raise FileExistsError(f"{p_out} already exists")
        print(f"opening {p_out} for writing to via NETCDF4 library ...")
        with Dataset(p_in, "r") as src, Dataset(p_out, "w", format="NETCDF4") as dst:
            # ----------------------------
            # dimensions
            # ----------------------------
            for dname, dim in src.dimensions.items():
                if dname == "time":
                    dst.createDimension(dname, None)
                else:
                    dst.createDimension(dname, len(dim))
                print(f"{dname} dim created")
            # ----------------------------
            # copy non-efreq variables first
            # ----------------------------
            for vname, var in src.variables.items():
                if vname == efreq_var:
                    continue
                fill_value = getattr(var, "_FillValue", None)
                kwargs = {}
                if fill_value is not None:
                    kwargs["fill_value"] = fill_value
                out_var = dst.createVariable(vname, var.dtype, var.dimensions,
                                             zlib      = (zlib if vname in ("TLON", "TLAT") else False),
                                             complevel = (complevel if vname in ("TLON", "TLAT") else 0), **kwargs)
                print(f"{vname} variable created")
                # copy attrs except _FillValue which is handled above
                for attr in var.ncattrs():
                    if attr == "_FillValue":
                        continue
                    out_var.setncattr(attr, getattr(var, attr))
                out_var[:] = var[:]
            # ----------------------------
            # create reordered efreq variable
            # input dims expected: (ni, nj, nfreq, time)
            # output dims written: (time, nfreq, nj, ni)
            # ----------------------------
            vin = src.variables[efreq_var]
            assert vin.dimensions == ("ni", "nj", "nfreq", "time"), (f"Unexpected input dims for {efreq_var}: {vin.dimensions}")
            fill_value = getattr(vin, "_FillValue", np.float32(0.0))
            vout = dst.createVariable(efreq_var, vin.dtype, ("time", "nfreq", "nj", "ni"),
                                      zlib       = zlib,
                                      complevel  = complevel,
                                      shuffle    = True,
                                      chunksizes = (1, len(src.dimensions["nfreq"]), 180, 180),
                                      fill_value = fill_value)
            print(f"{efreq_var} BIG var created")
            for attr in vin.ncattrs():
                if attr == "_FillValue":
                    continue
                vout.setncattr(attr, getattr(vin, attr))
            ntime = len(src.dimensions["time"])
            print(f"Reordering {efreq_var}: {p_in.name} -> {p_out.name}")
            print(f"ntime={ntime}, chunk={time_chunk}")
            for t0 in range(0, ntime, time_chunk):
                t1                     = min(ntime, t0 + time_chunk)
                print(f"  time slice {t0}:{t1}")
                # read chunk: (ni, nj, nfreq, chunk)
                arr                    = vin[:, :, :, t0:t1]
                # transpose to (chunk, nfreq, nj, ni)
                arr                    = np.transpose(arr, (3, 2, 1, 0))
                # sanitize
                arr                    = np.asarray(arr)
                arr[~np.isfinite(arr)] = fill_value
                arr[arr < 0]           = 0
                vout[t0:t1, :, :, :]   = arr
            # ----------------------------
            # global attrs
            # ----------------------------
            for attr in src.ncattrs():
                dst.setncattr(attr, getattr(src, attr))
            dst.setncattr("note", "efreq reordered from (ni,nj,nfreq,time) to (time,nfreq,nj,ni) for CICE wave forcing compatibility")
        return p_out

    # ------------------------------------------------------------------
    # APIs
    # ------------------------------------------------------------------
    def prepare_cawcr_wave_month(self,
                                 year             : int,
                                 month            : int,
                                 overwrite        : bool = False,
                                 overwrite_weights: bool = False) -> Path:
        """
        End-to-end monthly CAWCR -> native CICE spectral forcing pipeline.
        """
        ds_raw        = self.open_cawcr_month(year, month)
        efreq_station = self.collapse_directional_spectrum(ds_raw)
        efreq_grid    = self.regrid_station_spectra_to_cice(efreq_station, overwrite_weights=overwrite_weights)
        ds_out        = self.build_cice_wave_dataset(efreq_grid)
        return self.write_cice_wave_netcdf(ds_out, year, month, overwrite=overwrite)

    def prepare_cawcr_wave_range(self,
                                 dt0               : str,
                                 dtN               : str,
                                 overwrite         : bool = False,
                                 overwrite_weights : bool = False) -> list[Path]:
        """
        Process a monthly range from dt0 to dtN inclusive.
        """
        t0 = pd.Timestamp(dt0).to_period("M")
        t1 = pd.Timestamp(dtN).to_period("M")
        months = pd.period_range(t0, t1, freq="M")
        outputs = []
        for p in months:
            self.logger.info(f"Preparing CAWCR wave month: {p.year:04d}-{p.month:02d}")
            outputs.append( self.prepare_cawcr_wave_month(p.year, p.month, overwrite=overwrite, overwrite_weights=overwrite_weights))
        return outputs
