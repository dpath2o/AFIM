import os, gc, re, dask
import xarray            as xr
import numpy             as np
import pandas            as pd
import xesmf             as xe
from pathlib             import Path
from pyproj              import Transformer
from pyresample.geometry import AreaDefinition, SwathDefinition
from pyresample.kd_tree  import resample_nearest

__all__ = ["SeaIceGridWork"]

class SeaIceGridWork:

    def __init__(self, **kwargs):
        return

    def normalise_longitudes(self, lon, to="0-360", eps=1e-12):
        """
        Wrap longitudes to either [0, 360) or (-180, 180].
        Works with numpy arrays, scalars, and xarray objects.
        """
        # First get [0, 360)
        lon_wrapped = ((lon % 360) + 360) % 360  # safe for negatives, NaNs pass through
        if to == "0-360":
            # Collapse values extremely close to 360 back to 0
            if isinstance(lon_wrapped, xr.DataArray):
                lon_wrapped = xr.where(np.isclose(lon_wrapped, 360.0, atol=eps), 0.0, lon_wrapped)
            else:
                lon_wrapped = np.where(np.isclose(lon_wrapped, 360.0, atol=eps), 0.0, lon_wrapped)
            return lon_wrapped
        elif to == "-180-180":
            lon_180 = ((lon_wrapped + 180.0) % 360.0) - 180.0  # -> (-180, 180]
            # Prefer [-180, 180) by mapping exactly 180 to -180
            if isinstance(lon_180, xr.DataArray):
                lon_180 = xr.where(np.isclose(lon_180, 180.0, atol=eps), -180.0, lon_180)
            else:
                lon_180 = np.where(np.isclose(lon_180, 180.0, atol=eps), -180.0, lon_180)
            return lon_180
        else:
            raise ValueError("to must be '0-360' or '-180-180'")

    def reapply_landmask(self, DS, apply_unmodified=False):
        """

        Apply landmask to all spatial variables in the dataset.

        Uses the modified bathymetry (`kmt_mod`) field from the model grid to mask out land cells.
        Applies this mask to all variables with dimensions `("nj", "ni")`, ensuring that land areas
        are excluded from any subsequent analysis or output.

        INPUTS:
           DS : xarray.Dataset; Dataset containing sea ice fields to be masked.

        OUTPUTS:
           xarray.Dataset; Same dataset with land cells set to NaN for spatial variables.

        """
        if self.use_gi and not apply_unmodified:
            kmt_mask = xr.open_dataset(self.P_KMT_mod)['kmt'] == 1
        else:
            kmt_mask = xr.open_dataset(self.P_KMT_org)['kmt'] == 1
        for var in DS.data_vars:
            da = DS[var]
            if {"nj", "ni"}.issubset(da.dims):  # Only apply to spatial fields
                self.logger.debug(f"Masking land for variable: {var}")
                DS[var] = da.where(kmt_mask)
        self.logger.info("Applied landmask to rolled dataset")
        return DS

    def build_grid_corners(self, lat_rads, lon_rads, grid_res=0.25, source_in_radians=True):
        if source_in_radians:
            lon = self.radians_to_degrees(lon_rads)
            lat = self.radians_to_degrees(lat_rads)
        else:
            lon = lon_rads
            lat = lat_rads
        ny, nx               = lat.shape
        lat_b                = np.zeros((ny + 1, nx + 1))
        lon_b                = np.zeros((ny + 1, nx + 1))
        lat_b[1:-1, 1:-1]    = grid_res * (lat[:-1, :-1] + lat[:-1, 1:] + lat[1:, :-1] + lat[1:, 1:])
        lon_b[1:-1, 1:-1]    = grid_res * (lon[:-1, :-1] + lon[:-1, 1:] + lon[1:, :-1] + lon[1:, 1:])
        lat_b[[0, -1], 1:-1] = lat[[0, -1], :-1]
        lon_b[[0, -1], 1:-1] = lon[[0, -1], :-1]
        lat_b[1:-1, [0, -1]] = lat[:-1, [0, -1]]
        lon_b[1:-1, [0, -1]] = lon[:-1, [0, -1]]
        lat_b[0, 0]          = lat[0, 0]
        lat_b[0, -1]         = lat[0, -1]
        lat_b[-1, 0]         = lat[-1, 0]
        lat_b[-1, -1]        = lat[-1, -1]
        lon_b[0, 0]          = lon[0, 0]
        lon_b[0, -1]         = lon[0, -1]
        lon_b[-1, 0]         = lon[-1, 0]
        lon_b[-1, -1]        = lon[-1, -1]
        lon_b                = self.normalise_longitudes(lon_b)
        return lon_b, lat_b

    def define_cice_grid(self, grid_type='t', mask=False, build_grid_corners=False, slice_hem=False):
        """
        Build a minimal CICE grid dataset (lon/lat/area) with optional corners and mask.

        Parameters
        ----------
        grid_type : {'t','u'}, default 't'
            Which CICE grid to load from ``self.CICE_dict["P_G"]`` (e.g., `tlat/tlon/tarea`
            or `ulat/ulon/uarea`). Longitudes are normalised to 0–360°.
        mask : bool, default False
            If True, include a landmask variable `mask` derived from either the
            grounded-iceberg-modified mask (`self.P_KMT_mod`) when `self.use_gi` is True,
            or the original mask (`self.P_KMT_org`) otherwise.
        build_grid_corners : bool, default False
            If True, compute and attach corner coordinates `lon_b`, `lat_b` (dims
            `('nj_b','ni_b')`) using `build_grid_corners`.
        slice_hem : bool, default False
            If True, apply `slice_hemisphere()` to the returned dataset.

        Returns
        -------
        xr.Dataset
            Dataset with variables:
            - ``area`` (m²),
            - ``lon`` , ``lat`` (degrees),
            - optional ``lon_b``, ``lat_b`` (degrees),
            - optional ``mask`` (binary),
            and coordinates ``('nj','ni')`` (and ``'nj_b','ni_b'`` if corners requested).

        Notes
        -----
        - Dimension names come from ``self.CICE_dict['spatial_dims']`` for centres
        and are set to ``('nj_b','ni_b')`` for corners.
        - Set ``slice_hem=True`` if you want to stay within the configured hemisphere.
        """
        std_dim_names = self.CICE_dict['spatial_dims']
        G             = xr.open_dataset(self.CICE_dict["P_G"])
        lon_rads      = G[f"{grid_type}lon"].values
        lat_rads      = G[f"{grid_type}lat"].values
        ny, nx        = lon_rads.shape
        lon           = self.normalise_longitudes(self.radians_to_degrees(lon_rads))
        lat           = self.radians_to_degrees(lat_rads)
        area          = G[f'{grid_type}area'].values
        data_vars     = {"area": (std_dim_names, area),
                         "lon":  (std_dim_names, lon),
                         "lat":  (std_dim_names, lat)}
        coords        = {"nj": np.arange(ny),
                         "ni": np.arange(nx)}
        if build_grid_corners:
            crn_dim_names      = ("nj_b", "ni_b")
            lon_b, lat_b       = self.build_grid_corners(lon_rads, lat_rads)
            data_vars["lon_b"] = (crn_dim_names, lon_b)
            data_vars["lat_b"] = (crn_dim_names, lat_b)
            coords["nj_b"]     = np.arange(ny + 1)
            coords["ni_b"]     = np.arange(nx + 1)
        if mask:
            mask_data = xr.open_dataset(self.P_KMT_mod if self.use_gi else self.P_KMT_org).kmt.values
            data_vars["mask"] = (std_dim_names, mask_data)
        ds = xr.Dataset(data_vars=data_vars, coords=coords)
        if slice_hem:
            ds = self.slice_hemisphere(ds)
        return ds

    def load_bgrid(self, slice_hem=False):
        """
        Load and construct the B-grid datasets (t-grid and u-grid), including native and boundary coordinates,
        converted to degrees and standardized to [-180, 180]. Stores the resulting datasets in `self.G_t` and `self.G_u`.

        Parameters
        ----------
        slice_hem : bool, default=False
            If True, apply hemispheric slicing to the loaded datasets using the defined hemisphere.

        Sets
        ----
        self.G_t : xarray.Dataset
            Dataset containing t-grid variables and coordinates.

        self.G_u : xarray.Dataset
            Dataset containing u-grid variables and coordinates.

        self.bgrid_loaded : bool
            Flag indicating successful B-grid load.
        """
        G       = xr.open_dataset(self.CICE_dict['P_G'])
        KMT_org = xr.open_dataset(self.P_KMT_org).kmt.data
        KMT_mod = xr.open_dataset(self.P_KMT_mod).kmt.data if self.use_gi else KMT_org
        TLAT    = self.radians_to_degrees(G['tlat'].data)
        TLON    = self.radians_to_degrees(G['tlon'].data)
        ULAT    = self.radians_to_degrees(G['ulat'].data)
        ULON    = self.radians_to_degrees(G['ulon'].data)
        TLON_b, TLAT_b = self.build_grid_corners(TLAT, TLON)
        ULON_b, ULAT_b = self.build_grid_corners(ULAT, ULON)
        T_ANGLE = self.radians_to_degrees(G['angleT'].data)
        U_ANGLE = self.radians_to_degrees(G['angle'].data)
        TAREA   = G['tarea'].data
        UAREA   = G['uarea'].data
        j, i    = TLAT.shape
        jb, ib  = j + 1, i + 1
        nat_dim = self.CICE_dict["spatial_dims"]  # e.g., ("nj", "ni")
        ext_dim = tuple(f"{dim}_b" for dim in nat_dim)
        coords  = {nat_dim[0]: np.arange(j),
                   nat_dim[1]: np.arange(i),
                   ext_dim[0]: np.arange(jb),
                   ext_dim[1]: np.arange(ib)}
        G_t     = {'lat'     : (nat_dim, TLAT, {'units': 'degrees'}),
                   'lat_b'   : (ext_dim, TLAT_b, {'units': 'degrees'}),
                   'lon'     : (nat_dim, TLON, {'units': 'degrees'}),
                   'lon_b'   : (ext_dim, TLON_b, {'units': 'degrees'}),
                   'angle'   : (nat_dim, T_ANGLE, {'units': 'degrees'}),
                   'area'    : (nat_dim, TAREA, {'units': 'm^2'}),
                   'kmt_org' : (nat_dim, KMT_org, {'units'      : 'binary',
                                                   'description': '1=land, 0=ocean',
                                                   'long_name'  : 'original landmask on t-grid'}),
                   'kmt_mod' : (nat_dim, KMT_mod, {'units'      : 'binary',
                                                   'description': '1=land, 0=ocean',
                                                   'long_name'  : 'modified t-grid-landmask to simulate grounded icebergs'})}
        G_u     = {'lat'     : (nat_dim, ULAT, {'units': 'degrees'}),
                   'lat_b'   : (ext_dim, ULAT_b, {'units': 'degrees'}),
                   'lon'     : (nat_dim, ULON, {'units': 'degrees'}),
                   'lon_b'   : (ext_dim, ULON_b, {'units': 'degrees'}),
                   'angle'   : (nat_dim, U_ANGLE, {'units': 'degrees'}),
                   'area'    : (nat_dim, UAREA, {'units': 'm^2'}),
                   'kmt_org' : (nat_dim, KMT_org, {'units': 'binary',
                                                   'description': '1=land, 0=ocean',
                                                   'long_name': 'original landmask on t-grid'}),
                   'kmt_mod' : (nat_dim, KMT_mod, {'units': 'binary',
                                                   'description': '1=land, 0=ocean',
                                                   'long_name': 'modified t-grid-landmask to simulate grounded icebergs'})}
        self.G_t = xr.Dataset(data_vars=G_t, coords=coords)
        self.G_u = xr.Dataset(data_vars=G_u, coords=coords)
        if slice_hem:
            self.G_t = self.slice_hemisphere(self.G_t)
            self.G_u = self.slice_hemisphere(self.G_u)
        self.bgrid_loaded = True
        
    def compute_regular_grid_area(self, da):
        """
        Compute cell areas (km^2) for a regular lat/lon grid (1D coords).

        Assumes 1D, monotonically increasing latitude and longitude coordinates with
        uniform longitudinal spacing. Uses spherical geometry to integrate the area
        between latitude edges for each latitude band, broadcasting across longitudes.

        Parameters
        ----------
        da : xr.DataArray
            Any DataArray that carries 1D coordinates named ``'lat'`` and ``'lon'``.
            Values are not used; only the coordinates are read.

        Returns
        -------
        xr.DataArray
            2D area array with dims (``lat``, ``lon``) in **km²**, aligned to the
            provided coordinate vectors.

        Notes
        -----
        - Earth radius is fixed at 6,371,000 m.
        - For non-uniform lon spacing or curvilinear grids, use model grid areas
        (e.g., CICE `tarea`) rather than this helper.
        """
        R   = 6371000.0  # Earth radius in meters
        lat = np.deg2rad(da['lat'].values)
        lon = np.deg2rad(da['lon'].values)
        # Latitude edges
        lat_edges       = np.zeros(len(lat)+1)
        lat_edges[1:-1] = (lat[:-1] + lat[1:]) / 2
        lat_edges[0]    = lat[0] - (lat[1]-lat[0])/2
        lat_edges[-1]   = lat[-1] + (lat[-1]-lat[-2])/2
        # Longitude spacing
        dlon = lon[1] - lon[0]
        # 1D cell area per latitude
        dA_lat = (R**2) * dlon * (np.sin(lat_edges[1:]) - np.sin(lat_edges[:-1]))
        # Broadcast to 2D
        area_2d = np.tile(dA_lat[:, np.newaxis], (1, len(lon)))
        # Convert to km^2
        area_2d /= 1e6
        area_da = xr.DataArray(area_2d, dims=("lat","lon"), coords={"lat":da['lat'], "lon":da['lon']})
        return area_da