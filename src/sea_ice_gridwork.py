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
    
    def _circular_mean_lon_deg(self, *lon_deg_arrays):
        """
        Circular mean for longitudes in degrees. Inputs must be broadcastable to same shape.
        Returns degrees in [0, 360).
        """
        lon = np.stack(lon_deg_arrays, axis=0)  # (n, ...)
        lon = self.normalise_longitudes(lon, to="0-360")
        lonr = np.deg2rad(lon)
        x = np.cos(lonr).mean(axis=0)
        y = np.sin(lonr).mean(axis=0)
        out = np.rad2deg(np.arctan2(y, x))
        return self.normalise_longitudes(out, to="0-360")

    def build_grid_faces(self, lon, lat, source_in_radians=False):
        """
        Build C-grid geometry from *cell-centre* lon/lat arrays.

        Parameters
        ----------
        lon, lat : array-like
            2D arrays of cell-centre lon/lat, either in degrees or radians.
            Shape (nj, ni).
        source_in_radians : bool
            If True, lon/lat are radians and will be converted to degrees.

        Returns
        -------
        lon_b, lat_b : ndarray
            Corner coordinates, shape (nj+1, ni+1)
        lon_e, lat_e : ndarray
            Vertical face centres (E/W faces), shape (nj, ni+1)
            (dims can be interpreted as (nj, ni_b))
        lon_n, lat_n : ndarray
            Horizontal face centres (N/S faces), shape (nj+1, ni)
            (dims can be interpreted as (nj_b, ni))
        """
        if source_in_radians:
            lon = self.radians_to_degrees(lon)
            lat = self.radians_to_degrees(lat)
        lon = self.normalise_longitudes(lon, to="0-360")
        lat = np.asarray(lat)

        nj, ni = lat.shape

        # --- corners (nj+1, ni+1)
        lon_b = np.full((nj + 1, ni + 1), np.nan, dtype=float)
        lat_b = np.full((nj + 1, ni + 1), np.nan, dtype=float)

        # Interior corners: mean of surrounding 4 cell centres
        lat_b[1:-1, 1:-1] = 0.25 * (lat[:-1, :-1] + lat[:-1, 1:] + lat[1:, :-1] + lat[1:, 1:])
        lon_b[1:-1, 1:-1] = self._circular_mean_lon_deg(lon[:-1, :-1], lon[:-1, 1:], lon[1:, :-1], lon[1:, 1:])

        # Edges: mean of adjacent 2 cell centres (more consistent than copying)
        # top edge (j=0): between cells (0, i-1) and (0, i)
        lat_b[0, 1:-1] = 0.5 * (lat[0, :-1] + lat[0, 1:])
        lon_b[0, 1:-1] = self._circular_mean_lon_deg(lon[0, :-1], lon[0, 1:])

        # bottom edge (j=nj): between cells (nj-1, i-1) and (nj-1, i)
        lat_b[-1, 1:-1] = 0.5 * (lat[-1, :-1] + lat[-1, 1:])
        lon_b[-1, 1:-1] = self._circular_mean_lon_deg(lon[-1, :-1], lon[-1, 1:])

        # left edge (i=0): between cells (j-1,0) and (j,0)
        lat_b[1:-1, 0] = 0.5 * (lat[:-1, 0] + lat[1:, 0])
        lon_b[1:-1, 0] = self._circular_mean_lon_deg(lon[:-1, 0], lon[1:, 0])

        # right edge (i=ni): between cells (j-1,ni-1) and (j,ni-1)
        lat_b[1:-1, -1] = 0.5 * (lat[:-1, -1] + lat[1:, -1])
        lon_b[1:-1, -1] = self._circular_mean_lon_deg(lon[:-1, -1], lon[1:, -1])

        # Four outer corners: nearest cell centre
        lat_b[0, 0] = lat[0, 0]
        lat_b[0, -1] = lat[0, -1]
        lat_b[-1, 0] = lat[-1, 0]
        lat_b[-1, -1] = lat[-1, -1]

        lon_b[0, 0] = lon[0, 0]
        lon_b[0, -1] = lon[0, -1]
        lon_b[-1, 0] = lon[-1, 0]
        lon_b[-1, -1] = lon[-1, -1]

        lon_b = self.normalise_longitudes(lon_b, to="0-360")

        # --- face centres from corners
        # Vertical faces (E/W): average corners above/below along a vertical edge
        # shape (nj, ni+1) = (nj, ni_b)
        lat_e = 0.5 * (lat_b[:-1, :] + lat_b[1:, :])
        lon_e = self._circular_mean_lon_deg(lon_b[:-1, :], lon_b[1:, :])

        # Horizontal faces (N/S): average corners left/right along a horizontal edge
        # shape (nj+1, ni) = (nj_b, ni)
        lat_n = 0.5 * (lat_b[:, :-1] + lat_b[:, 1:])
        lon_n = self._circular_mean_lon_deg(lon_b[:, :-1], lon_b[:, 1:])

        return lon_b, lat_b, lon_e, lat_e, lon_n, lat_n


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

    def load_cice_grid(self, slice_hem=False, build_faces=True):
        """
        Unified grid loader.

        Sets class attributes (no return):
        - self.G_t : t-grid dataset (centres + corners, mask, angle, area)
        - self.G_u : legacy u-grid dataset (centres + corners, mask, angle, area)
        - self.G_e : C-grid vertical-face coords derived from t-grid corners (lon_e/lat_e)
        - self.G_n : C-grid horizontal-face coords derived from t-grid corners (lon_n/lat_n)

        Parameters
        ----------
        slice_hem : bool
            If True, applies self.slice_hemisphere() to all constructed grid datasets.
        build_faces : bool
            If True, compute C-grid face coordinate datasets (G_e, G_n).
        """
        if self.grid_loaded:
            self.logger.info("CICE grid previously loaded ... returning")
            return
        nat_dim = self.CICE_dict["spatial_dims"]      # e.g. ("nj","ni")
        ext_dim = tuple(f"{d}_b" for d in nat_dim)    # e.g. ("nj_b","ni_b")
        with xr.open_dataset(self.CICE_dict["P_G"]) as G:
            # land masks
            kmt_org = xr.open_dataset(self.P_KMT_org).kmt.data
            kmt_mod = xr.open_dataset(self.P_KMT_mod).kmt.data if self.use_gi else kmt_org
            # --- T grid (centres)
            TLAT = self.radians_to_degrees(G["tlat"].data)
            TLON = self.normalise_longitudes(self.radians_to_degrees(G["tlon"].data), to="0-360")
            TAREA = G["tarea"].data
            T_ANGLE = self.radians_to_degrees(G["angleT"].data)
            # corners + faces from t-grid centres
            ULON_b, ULAT_b, lon_e, lat_e, lon_n, lat_n = self.build_grid_faces(TLON, TLAT, source_in_radians=False)
            # --- U grid (legacy / B-grid velocity points)
            ULAT = self.radians_to_degrees(G["ulat"].data)
            ULON = self.normalise_longitudes(self.radians_to_degrees(G["ulon"].data), to="0-360")
            UAREA = G["uarea"].data
            U_ANGLE = self.radians_to_degrees(G["angle"].data)
        # dims
        nj, ni = TLAT.shape
        coords_nat = {nat_dim[0]: np.arange(nj),
                      nat_dim[1]: np.arange(ni)}
        coords_ext = {nat_dim[0]: np.arange(nj),
                      nat_dim[1]: np.arange(ni),
                      ext_dim[0]: np.arange(nj + 1),
                      ext_dim[1]: np.arange(ni + 1)}
        # landmasks
        self.kmt_org = xr.Dataset(data_vars = {"kmt_org": (nat_dim, kmt_org, {"units"    : "binary", 
                                                                              "long_name": "original landmask on t-grid"})},
                                  coords    = coords_nat)
        # Grounded icebergs
        if self.use_gi:
            self.kmt_mod = xr.Dataset(data_vars = {"kmt_mod": (nat_dim, kmt_mod, {"units"    : "binary", 
                                                                                "long_name": "modified landmask on t-grid (grounded icebergs)"})},
                                    coords    = coords_nat)
            # Difference: grounded icebergs are cells that changed from ocean (1) to land (0)
            GI_mask   = (kmt_org == 1) & (kmt_mod == 0)
            GI_area   = TAREA[GI_mask]
            self.G_GI = xr.Dataset(data_vars = {'mask' : (nat_dim, GI_mask.astype('float32')),
                                                'area' : (nat_dim, GI_area),
                                                'lat'  : (nat_dim, TLAT[GI_mask]),
                                                'lon'  : (nat_dim, TLON[GI_mask])},
                                    coords   = coords_nat)
            # Get coordinates of affected cells (shifted west by one ni index to match B-grid layout)
            nj_idx, ni_idx = np.where(GI_mask)
            ni_idx_shifted = ni_idx - 1
            valid          = ni_idx_shifted >= 0
            nj_idx         = nj_idx[valid]
            ni_idx         = ni_idx_shifted[valid]
            # Save 1D arrays with iceberg IDs
            nGI              = np.arange(len(TLAT[nj_idx, ni_idx]))
            self.G_GI['lat_nGI'] = (('nGI',), TLAT[nj_idx, ni_idx])
            self.G_GI['lon_nGI'] = (('nGI',), TLON[nj_idx, ni_idx])
            self.G_GI['lat_nGI'].attrs.update({'long_name': 'Grounded iceberg latitude'})
            self.G_GI['lon_nGI'].attrs.update({'long_name': 'Grounded iceberg longitude'})
        # --- build datasets
        self.G_t = xr.Dataset(data_vars = {"lat"    : (nat_dim, TLAT, {"units": "degrees"}),
                                        "lon"    : (nat_dim, TLON, {"units": "degrees"}),
                                        "angle"  : (nat_dim, T_ANGLE, {"units": "degrees"}),
                                        "area"   : (nat_dim, TAREA, {"units": "m^2"})},
                            coords    = coords_nat)
        self.G_u = xr.Dataset(data_vars = {"lat"    : (nat_dim, ULAT, {"units": "degrees"}),
                                        "lon"    : (nat_dim, ULON, {"units": "degrees"}),
                                        "lat_b"  : (ext_dim, ULAT_b, {"units": "degrees"}),
                                        "lon_b"  : (ext_dim, ULON_b, {"units": "degrees"}),
                                        "angle"  : (nat_dim, U_ANGLE, {"units": "degrees"}),
                                        "area"   : (nat_dim, UAREA, {"units": "m^2"})},
                            coords    = coords_ext)
        if build_faces:
            # Use ext_dim/nat_dim naming so shapes are explicit:
            # lon_e/lat_e: (nj, ni_b) == (nat_dim[0], ext_dim[1])
            # lon_n/lat_n: (nj_b, ni) == (ext_dim[0], nat_dim[1])
            self.G_e = xr.Dataset(data_vars = {"lat": ((nat_dim[0], ext_dim[1]), lat_e, {"units": "degrees"}),
                                               "lon": ((nat_dim[0], ext_dim[1]), lon_e, {"units": "degrees"})},
                                  coords    = coords_ext)
            self.G_n = xr.Dataset(data_vars = {"lat": ((ext_dim[0], nat_dim[1]), lat_n, {"units": "degrees"}),
                                               "lon": ((ext_dim[0], nat_dim[1]), lon_n, {"units": "degrees"})},
                                  coords    = coords_ext)
        else:
            self.G_e = None
            self.G_n = None
        if slice_hem:
            self.G_t = self.slice_hemisphere(self.G_t)
            self.G_u = self.slice_hemisphere(self.G_u)
            if self.G_e is not None:
                self.G_e = self.slice_hemisphere(self.G_e)
            if self.G_n is not None:
                self.G_n = self.slice_hemisphere(self.G_n)
        # flags for downstream
        self.bgrid_loaded = True
        self.cgrid_loaded = bool(build_faces)
        self.grid_loaded  = True

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