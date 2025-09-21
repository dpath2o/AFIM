import gc, re, dask
import xarray      as xr
import numpy       as np
import pandas      as pd
from pathlib       import Path

class SeaIceGridWork:

    def __init__(self, **kwargs):
        return

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

    def define_reG_weights(self):
        """
        Define and store an xESMF regridder to remap CICE B-grid (U-point) data to the T-grid.

        This method constructs two CICE grids—one at the U-point (B-grid) and one at the T-point (centered)—using 
        internally defined grid definitions. It then either reuses existing xESMF regridding weights or creates 
        them if they do not already exist.

        The resulting regridder is stored as `self.reG`, and a flag `self.reG_weights_defined` is set to True 
        upon successful creation or reuse of the regridder.

        Grid details:
        - The source grid (`G_u`) is built using U-point coordinates (`grid_type='u'`).
        - The target grid (`G_t`) is built using T-point coordinates (`grid_type='t'`) and includes a land-sea mask.
        - Both grids include corner information to support conservative or bilinear regridding.

        Regridding parameters:
        - Method: "bilinear"
        - Periodic: True (assumes global grid)
        - Degenerate cells: ignored
        - Extrapolation: nearest-neighbor (source to destination)
        - Weight reuse: enabled if existing file is found at path `self.CICE_dict["P_reG_u2t_weights"]`

        Logging:
        - Logs whether existing weights are reused or new weights are created.

        Returns
        -------
        None
            The regridder is stored as `self.reG` and is not returned explicitly.
        """
        G_u           = self.define_cice_grid( grid_type='u'             , build_grid_corners=True )
        G_t           = self.define_cice_grid( grid_type='t' , mask=True , build_grid_corners=True )
        F_weights     = self.CICE_dict["P_reG_u2t_weights"]
        weights_exist = os.path.exists(F_weights)
        self.logger.info(f"{'Reusing' if weights_exist else 'Creating'} regrid weights: {F_weights}")
        self.reG = xe.Regridder(G_u, G_t,
                                method            = "bilinear",
                                periodic          = True,
                                ignore_degenerate = True,
                                extrap_method     = "nearest_s2d",
                                reuse_weights     = weights_exist,
                                filename          = F_weights)
        self.reG_weights_defined = True

    def reG_bgrid_to_tgrid_xesmf(self, da, coord_names=None):
        """
        Regrid a single B-grid DataArray to the T-grid using pre-defined xESMF regridder.

        If coord_names are not provided then assumes coordinate names provided in JSON configuration file,
        which are like ["ULON","ULAT"]

        INPUTS:
        da : xr.DataArray; b-grid variable with coordinates coord_names

        OUTPUTS:
        xr.DataArray; re-gridded DataArray on the T-grid with dimensions (time, nj, ni).
        """
        coord_names = coord_names if coord_names is not None else self.CICE_dict["bcoord_names"]
        if not set(coord_names).issubset(set(da.coords)):
            self.logger.error(f"Cannot regrid: as {coord_names} not found in coordinates.")
            return None
        coord_map = {}
        for name in coord_names:
            if "LAT" in name.upper():
                coord_map[name] = "lat"
            elif "LON" in name.upper():
                coord_map[name] = "lon"
        if set(coord_map.values()) != {"lat", "lon"}:
            self.logger.error(f"Could not identify lat/lon from coord_names: {coord_names}")
            return None
        da_tmp = da.rename(coord_map)
        try:
            da_reG = self.reG(da_tmp)
        except Exception as e:
            self.logger.error(f"Regridding failed: {e}")
            return None
        return da_reG

    def simple_spatial_averaging_bgrid_to_tgrid(self, var):
        """
        Dask-safe 4-point unweighted average from B-grid to T-grid.

        Uses efficient array shifting and avoids costly concatenation over new dimensions.

        Parameters
        ----------
        var : xr.DataArray
            2D or 3D (time, nj, ni) array on B-grid.

        Returns
        -------
        xr.DataArray
            Averaged field on T-grid with shape (time, nj, ni).
        """
        x_dim = self.CICE_dict["x_dim"]
        y_dim = self.CICE_dict["y_dim"]
        x_len = self.CICE_dict["x_dim_length"]
        y_len = self.CICE_dict["y_dim_length"]
        self.logger.info(f"input shape to spatial averaging: {var.shape}")
        self.logger.info("  → Slicing corner points for averaging...")
        v00 = var.isel({y_dim: slice(None, -1), x_dim: slice(None, -1)})
        v01 = var.isel({y_dim: slice(None, -1), x_dim: slice(1, None)})
        v10 = var.isel({y_dim: slice(1, None), x_dim: slice(None, -1)})
        v11 = var.isel({y_dim: slice(1, None), x_dim: slice(1, None)})
        self.logger.info("  → Computing mean of four corners...")
        avg = (v00 + v01 + v10 + v11) / 4.0
        self.logger.info("  → Padding with NaNs to restore original grid size...")
        pad_y = max(y_len - avg.sizes.get(y_dim, 0), 0)
        pad_x = max(x_len - avg.sizes.get(x_dim, 0), 0)
        avg = avg.pad({y_dim: (0, pad_y), x_dim: (0, pad_x)}, constant_values=np.nan)
        self.logger.info("  → Applying cyclic wrap for last column...")
        if avg.sizes.get(x_dim, 0) > 1:
            avg[{x_dim: -1}] = avg.isel({x_dim: 0})
        # Force re-slicing to expected grid size to ensure consistency
        avg = avg.isel({y_dim: slice(0, y_len), x_dim: slice(0, x_len)})
        if "time" in var.coords:
            avg = avg.assign_coords(time=var["time"])
            self.logger.info("  → Time coordinate restored.")
        assert avg.sizes[y_dim] == y_len, f"{y_dim} mismatch: got {avg.sizes[y_dim]}, expected {y_len}"
        assert avg.sizes[x_dim] == x_len, f"{x_dim} mismatch: got {avg.sizes[x_dim]}, expected {x_len}"
        for dim in [y_dim, x_dim]:
            if dim in avg.indexes:
                avg = avg.drop_indexes(dim)
        return avg

    def pygmt_regrid(self, da, lon, lat, grid_res=None, region=None, search_radius="200k"):
        """
        Regrid a 2D data array using PyGMT's nearneighbor interpolation.

        This method applies PyGMT's `nearneighbor` algorithm to interpolate scattered
        data values (`da`) onto a regular grid based on specified longitude and latitude 
        arrays. The input is masked to ignore NaNs or non-finite values.

        Parameters
        ----------
        da : xarray.DataArray
            2D array of data values to interpolate (e.g., sea ice thickness).
        lon : xarray.DataArray or np.ndarray
            Longitude values corresponding to `da`, same shape.
        lat : xarray.DataArray or np.ndarray
            Latitude values corresponding to `da`, same shape.
        grid_res : str or float, optional
            Grid spacing for the output grid (e.g., "0.5", "10m"). Required by PyGMT.
        region : list or tuple, optional
            Bounding box for the output grid in the form [west, east, south, north].
        search_radius : str or float, default "200k"
            Search radius for PyGMT's nearneighbor (e.g., "100k" for 100 km).

        Returns
        -------
        gridded : xarray.DataArray
            Gridded output with interpolated values over the defined region.

        Notes
        -----
        - All non-finite values in `da` are excluded prior to interpolation.
        - PyGMT must be properly installed and configured with GMT for this to work.
        """
        import pygmt
        mask = np.isfinite(da)
        df   = pd.DataFrame({"longitude": lon.values[mask].ravel(), 
                             "latitude" : lat.values[mask].ravel(),
                             "z"        : da.values[mask].ravel()})
        return pygmt.nearneighbor(data          = df,
                                  spacing       = grid_res,
                                  region        = region,
                                  search_radius = search_radius)
        
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